#include "cfr.h"
#include "strategy_io.h"
#include <filesystem>
#include <pokerstove/peval/CardSet.h>
#include <pokerstove/peval/Card.h>
#include <future>

int CFR::Solve(Strategy *blueprint,
               Strategy *strategy,
               const sCFRProgress &convergence,
               const std::atomic_bool &cancelled,
               int starting_checkpoint)
{
    auto num_thread = cfr_param_.num_threads;
    logger::info("   [CFR]: solving begins with [%d threads] [round %d] [checkpoint %d] [depth limited %d]",
                 num_thread,
                 strategy->ag_->root_node_->GetRound(),
                 starting_checkpoint,
                 cfr_param_.depth_limited
    );

    //more safegaurding codes
    if (strategy->ag_->root_node_->GetRound() != HOLDEM_ROUND_RIVER) {
        if (strategy->ag_->root_hand_beliefs_for_all_[0].NonZeroBeliefCount() <= 30
            || strategy->ag_->root_hand_beliefs_for_all_[1].NonZeroBeliefCount() <= 30) {
            logger::warn(
                    "🚨the new ag is not at river and the range width is too small. It may all become 0 in cfr solving"
            );
        }
    }

    SimpleTimer timer;
    //prep
    auto local_commands = commands_;
    //copy local param
    auto local_cfr_param = cfr_param_;

    //allocate flops if running hierarchical bucket and solve from round 0
    auto pub_bucket_flop_boards = new std::vector<Board_t>[HIER_PUB_BUCKET];
    auto thread_board = new std::vector<int>[num_thread];
    if (strategy->ag_->root_node_->GetRound() == 0
        && (strategy->ag_->bucket_reader_.buckets_[1]->type_ == HIERARCHICAL_BUCKET
            || strategy->ag_->bucket_reader_.buckets_[2]->type_ == HIERARCHICAL_BUCKET
            || strategy->ag_->bucket_reader_.buckets_[3]->type_ == HIERARCHICAL_BUCKET)) {
        AllocateFlops(pub_bucket_flop_boards, thread_board, num_thread, strategy->ag_->game_.numBoardCards[1]);
    }

    //prepare for vector cfr
    if (local_cfr_param.cfr_mode_ == CFR_VECTOR_PAIRWISE_SOLVE ||
        local_cfr_param.cfr_mode_ == CFR_VECTOR_ALTERNATE_SOLVE) {
        logger::trace("allocate memory for adjustment");
        strategy->ag_->root_node_->InitAdjustmentFactor_r(strategy->ag_->kernel_->bmax_by_r_);
    }

    // NOTE(kwok): The Grand CFR loop.

    sCFRProgress thread_local_current_progress;

    /*
     * Run all commands by default,
     *      - unless the cfr_solving termininate it naturally.
     *      - unless it is cancalled from the outside
     */
    bool keep_solving = true;
    thread_local_current_progress.iteration = starting_checkpoint;
    while (!local_commands.empty() && keep_solving && !cancelled) {
#if DEV > 1
        SimpleTimer cmd_timer;
#endif
        auto cmd = local_commands.front();
        if (cmd.trigger_iter_ >= starting_checkpoint) {
            // only print this in blueprint training for now. Not using sgs in preflop.
            if (strategy->ag_->root_node_->GetRound() == HOLDEM_ROUND_PREFLOP) {
                cmd.Print();
            }
            switch (cmd.type_) {
                case CMD_SCALAR_SOLVE :
                case CMD_VECTOR_ALTERNATE_SOLVE :
                case CMD_VECTOR_PAIRWISE_SOLVE : {
                    // todo: in fact we are not using the convergence state at all.
                    logger::info("🚦current state iteration = %d, convergence iteration = %d",
                                 thread_local_current_progress.iteration, convergence.iteration
                    );
                    if (thread_local_current_progress < convergence) {
                        sTotalThreadOutput total_result;
                        ThreadedCfrSolve(blueprint,
                                         strategy,
                                         local_cfr_param,
                                         cmd.steps_,
                                         total_result,
                                         thread_board,
                                         pub_bucket_flop_boards,
                                         thread_pool_,
                                         cancelled
                        );
                        if (!cancelled) {
                            // only update the current state if not abruptly cancelled from outside
                            thread_local_current_progress.UpdateState(cmd.steps_, timer.GetLapseFromBegin(), total_result.avg_);
                        }
                    } else {
                        // break all loops
                        keep_solving = false;
                    }
                    break;
                }
                case CMD_VECTOR_PROFILING : {
                    auto profiler_cfr_param = sCfrParam();
                    sTotalThreadOutput br;
                    profiler_cfr_param = local_cfr_param;
                    profiler_cfr_param.ProfileModeOn();
                    if (!cfr_param_.depth_limited) {
                        ThreadedCfrSolve(blueprint,
                                         strategy,
                                         profiler_cfr_param,
                                         cmd.steps_,
                                         br,
                                         thread_board,
                                         pub_bucket_flop_boards,
                                         thread_pool_,
                                         cancelled
                        );
                        auto br_tuple = ExplTuple{cmd.trigger_iter_, 1, br.avg_, br.max_, br.min_,
                                                  br.avg_ + br.std_dev_, br.avg_ - br.std_dev_};
                        profiling_writer_.WriteToFile(br_tuple);
                    }
                    logger::debug("PROFILER: br[%f] at iter [%d]", br.avg_, thread_local_current_progress.iteration);
                    // In the case of blueprint training, inspect the strategy over time
                    if (strategy->ag_->root_node_->GetRound() == HOLDEM_ROUND_PREFLOP) {
                        auto name = strategy->ag_->name_ + "_" + local_cfr_param.name.substr(0, 7) + "_"
                                    + std::to_string(thread_local_current_progress.iteration);
                        strategy->InspectPreflopBets(name, profiler_cfr_param.profiling_strategy_);
                    }
                    break;
                }
                case CMD_VECTOR_PRUNING_ON : {
                    local_cfr_param.pruning_on = true;
                    local_cfr_param.rollout_prune_thres *= local_cfr_param.rm_floor; //make flooring scale in variant
                    break;
                }
                case CMD_AVG_UPDATE_ON : {
                    local_cfr_param.rm_avg_update = AVG_CLASSIC;
                    break;
                }
                case CMD_DISCOUNTING : {
                    int step = (thread_local_current_progress.iteration / local_cfr_param.rm_disc_interval) + 1;
                    //def of LCFR, (3/2, 2)
                    auto regret_factor =
                            (double) pow(step, local_cfr_param.lcfr_alpha) /
                            (double) (pow(step, local_cfr_param.lcfr_alpha) + 1);
                    strategy->DiscountStrategy(STRATEGY_REG, regret_factor);

                    auto wavg_factor = (double) step / (double) (step + 1);
                    if (local_cfr_param.rm_avg_update == AVG_CLASSIC)
                        strategy->DiscountStrategy(STRATEGY_WAVG, wavg_factor);
                    break;
                }
                case CMD_SAVE_REG_WAVG : {
                    // Assuming CFR name beginning with either "cfrv" or "cfrs"
                    auto name = strategy->ag_->name_ + "_" + local_cfr_param.name.substr(0, 7) + "_"
                                + std::to_string(thread_local_current_progress.iteration);
                    SaveAG(strategy, name);
                    local_cfr_param.SaveCFRConfig(name);
                    // For now, we don't need checkpoints to do continuous training, hence no need for regret
                    // checkpoints to merge strategy. Just mkake it simple here.
                    strategy->InitMemory(STRATEGY_ZIPAVG, CFR_SCALAR_SOLVE);
                    // FIXME(kwok): This command potentially triggers the "b_max %d must >= num_threads %d" fatal error.
                    strategy->ConvertWavgToZipAvg(thread_pool_, num_thread);
                    SaveStrategy(strategy, STRATEGY_ZIPAVG, name);
                    strategy->ClearZipAvgMemory();
                    // SaveRegret(strategy, name);
                    // if (local_cfr_param.rm_avg_update != AVG_OFF) {
                    //     SaveWeightedAvg(strategy, name);
                    // }
                    break;
                }
                case CMD_PRINT_ROOT_STG: {
                    STRATEGY_TYPE calc_mode = STRATEGY_WAVG;
                    logger::debug("🔬%s printing the inspection of the root subgame with %s",
                                  profiling_writer_.prefix_,
                                  StrategyToNameMap[calc_mode]
                    );
                    strategy->InspectNode(strategy->ag_->root_node_, profiling_writer_.prefix_, calc_mode);
                    break;
                }
                case CMD_DUMMY: {
                    break;
                }
                default: {
                    logger::critical("no default command! check command pipeline");
                }
            }
#if DEV > 1
            cmd.Log(thread_local_current_progress.iteration.has_value() ? thread_local_current_progress.iteration.value() : -1,
                    cmd_timer.GetLapseFromBegin(),
                    timer.GetLapseFromBegin());
#endif
        }
        // pop outside the checkpoint_check loop
        local_commands.pop_front();
    }

    delete[] pub_bucket_flop_boards;
    delete[] thread_board;

    int return_code;
    if (cancelled) {
        //return the iteration checkpoint;
        return_code = thread_local_current_progress.iteration >= CFR_SOLVING_TERMINATED_ASYNC_CHECKPOINT
                      ? thread_local_current_progress.iteration
                      : CFR_SOLVING_TERMINATED_ASYNC_CHECKPOINT;
    } else if (!keep_solving) {
        return_code = CFR_SOLVING_TERMINATED_EARLY;
    } else {
        return_code = CFR_SOLVING_TERMINATED_ALL_COMMANDS;
    }

    logger::info("   [CFR]: ⏳solving took %d (ms) || ends at iter %d || end code %s",
                 timer.GetLapseFromBegin(),
                 thread_local_current_progress.iteration,
                 PrintCfrResultCode(return_code));
    return return_code;
}

/**
 * Life Cycle:
 *
 * thread_prepare
 *      main_loop
 *          iter_prepare
 *              alternate_training_prepare
 *              alternate_training_clean_up
 *          iter_clean_up
 * thread_clean_Up
*/
void *CFR::CfrSolve(void *thread_args)
{
    std::ostringstream thread_id_oss;
    thread_id_oss << std::this_thread::get_id();

    std::string thread_id = thread_id_oss.str();

    auto *args = (sThreadInput *) thread_args;
    auto ag = args->strategy_->ag_;

    // If the current is round pre-flop, assume hierarchically bucketing is used by default
    std::vector<Board_t> my_flops;
    if (ag->root_node_->GetRound() == 0
        && (ag->bucket_reader_.buckets_[1]->type_ == HIERARCHICAL_BUCKET
            || ag->bucket_reader_.buckets_[2]->type_ == HIERARCHICAL_BUCKET
            || ag->bucket_reader_.buckets_[3]->type_ == HIERARCHICAL_BUCKET)) {
        // get the assigned boards
        for (int &pub_board_idx: args->thread_board_) {
            auto v = args->pub_bucket_flop_boards_[pub_board_idx];
            my_flops.insert(my_flops.end(), v.begin(), v.end());
        }

        logger::info("[🧵thread %d] has been assigned %d flops", args->thread_idx_, my_flops.size());

        // Shuffle the flop order, as the original order is according to public bucket
        auto rd = std::random_device{};
        auto rng = std::default_random_engine{rd()};
        std::shuffle(my_flops.begin(), my_flops.end(), rng);
    }

    auto remaining_iter = args->iterations_;
    args->output_->Prepare(remaining_iter);

    CfrWorker *worker;
    switch (args->cfr_param_.cfr_mode_) {
        case CFR_VECTOR_PAIRWISE_SOLVE:
        case CFR_VECTOR_ALTERNATE_SOLVE:
            worker = new VectorCfrWorker(args->blueprint_, args->strategy_, &args->cfr_param_, my_flops, args->seed_);
            worker->SetWalkingMode(args->cfr_param_.cfr_mode_);
            break;
        case CFR_SCALAR_SOLVE:
            worker = new ScalarCfrWorker(args->blueprint_, args->strategy_, &args->cfr_param_, my_flops, args->seed_);
            break;
        default:
            logger::critical("unsupported cfr type");
    }

    auto cur_flop_idx = 0;
    while (remaining_iter-- && !args->cancelled_token_) {
        Board_t board{};
        // NOTE(kwok): Here we're sampling the public chance events, i.e. the public cards. Required
        // private chance events will be sampled within ScalaraCfrWorker::Solve in the case of the
        // primitive Chance-Sampled variant, or be exhaustedly considered within VectorCfrWorker::Solve
        // in the case of the Public Chance Sampling variant.
        //
        // Note that we are not going to sample the public chance events one by one as the game proceeds.
        // Rather, we're sampling all required public chance events in one breath.
        SampleSequentialFullBoard(ag->root_state_, &ag->game_, board, cur_flop_idx, worker->my_flops_);
        // board.Print();
        double local_util = worker->Solve(board);
        args->output_->AddIterResult(local_util);
    }

    args->output_->remaining_iter = remaining_iter;
    args->output_->Process();

    // threads clean-up
    delete args;
    delete worker;
    return nullptr;
}

CFR::CFR(const char *config_file)
{
    std::filesystem::path dir(AUTODIDACT_DIR_CFG_ENG);
    std::filesystem::path filename(config_file);
    std::ifstream cfr_file(dir / filename);
    auto config_str = std::string(config_file);
    cfr_param_.name = config_str.substr(0, config_str.length() - 5);
    web::json::value cfr_config;
    if (cfr_file.good()) {
        std::stringstream buffer;
        buffer << cfr_file.rdbuf();
        cfr_config = web::json::value::parse(buffer);
        Config(cfr_config.at("cfr"));
        cfr_file.close();
        cfr_param_.raw_ = cfr_config;
        auto num_threads = cfr_param_.num_threads;
        thread_pool_ = new pthread_t[num_threads];
    } else {
        logger::error("unable to open cfr_config file %s", config_file);
    }
}

// NOTE(kwok): Only useful for hierarchical buckets.
void CFR::AllocateFlops(std::vector<Board_t> *pub_flop_boards,
                        std::vector<int> *thread_board,
                        int num_thread,
                        int num_flop_cards)
{
    //only thread by flop boards if training blueprint
    Bucket bucket;
    bucket.LoadHierarchicalPublic();

    //put raw boards into pub board
    pokerstove::combinations cards(52, num_flop_cards);
    int num_combinations = 0;
    do {
        Board_t board{};
        pokerstove::CardSet hand;
        for (auto i = 0; i < num_flop_cards; i++) {
            board.cards[i] = cards[i];
            hand.insert(pokerstove::Card(cards[i]));
        }
        auto public_bucket = bucket.GetPublicBucket(hand.canonize().colex());
        pub_flop_boards[public_bucket].push_back(board);
        num_combinations++;
    } while (cards.next());

    std::vector<std::pair<int, int>> pub_flop_board_sortable;
    pub_flop_board_sortable.reserve(HIER_PUB_BUCKET);
    for (int i = 0; i < HIER_PUB_BUCKET; i++) {
        pub_flop_board_sortable.emplace_back(i, pub_flop_boards[i].size());
    }

    //sort in decending order of number of flops in public bucket
    std::sort(pub_flop_board_sortable.begin(), pub_flop_board_sortable.end(), [](const auto &l, const auto &r) -> bool
              {
                  return l.second > r.second;
              }
    );

    if (num_thread >= 60 && num_thread % 30 != 0) {
        logger::critical("if thread >= 60, it has be be multiple of 30");
    }

    if (num_thread > 1) {
        int cluster = num_thread > 60 ? 30 : num_thread;
        //evenly distribute work
        int cap = num_combinations / cluster;
        int idx = 0;
        int tally[cluster];
        std::memset(tally, 0, sizeof(int) * cluster);
        for (auto &[pub, size]: pub_flop_board_sortable) {
            bool suboptimal_assignment_flag = true;
            for (int i = 0; i < cluster; i++) {
                int new_size = tally[i] + size;
                if (cap - new_size >= 0) {
                    thread_board[i].push_back(pub);
                    tally[i] = new_size;
                    auto min_idx = std::min_element(tally, tally + cluster);
                    idx = std::distance(tally, min_idx);
                    suboptimal_assignment_flag = false;
                    break;
                }
            }
            if (suboptimal_assignment_flag) {
                logger::debug("all threads have been assigned full amount, adding to thread %d (with least flops %d)",
                              idx,
                              tally[idx]
                );
                thread_board[idx].push_back(pub);
                tally[idx] += size;
                auto min_idx = std::min_element(tally, tally + cluster);
                idx = std::distance(tally, min_idx);
            }
        }
        // copy over cluster by cluster
        if (num_thread > 60) {
            int rep = num_thread / cluster;
            for (int c = 0; c < cluster; c++) {
                for (int r = 1; r < rep; r++) {
                    thread_board[c + r * cluster] = thread_board[c];
                }
            }
        }
    } else {
        for (int i = 0; i < HIER_PUB_BUCKET; i++) {
            thread_board[0].push_back(i);
        }
    }

    for (int t = 0; t < num_thread; t++) {
        logger::debug("thread %d was assigned %d public flops", t, thread_board[t].size());
    }
}

void CFR::ThreadedCfrSolve(Strategy *blueprint,
                           Strategy *strategy,
                           sCfrParam &cfr_param,
                           int steps,
                           sTotalThreadOutput &total_result,
                           std::vector<int> *thread_board,
                           std::vector<Board_t> *pub_bucket_flop_boards,
                           pthread_t *thread_pool,
                           const std::atomic_bool &cancelled)
{

    int num_threads = cfr_param.num_threads;
    int effective_thread = num_threads;

    // NOTE(kwok): it is typical to profile with just 1 thread
    if (steps < num_threads) {
        effective_thread = 1;
        if (cfr_param.cfu_compute_acting_playing == WEIGHTED_RESPONSE) {
            logger::warn("🚨running single thread");
        }
    }

    sThreadOutput thread_output[effective_thread];
    if (effective_thread > 1) {
        int iter_avg = steps / effective_thread;

        /* Deal with inputs to the threads and firing them */

        int last_distinct_quota = -1;
        // NOTE(kwok): fire threads
        for (auto i = 0; i < effective_thread; i++) {
            int quota = iter_avg;
            if (i == effective_thread - 1) {
                auto remain = steps - (iter_avg * (effective_thread - 1));
                quota = remain;
            }
            auto thread_input = new sThreadInput(blueprint,
                                                 strategy,
                                                 pub_bucket_flop_boards,
                                                 i,
                                                 cfr_param,
                                                 thread_board[i],
                                                 &thread_output[i],
                                                 quota,
                                                 cancelled,
                                                 std::random_device()());
            if (quota != last_distinct_quota) {
                logger::info("[🧵thread %d] assigned %d iterations", i, quota);
                logger::info("\t\tthreads being assigned with the same quota of %d were omitted here ...", quota);
                last_distinct_quota = quota;
            }
            if (pthread_create(&thread_pool[i], nullptr, CFR::CfrSolve, thread_input)) {
                logger::critical("failed to launch threads.");
            }
        }

        /* Deal with joining the threads and outputs from them */

        int last_distinct_unfinished_quota = -1;
        // NOTE(kwok): block the spawning thread while waiting for all the spawned threads to finish
        for (int i = 0; i < effective_thread; ++i) {
            if (pthread_join(thread_pool[i], nullptr)) {
                logger::error("Couldn't join to thread %d", i);
            }
            int unfinished_quota = thread_output[i].remaining_iter;
            if (unfinished_quota != last_distinct_unfinished_quota) {
                logger::info("[🧵thread %d] unfinished quota = %d", i, unfinished_quota);
                logger::info("\t\tthreads ending with the same unfinished quota of %d were omitted here ...",
                             unfinished_quota
                );
                last_distinct_unfinished_quota = unfinished_quota;
            }
        }
    } else {
        // NOTE(kwok): solo thread
        auto thread_input = new sThreadInput(blueprint,
                                             strategy,
                                             pub_bucket_flop_boards,
                                             0,
                                             cfr_param,
                                             thread_board[0],
                                             &thread_output[0],
                                             steps,
                                             cancelled,
                                             std::random_device()());
        CFR::CfrSolve(thread_input);
    }

    // Also required for thread = 1, so that the same variables can be accessed similarly
    total_result.MergeThreadOutputs(thread_output, effective_thread);
}

int CFR::AsyncCfrSolving(CFR *cfr,
                         Strategy *new_strategy,
                         Strategy *blueprint,
                         sCFRProgress *convergence,
                         const std::atomic_bool &cancelled,
                         int cfr_checkpoint)
{
    logger::info("🟢starting async MCCFR solving");
    return cfr->Solve(blueprint, new_strategy, *convergence, cancelled, cfr_checkpoint);
}

void CFR::Config(web::json::value data)
{
    bool check_must_opts = data.has_field("num_threads")
                           && data.has_field("algo")
                           && data.has_field("regret_matching");
    if (!check_must_opts) {
        logger::critical("please fill in the must-have for cfr.");
    }

    //check and set meta cfr params
    int threads = data.at("num_threads").as_integer();
    cfr_param_.num_threads = threads == -1 ? std::thread::hardware_concurrency() : threads;

    cfr_param_.cfr_mode_ = CfrModeMap[data.at("algo").as_string()];
    if (cfr_param_.cfr_mode_ == CFR_UNKNOWN) {
        logger::critical("unknown cfr algo");
    }

    if (data.has_field("rollout")) {
        auto rollout = data.at("rollout");
        if (rollout.has_field("pruning")) {
            //check and set overall cfr params
            auto rollout_prune = rollout.at("pruning");
            if (rollout_prune.has_field("regret_thres")) {
                cfr_param_.rollout_prune_thres = rollout_prune.at("regret_thres").as_double();
                //never prune at 0
                if (cfr_param_.rollout_prune_thres == 0.0) {
                    cfr_param_.rollout_prune_thres = -1.0;
                }
            }
            if (rollout_prune.has_field("prob")) {
                cfr_param_.rollout_prune_prob = rollout_prune.at("prob").as_double();
            }
        }
    }

    /*
     * rollin
     */
    if (data.has_field("rollin")) {
        auto rollin_estimator = data.at("rollin").at("estimator");
        if (rollin_estimator.has_field("my")) {
            cfr_param_.cfu_compute_acting_playing = RollinEstimatorMap[rollin_estimator.at("my").as_string()];
        }
        if (rollin_estimator.has_field("opp")) {
            cfr_param_.cfu_compute_opponent = RollinEstimatorMap[rollin_estimator.at("opp").as_string()];
        }
    }

    /*
     * regret matching
     */
    auto regret_matching = data.at("regret_matching");
    if (regret_matching.has_field("floor")) {
        cfr_param_.rm_floor = regret_matching.at("floor").as_double();
    }
    if (regret_matching.has_field("side_walk")) {
        cfr_param_.avg_side_update_ = regret_matching.at("side_walk").as_bool();
    }
    if (regret_matching.has_field("discounting")) {
        auto rm_discounting = regret_matching.at("discounting");
        if (!rm_discounting.has_field("interval")) {
            logger::critical("discounting must have interval value specified.");
        }
        cfr_param_.rm_disc_interval = rm_discounting.at("interval").as_integer();
        if (rm_discounting.has_field("alpha")) {
            cfr_param_.lcfr_alpha = rm_discounting.at("alpha").as_double();
        }
    }

    cfr_param_.rm_floor *= 100 * REGRET_SCALER; //todo hack
    //  logger::debug("cfr solving with regret matching floor %f", cfr_param_.rm_floor);

    if (data.has_field("convergence")) {
        auto conv_config = data.at("convergence");
        cfr_param_.iteration =
                conv_config.has_field("max_iter") ? conv_config.at("max_iter").as_integer() : cfr_param_.iteration;
        cfr_param_.timeout_ms =
                conv_config.has_field("timeout_ms") ? conv_config.at("timeout_ms").as_integer() : cfr_param_.timeout_ms;
        if (cfr_param_.timeout_ms.has_value() && cfr_param_.timeout_ms.value() <= 1000) {
            logger::critical("[CFR %s]: timeout_ms = %llu is way too small",
                             cfr_param_.name,
                             cfr_param_.timeout_ms.value());
        }
    }

    // depth limited searching
    if (data.has_field("depth_limited")) {
        auto dls_conf = data.at("depth_limited");
        cfr_param_.depth_limited = dls_conf.at("on").as_bool();
        //only cfrs and sidewalk
        if (cfr_param_.depth_limited) {
            if (cfr_param_.cfr_mode_ != CFR_SCALAR_SOLVE) {
                logger::critical("depth limit solving only supports scalar mode");
            }
            // depth limit does not support sidewalk avg update mode. it would be too slow
            cfr_param_.avg_side_update_ = false;
            if (dls_conf.has_field("reps")) {
                cfr_param_.depth_limited_rollout_reps_ = dls_conf.at("reps").as_integer();
            }
            if (dls_conf.has_field("cache")) {
                cfr_param_.depth_limited_cache_ = dls_conf.at("cache").as_bool();
            }
            if (cfr_param_.depth_limited) {
                logger::debug("depth limit cfr [%d][reps %d]",
                              cfr_param_.depth_limited,
                              cfr_param_.depth_limited_rollout_reps_
                );
            }
        }
    }
}

void CFR::BuildCMDPipeline()
{
    auto cfr_config = cfr_param_.raw_.at("cfr");

    std::deque<CfrCommandWrapper> initial_commands{};

    //discounting
    auto regret_matching = cfr_config.at("regret_matching");
    if (regret_matching.has_field("discounting")) {
        auto rm_discounting = regret_matching.at("discounting");
        auto first_iter = rm_discounting.at("first_iter").as_integer();
        auto steps = rm_discounting.at("iterations").as_integer();
        //do discount with interval in [first_iter, last_iter)
        for (int i = 0; i < steps; i++) {
            if (first_iter + i * cfr_param_.rm_disc_interval < cfr_param_.iteration) {
                initial_commands.emplace_back(first_iter + i * cfr_param_.rm_disc_interval, CMD_DISCOUNTING);
            }
        }
    }

    //set avg
    if (regret_matching.has_field("avg_update_on")) {
        auto start_iter = regret_matching.at("avg_update_on").as_integer();
        if (start_iter != -1) {
            initial_commands.emplace_back(start_iter, CMD_AVG_UPDATE_ON);
        }
    }

    if (cfr_config.has_field("checkpoint")) {
        auto cp_config = cfr_config.at("checkpoint");
        /*
         * profiler
         */
        if (cp_config.has_field("profiler")) {
            auto profiler_config = cp_config.at("profiler");
            auto interval = profiler_config.at("interval").as_integer();
            auto profiler_step = profiler_config.at("steps").as_integer();
            if (profiler_config.has_field("strategy")) {
                cfr_param_.profiling_strategy_ = StrategyMap[profiler_config.at("strategy").as_string()];
                if (cfr_param_.profiling_strategy_ == UNKNOWN) {
                    logger::critical("cfr config error || unknown profiling strategy");
                }
            }

            CFR_COMMAND exec_cmd = CMD_VECTOR_PROFILING;

            //starting from where we trun on wavg, if the profiler strategy is wavg
            int offset = 0;
            if (cfr_param_.profiling_strategy_ == STRATEGY_WAVG) {
                if (regret_matching.has_field("avg_update_on")) {
                    offset = regret_matching.at("avg_update_on").as_integer();
                }
            }

            for (int i = offset + interval; i <= cfr_param_.iteration; i += interval) {
                initial_commands.emplace_back(i, exec_cmd, profiler_step);
            }
        }
        if (cp_config.has_field("strategy_checkpoint")) {
            auto interval = cp_config.at("strategy_checkpoint").as_integer();
            for (int i = interval; i < cfr_param_.iteration; i += interval) {
                initial_commands.emplace_back(i, CMD_SAVE_REG_WAVG);
            }
        }
        if (cp_config.has_field("save_final")) {
            initial_commands.emplace_back(cfr_param_.iteration, CMD_SAVE_REG_WAVG);
        }
        if (cp_config.has_field("save_root_prob")) {
            initial_commands.emplace_back(cfr_param_.iteration, CMD_PRINT_ROOT_STG);
        }
    }

    //do pruning (first_iter, inf]
    if (cfr_config.has_field("rollout")) {
        auto rollout = cfr_config.at("rollout");
        if (rollout.has_field("pruning")) {
            auto rollout_prune = rollout.at("pruning");
            if (rollout_prune.has_field("first_iter")) {
                auto pruning = rollout_prune.at("first_iter").as_integer();
                if (pruning < cfr_param_.iteration) {
                    initial_commands.emplace_back(pruning, CMD_VECTOR_PRUNING_ON);
                }
            }
        }
    }

    std::sort(initial_commands.begin(), initial_commands.end());

    // insert cfr solve
    CFR_COMMAND cfr_command;
    switch (cfr_param_.cfr_mode_) {
        case CFR_VECTOR_PAIRWISE_SOLVE:
            cfr_command = CMD_VECTOR_PAIRWISE_SOLVE;
            break;
        case CFR_VECTOR_ALTERNATE_SOLVE:
            cfr_command = CMD_VECTOR_ALTERNATE_SOLVE;
            break;
        case CFR_SCALAR_SOLVE:
            cfr_command = CMD_SCALAR_SOLVE;
            break;
        default:
            logger::critical("cfr config error || unknown algo");
            break;
    }

    initial_commands.push_front(CfrCommandWrapper(0, CMD_DUMMY));
    while (!initial_commands.empty()) {
        auto first = initial_commands.front();
        initial_commands.pop_front();
        if (!initial_commands.empty()) {
            auto second = initial_commands.front();
            auto diff = second.trigger_iter_ - first.trigger_iter_;
            if (diff > 0) {
                commands_.emplace_back(first.trigger_iter_, cfr_command, diff);
                commands_.push_back(second);
            } else {
                commands_.push_back(second);
            }
        }
    }

    if (commands_.empty()) {
        //only add solve
        commands_.emplace_back(0, cfr_command, cfr_param_.iteration);
    }

    auto last_cmd = commands_.end() - 1;
    // todo: hard coded condition, assuming only cfr and profiling has step!=0
    auto end_steps = last_cmd->type_ == CMD_VECTOR_PROFILING ? 0 : last_cmd->steps_;
    auto ending_iter = last_cmd->trigger_iter_ + end_steps;
    if (ending_iter != cfr_param_.iteration) {
        // last cfr param
        commands_.emplace_back(ending_iter, cfr_command, cfr_param_.iteration - ending_iter);
    }

#if DEV > 1
    int count = 0;
    for (auto it = commands_.begin(); it != commands_.end() - 1; it++) {
        count++;
        it->Print();
        auto next = it + 1;
        auto step = it->type_ == CMD_VECTOR_PROFILING ? 0 : it->steps_;
        bool check = (it->trigger_iter_ + step == next->trigger_iter_);
        if (!check) {
            logger::critical("command pipeline error");
        }
    }
    count++;
    (commands_.end() - 1)->Print();
    if (count != (int) commands_.size()) {
        logger::critical("command pipeline failed size check");
    }
#endif
}

std::string PrintCfrResultCode(int code)
{
    std::string msg;
    if (code >= CFR_SOLVING_TERMINATED_ASYNC_CHECKPOINT) {
        code = CFR_SOLVING_TERMINATED_ASYNC;
    }
    return CFR_RESULT_MAP[code];
}
