#include "cfr_worker.h"

double ScalarCfrWorker::Solve(Board_t board, bool calc_bru_explo, double *out_bru_explo)
{
    auto ag = strategy->ag_;
    auto active_players = AbstractGame::GetActivePlayerNum();
    auto hands_info = sPrivateHandsInfo(active_players, board, gen);

    // NOTE(kwok): assemble local root private hand beliefs based on the given board
    // FIXME(kwok): The number of players is not supposed to be fixed to 2.
    std::array<sPrivateHandBelief *, 2> root_beliefs{};
    for (int p = 0; p < active_players; p++) {
        root_beliefs[p] = new sPrivateHandBelief(&ag->root_hand_beliefs_for_all_[p]);
        root_beliefs[p]->NormalizeExcludeBoard(board);
    }

    double avg_cfu = 0.0;

    // FIXME(kwok): What about blueprint computing?
    // TODO(kwok): See `n_priv_hand_samples` also as a hyper-parameter.
    static const int n_priv_hand_samples = 1;

    double *brus = nullptr;
    size_t brus_size = n_priv_hand_samples * 2;
    if (calc_bru_explo) {
        // FIXME(kwok): The number of players is not supposed to be fixed to 2.
        brus = new double[brus_size];
        for (int i = 0; i < n_priv_hand_samples; i++) {
            brus[i] = std::numeric_limits<double>::infinity();
        }
    }

    for (int i = 0; i < n_priv_hand_samples; i++) {
        // NOTE(kwok): On each iteration, we start by sampling all of chance’s actions: the public chance
        // events visible to each player, as well as the private chance events that are visible to only a
        // subset of the players. In poker, this corresponds to randomly choosing the public cards revealed
        // to the players, and the private cards that each player is dealt.
        //
        // NOTE(kwok): Here we're sampling all the private chance events, i.e. the private hands of each player.
        // The public chance event, i.e. the public cards, has already been sampled by the outside invoker of
        // ScalarCfrWorker::Solve.
        hands_info.SamplePrivateHandsForAll(ag, root_beliefs);
        if (cfr_param_->pruning_on && cfr_param_->rm_floor < 0) {
            iter_prune_flag = GenRndNumber(1, 100) <= cfr_param_->rollout_prune_prob * 100;
        }
        // NOTE(kwok): Next, we recursively traverse the portion of the game tree that is reachable given
        // the sampled chance events, and explore all the players’ actions.
        //
        // * ⬇ On the way from the root to the leaves, we pass forward two SCALAR values: the probability that
        // each player would take actions to reach their respective information sets, given their current
        // strategy and their private information.
        //
        // * ⬆ On the way back from the leaves to the root, we return a single SCALAR value: the sampled
        // counterfactual value v ̃i(σ, I) for player i. At each choice node for player i, these values are all
        // that is needed to calculate the regret for each action and update the strategy.
        //
        // Note that at a terminal node z ∈ Z, it takes O(1) work to determine the avg_cfu for player i, u_i(z).

        // walk down the training tree alternatively.
        for (int trainee = 0; trainee < active_players; trainee++) {
            avg_cfu += WalkTree(ag->root_node_,
                                trainee, hands_info, strategy,
                                std::optional<CFU_COMPUTE_MODE>(),
                                std::optional<STRATEGY_TYPE>(),
                                true
            );
        }
        if (calc_bru_explo) {
            for (int trainee = 0; trainee < active_players; trainee++) {
                // FIXME(kwok): The number of players is not supposed to be fixed to 2.
                WalkTree(ag->root_node_, trainee, hands_info, strategy,
                         BEST_RESPONSE, std::optional<STRATEGY_TYPE>(), true
                );
            }
            for (int trainee = 0; trainee < active_players; trainee++) {
                // FIXME(kwok): The number of players is not supposed to be fixed to 2.
                double bru = WalkTree(
                        ag->root_node_, trainee, hands_info, strategy,
                        std::optional<CFU_COMPUTE_MODE>(), STRATEGY_BEST_RESPONSE, false
                );
                brus[i * active_players + trainee] = bru;
            }
        }
    }

    // conduct a sidewalk for updating WAVG
    if (cfr_param_->avg_side_update_ && cfr_param_->rm_avg_update == AVG_CLASSIC) {
        int n_sidewalk_iters = 50;
        for (int i = 0; i < n_sidewalk_iters; i++) {
            hands_info.SamplePrivateHandsForAll(ag, root_beliefs);
            for (int trainee = 0; trainee < active_players; trainee++) {
                WavgUpdateSideWalk(ag->root_node_, trainee, hands_info, strategy);
            }
        }
    }

    // FIXME(kwok): The number of players is not supposed to be fixed to 2.
    avg_cfu /= (n_priv_hand_samples * 2);

    for (auto b: root_beliefs) {
        delete b;
    }

    if (calc_bru_explo) {
        assert(brus != nullptr);
        double bru_per_sample = 0.0;
        for (int i = 0; i < n_priv_hand_samples; i++) {
            double bru_per_trainee = 0.0;
            for (int trainee = 0; trainee < active_players; trainee++) {
                bru_per_trainee += brus[i * active_players + trainee];
            }
            bru_per_sample += bru_per_trainee / active_players;
        }
        *out_bru_explo = bru_per_sample / n_priv_hand_samples;
        delete[] brus;
        fprintf(stderr, "avg_cfu = %g, bru_explo = (%g + %g)/2 = %g\n", avg_cfu, brus[0], brus[1], *out_bru_explo);
    }

    return avg_cfu;
}

double ScalarCfrWorker::WalkTree(Node *this_node, int trainee, sPrivateHandsInfo &hands_info, Strategy *target_strategy,
                                 std::optional<CFU_COMPUTE_MODE> trainee_cfu_compute_mode_hint,
                                 std::optional<STRATEGY_TYPE> trainee_strategy_type_hint, bool learn)
{
    if (this_node->IsTerminal()) {
        return EvalTermNode(this_node, trainee, hands_info);
    }
    // depth limited solving, till the local root of next round (rather than the local root of this round)
    if (this_node->IsLeafNode()) {
        // SimpleTimer timer;
        auto result = EvalLeafRootNode(this_node, trainee, hands_info);
        // timer.Checkpoint("🦃all " + std::to_string(cfr_param_->depth_limited_rollout_reps_) + " rollouts for round " +
        //                  std::to_string(this_node->GetRound() - 1));
        return result;
    }
    return EvalChoiceNode(
            this_node, trainee, hands_info, target_strategy,
            trainee_cfu_compute_mode_hint, trainee_strategy_type_hint, learn
    );
}

double ScalarCfrWorker::EvalTermNode(Node *this_node, int trainee, sPrivateHandsInfo &hands_info)
{
    if (this_node->IsShowdown()) {
        int stake = this_node->GetStake(trainee);
        int payoff = hands_info.payoff_[trainee]; // payoff ⊃ {-1, 0, 1}, where 0 denotes tie, 1 win, and, -1 lose
        stake *= payoff;
        return (double) stake;
    } else {
        // fold
        int stake = this_node->GetStake(trainee);
        return (double) stake;
    }
}

/*
 * v0: lazy rollout for traverser, 3 reps avg, both using non-biased strategy
 * v1: lazy rollout, 3 reps, alternating, external sampling, avg cfu.
 * v2: profile it before you cache anything. cache cfu for the leaf node (by node ID and hands), and optimize leaf node action caching
 * v3: cache preflop as file
 * v4: ...
 */
double ScalarCfrWorker::EvalLeafRootNode(Node *leaf_root_node, int trainee, sPrivateHandsInfo &hands_info)
{
    // should be the last round. it's the pluribus way
    auto r = leaf_root_node->GetRound() - 1;

    // FIXME(kwok): The number of players is not supposed to be fixed to 2.
    Bucket_t b0 = 0;
    Bucket_t b1 = 0;
    if (cfr_param_->depth_limited_cache_) {
        // read from the cache if exists: (PREFLOP) RAM cache from file, (FLOP) rollout/NN cache if needed
        leaf_root_node->CreateValueCacheIfNull();
        if (r > HOLDEM_ROUND_FLOP) {
            logger::critical("We don't do DLS for post-flop");
        }
        // FIXME(kwok): The number of players is not supposed to be fixed to 2.
        b0 = hands_info.buckets_[0][r];
        b1 = hands_info.buckets_[1][r];
        if (leaf_root_node->value_cache_->Exists(b0, b1)) {
            auto player0_cfu = leaf_root_node->value_cache_->GetValue(b0, b1);
            // FIXME(kwok): The number of players is not supposed to be fixed to 2.
            if (trainee == 1) {
                player0_cfu *= -1;
            }
            return player0_cfu;
        }
    }

    /* rollout for the current trainee rather than pairwise or alternatively */
    double player0_rollout_cfu = RolloutLeafRootNode(leaf_root_node, hands_info);

    /* or use NN */

    // insert to cache
    if (cfr_param_->depth_limited_cache_) {
        leaf_root_node->value_cache_->SetValue(b0, b1, player0_rollout_cfu);
    }

    // FIXME(kwok): The number of players is not supposed to be fixed to 2.
    // only in terms of player 0
    return trainee == 0 ? player0_rollout_cfu : -player0_rollout_cfu;
}

void ScalarCfrWorker::ComputeCfu(Node *this_node, const double *children_cfus, double &out_this_node_cfu,
                                 CFU_COMPUTE_MODE cfu_compute_mode, const float *distr_rnb,
                                 const bool *prune_flag) const
{
    auto a_max = this_node->GetAmax();
    double expected_utility = 0.0;
    switch (cfu_compute_mode) {
        case WEIGHTED_RESPONSE: {
            for (auto a = 0; a < a_max; a++) {
                if (prune_flag[a]) {
                    // pruned node won't be taken into account when calculating the final CFU of the current node.
                    continue;
                }
                expected_utility += distr_rnb[a] * children_cfus[a];
            }
            break;
        }
        case SUM_RESPONSE: {
            for (int a = 0; a < a_max; a++) {
                if (prune_flag[a]) {
                    continue;
                }
                expected_utility += children_cfus[a];
            }
            break;
        }
        case BEST_RESPONSE: {
            // the best response is a strategy for a player that is optimal against the opponent strategy
            expected_utility = -std::numeric_limits<double>::infinity();
            for (int a = 0; a < a_max; a++) {
                if (distr_rnb[a] > 0) {
                    // pruning will never happen in the BEST_RESPONSE mode.
                    auto utility = children_cfus[a];
                    if (utility > expected_utility) {
                        expected_utility = utility;
                    }
                }
            }
            break;
        }
    }
    out_this_node_cfu = expected_utility;
}

void ScalarCfrWorker::CollectChildBRUs(Node *this_node, sPrivateHandsInfo &hands_info, const double *children_brus,
                                       const double &this_node_bru, Strategy *target_strategy, const bool *prune_flag)
{
    int acting_player = this_node->GetActingPlayer();
    auto a_max = this_node->GetAmax();
    auto r = this_node->GetRound();
    auto n = this_node->GetN();
    auto b = hands_info.buckets_[acting_player][r];
    auto rnb0 = target_strategy->ag_->kernel_->hash_rnba(r, n, b, 0);
    for (int a = 0; a < a_max; a++) {
        if (prune_flag[a]) {
            continue;
        }
        double bru_child_a = children_brus[a];
        target_strategy->double_bru->insert(rnb0 + a);
        target_strategy->double_bru->update_fn(
                rnb0 + a,
                [&](auto &bru)
                {
                    bru = bru_child_a;
                }
        );
    }
}

void ScalarCfrWorker::CollectRegrets(Node *this_node, sPrivateHandsInfo &hands_info, const double *children_cfus,
                                     const double &this_node_cfu, Strategy *target_strategy, const bool *prune_flag)
{
    int acting_player = this_node->GetActingPlayer();
    auto a_max = this_node->GetAmax();
    auto r = this_node->GetRound();
    auto n = this_node->GetN();
    auto b = hands_info.buckets_[acting_player][r];
    auto rnb0 = target_strategy->ag_->kernel_->hash_rnba(r, n, b, 0);
    for (auto a = 0; a < a_max; a++) {
        if (prune_flag[a]) {
            continue;
        }
        int diff = (int) round(children_cfus[a] - this_node_cfu);
#ifdef DEBUG_EAGER_LOOKUP
        // TODO(kwok): 🦊
            double temp_reg = target_strategy->eager_int_regret_[rnb0 + a] + diff; // NOTE(kwok): Accumulate regret values
            // clamp it
            int new_reg = (int) std::fmax(temp_reg, cfr_param_->rm_floor);
            // total regret should have a ceiling
            if (new_reg > 2107483647) { // 2147483647 - 4 * 10^7
                logger::critical(
                        "if the regret overflowed, think about the possibility of the regret %d not being converging. [temp_reg %f][diff %f][cfu_a %f][this_node_cfu %f]",
                        new_reg,
                        temp_reg,
                        diff,
                        children_cfus[a],
                        this_node_cfu
                );
            }
            target_strategy->eager_int_regret_[rnb0 + a] = new_reg;
#endif
        // TODO(kwok): 🐦
        target_strategy->int_regret_->insert(rnb0 + a);
        target_strategy->int_regret_->update_fn(
                rnb0 + a,
                [&](auto &regret)
                {
                    double temp_reg = regret + diff; // accumulate regret values
                    int new_reg = (int) std::fmax(temp_reg, cfr_param_->rm_floor); // clamp it
                    // the total regret should have a ceiling
                    if (new_reg > 2107483647) { // 2147483647 - 4 * 10^7
                        logger::critical(
                                "if the regret overflowed, think about the possibility of the regret %d not being converging. [temp_reg %f][diff %f][cfu_a %f][this_node_cfu %f]",
                                new_reg,
                                temp_reg,
                                diff,
                                children_cfus[a],
                                this_node_cfu
                        );
                    }
                    regret = new_reg;
                }
        );
    }
}

double
ScalarCfrWorker::EvalChoiceNode(Node *this_node, int trainee, sPrivateHandsInfo &hands_info, Strategy *target_strategy,
                                std::optional<CFU_COMPUTE_MODE> trainee_cfu_compute_mode_hint,
                                std::optional<STRATEGY_TYPE> trainee_strategy_type_hint,
                                bool learn)
{
    int acting_player = this_node->GetActingPlayer();
    bool is_trainee_turn = acting_player == trainee;

    CFU_COMPUTE_MODE resolved_cfu_compute_mode;
    if (is_trainee_turn && trainee_cfu_compute_mode_hint.has_value()) {
        resolved_cfu_compute_mode = trainee_cfu_compute_mode_hint.value();
    } else if (is_trainee_turn) {
        resolved_cfu_compute_mode = cfr_param_->cfu_compute_acting_playing;
    } else {
        resolved_cfu_compute_mode = cfr_param_->cfu_compute_opponent;
    }

    STRATEGY_TYPE resolved_strategy_type;
    if (is_trainee_turn && trainee_strategy_type_hint.has_value()) {
        resolved_strategy_type = trainee_strategy_type_hint.value();
    } else {
        resolved_strategy_type = cfr_param_->strategy_cal_mode_;
    }

    auto a_max = this_node->GetAmax();
    auto r = this_node->GetRound();
    auto n = this_node->GetN();
    auto b = hands_info.buckets_[acting_player][r];
    auto rnb0 = target_strategy->ag_->kernel_->hash_rnba(r, n, b, 0);

    if (is_trainee_turn) {
        double children_cfus[a_max];
        bool prune_flag[a_max];
        for (auto a = 0; a < a_max; a++) {
            prune_flag[a] = false;
            auto next_node = this_node->children[a];
            // do pruning if the flag is set. skip river nodes and nodes leading to terminal.
            if (iter_prune_flag && !next_node->IsTerminal() && next_node->GetRound() != HOLDEM_ROUND_RIVER) {
                // if (iter_prune_flag && !next_node->IsTerminal()) {
#ifdef DEBUG_EAGER_LOOKUP
                // TODO(kwok): 🦊
                if (target_strategy->eager_int_regret_[rnb0 + a] <= cfr_param_->rollout_prune_thres) {
                    prune_flag[a] = true;
                    continue;
                }
#endif
                // TODO(kwok): 🐦
                // TODO(kwok): Test if this operation blocks.
                target_strategy->int_regret_->insert(rnb0 + a);
                target_strategy->int_regret_->find_fn(
                        rnb0 + a,
                        [&](const auto &regret)
                        {
                            if (regret <= cfr_param_->rollout_prune_thres) {
                                prune_flag[a] = true;
                            }
                        }
                );
                if (prune_flag[a]) {
                    continue;
                }
            }
            children_cfus[a] = WalkTree(
                    next_node, trainee, hands_info, target_strategy,
                    trainee_cfu_compute_mode_hint, trainee_strategy_type_hint, learn
            );
        }

        // only WEIGHTED_RESPONSE is supported
        // if (cfr_param_->cfu_compute_acting_playing != WEIGHTED_RESPONSE) {
        //     logger::critical("scalar does not support BEST_RESPONSE eval");
        // }

        float distr_rnb[a_max];
        target_strategy->ComputeStrategy(this_node, b, distr_rnb, resolved_strategy_type);

        double this_node_cfu = 0.0;
        ComputeCfu(this_node, children_cfus, this_node_cfu, resolved_cfu_compute_mode, distr_rnb, prune_flag);

        if (is_trainee_turn && learn) {
            switch (resolved_cfu_compute_mode) {
                case WEIGHTED_RESPONSE:
                case SUM_RESPONSE:
                    CollectRegrets(this_node, hands_info, children_cfus, this_node_cfu, target_strategy, prune_flag);
                    break;
                case BEST_RESPONSE:
                    CollectChildBRUs(this_node, hands_info, children_cfus, this_node_cfu, target_strategy, prune_flag);
                    break;
            }
        }
        return this_node_cfu;

    } else {
        // non-trainee's turn. sample an action using its own strategy.
        float distr_rnb[a_max];
        target_strategy->ComputeStrategy(this_node, b, distr_rnb, resolved_strategy_type);
        int sampled_a = RndXorShift<float>(distr_rnb, a_max, x, y, z, (1 << 16));
        if (sampled_a == -1) {
            target_strategy->PrintNodeStrategy(this_node, b, resolved_strategy_type);
            logger::critical("sampled_a == -1");
        }
        if (!cfr_param_->avg_side_update_ /* set by "side_walk" */ && cfr_param_->rm_avg_update == AVG_CLASSIC) {
            // update wavg. don't use it for blueprint training as it explodes quickly
            //    for (auto a = 0; a < a_max; a++) {
            //      target_strategy->uint_wavg_[rnb0 + a] += distr_rnb[a] * 1000;
            //    }
#ifdef DEBUG_EAGER_LOOKUP
            // TODO(kwok): 🦊
            target_strategy->eager_uint_wavg_[rnb0 + sampled_a] += 1;
#endif
            // TODO(kwok): 🐦
            target_strategy->uint_wavg_->upsert(rnb0 + sampled_a, [&](auto &n) { n++; }, 1);
        }
        return WalkTree(
                this_node->children[sampled_a], trainee, hands_info, target_strategy,
                trainee_cfu_compute_mode_hint, trainee_strategy_type_hint, learn
        );
    }
}

double ScalarCfrWorker::RolloutWalkLeafTreeWithBiasFavor(Node *this_node, int trainee, sPrivateHandsInfo &hands_info,
                                                         int *bias_favors_for_all)
{
#if 0
    this_node->PrintState("leaf node: ");
#endif
    if (this_node->IsTerminal()) {
        return EvalTermNode(this_node, trainee, hands_info);
    }
    return RolloutLeafInterNodeWithBiasFavor(this_node, trainee, hands_info, bias_favors_for_all);
}

double ScalarCfrWorker::RolloutLeafRootNode(Node *leaf_root_node, sPrivateHandsInfo &hands_info)
{
    // NOTE(kwok): Map leaf_root_node to a node from the blueprint. We only consider the blueprint during
    // the rollout.
    NodeMatchResult condition;
    // todo: change the strategy match node
    blueprint_->ag_->MapStateToNode(leaf_root_node->state_, condition);

    auto *matched_node = condition.matched_node_;

    // NOTE(kwok): rollout for `n_rollout_iters` times starting from the matched node till we hit terminals
    sPrivateHandsInfo subgame_priv_hands_info(hands_info.num_players, hands_info.external_sampled_board_, gen);
    // FIXME(kwok): The number of players is not supposed to be fixed to 2.
    subgame_priv_hands_info.internal_sampled_priv_hands_[0] = hands_info.internal_sampled_priv_hands_[0];
    subgame_priv_hands_info.internal_sampled_priv_hands_[1] = hands_info.internal_sampled_priv_hands_[1];

    // NOTE(kwok): fill the board according to the round we are currently at
    auto r = leaf_root_node->GetRound();
    int n_init_board_cards = r == HOLDEM_ROUND_PREFLOP ? 0 : 3;
    for (int c = n_init_board_cards; c < HOLDEM_MAX_BOARD; c++) {
        // fill the remainder of the board cards array with placeholders
        subgame_priv_hands_info.external_sampled_board_.cards[c] = IMPOSSIBLE_CARD;
    }

    int n_rollout_iters = cfr_param_->depth_limited_rollout_reps_;

    // Allocate regrets for each four strategy. Should the regrets be global?
    double all_players_regrets[2][MAX_META_STRATEGY]; // FIXME(kwok): The number of players is not supposed to be fixed to 2.
    for (auto &player_regrets: all_players_regrets) {
        for (auto &regret: player_regrets) {
            regret = 0;
        }
    }

    double rollout_final_cfus[2]; // FIXME(kwok): The number of players is not supposed to be fixed to 2.

    // NOTE(kwok): The fundamental structure of the following giant loop:
    //      3 rollout reps (assumed here)
    //          2 trainees FIXME(kwok): The number of players is not supposed to be fixed to 2.
    //              4 continuation strategies (pre-computed, along with biased towards folding, calling, and raising respectively)
    // doing 3 x 2 x 4 = 24 probes in total
    // TODO(kwok): Separate this loop body into several testable functions.
    for (int rollout_i = 0; rollout_i < n_rollout_iters; rollout_i++) {
        HoldemDeck deck{subgame_priv_hands_info.external_sampled_board_}; // excluding existing public cards
        deck.Shuffle();

        // NOTE(kwok): sample a board for every rollout iteration
        int n_curr_board_cards = n_init_board_cards;
        int deck_cursor = 0;
        while (n_curr_board_cards <= HOLDEM_MAX_BOARD) {
            auto sampled_public_card = deck.cards_[deck_cursor++];
            // FIXME(kwok): The number of players is not supposed to be fixed to 2.
            if (VectorIdxCrashesWithCard(subgame_priv_hands_info.internal_sampled_priv_hands_[0], sampled_public_card)
                || VectorIdxCrashesWithCard(subgame_priv_hands_info.internal_sampled_priv_hands_[1], sampled_public_card
            )) {
                // NOTE(kwok): the sampled board must not crash with the current private hands
                continue;
            }
            subgame_priv_hands_info.external_sampled_board_.cards[n_curr_board_cards++] = sampled_public_card;
        }

        // subgame_priv_hands_info.external_sampled_board_.Print();
        subgame_priv_hands_info.SetBucketAndPayoff(blueprint_->ag_);

        // FIXME(kwok): The number of players is not supposed to be fixed to 2.
        for (int trainee = 0; trainee < 2; trainee++) {
            double trainee_bias_favor_cfus[MAX_META_STRATEGY];

            // NOTE(kwok): This iteration is all about the current trainee. For the opponents, simply pick
            // a bias favor based on their regrets accumulated during the rollout.
            float opp_distr[MAX_META_STRATEGY];
            GetPolicy<double>(opp_distr, MAX_META_STRATEGY, all_players_regrets[1 - trainee]);
            auto opp_bias_favor = RndXorShift<float>(opp_distr, MAX_META_STRATEGY, x, y, z, (1 << 16));
            if (opp_bias_favor == -1) {
                logger::warn("🚨depth limit meta strategy regret problem");
                for (int g = 0; g < MAX_META_STRATEGY; g++) {
                    logger::warn("action %d with the regret of %f", g, all_players_regrets[1 - trainee][g]);
                }
            }

            // NOTE(kwok): probe with each of the four continuation strategies
            for (int trainee_bias_favor = 0; trainee_bias_favor < MAX_META_STRATEGY; trainee_bias_favor++) {
                // FIXME(kwok): The number of players is not supposed to be fixed to 2.
                int bias_favors_for_all[2];
                bias_favors_for_all[trainee] = trainee_bias_favor;
                bias_favors_for_all[1 - trainee] = opp_bias_favor;
                // NOTE(kwok): Fire a probe with the selected bias favor until reached a terminal. All nodes
                // on the decision chain will be biased towards the specified bias favor.
                trainee_bias_favor_cfus[trainee_bias_favor] = RolloutWalkLeafTreeWithBiasFavor(matched_node,
                                                                                               trainee,
                                                                                               subgame_priv_hands_info,
                                                                                               bias_favors_for_all
                );
            }

            float trainee_distr[MAX_META_STRATEGY];
            GetPolicy<double>(trainee_distr, MAX_META_STRATEGY, all_players_regrets[trainee]);
            double trainee_cfu = 0.0;
            for (int bias_favor = 0; bias_favor < MAX_META_STRATEGY; bias_favor++) {
                trainee_cfu += trainee_distr[bias_favor] * trainee_bias_favor_cfus[bias_favor];
            }

            if (rollout_i == n_rollout_iters - 1) {
                // if at the final rollout iteration
                rollout_final_cfus[trainee] = trainee_cfu;
            } else {
                // update the regrets for the current player
                for (int bias_favor = 0; bias_favor < MAX_META_STRATEGY; bias_favor++) {
                    double diff = trainee_bias_favor_cfus[bias_favor] - trainee_cfu;
                    all_players_regrets[trainee][bias_favor] += diff;
                }
            }
        }
    }

    return rollout_final_cfus[0]; // only for player 0
}

double ScalarCfrWorker::RolloutLeafInterNodeWithBiasFavor(Node *this_node, int trainee, sPrivateHandsInfo &hands_info,
                                                          int *bias_favors_for_all)
{
    auto r = this_node->GetRound();
    auto acting_player = this_node->GetActingPlayer();
    auto b = hands_info.buckets_[acting_player][r];
    int a_max = this_node->GetAmax();

    float continuation_distr_rnb[a_max]; // NOTE(kwok): one of the four continuation strategies
    blueprint_->ComputeStrategy(this_node, b, continuation_distr_rnb, STRATEGY_ZIPAVG);
    auto bias_favor = bias_favors_for_all[acting_player];

    // find out the children that call
    int call_sibling_idx = -1;
    for (auto *c: this_node->children) {
        if (c->GetLastAction().type == a_call) {
            call_sibling_idx = c->sibling_idx_;
            break;
        }
    }

    if (call_sibling_idx == -1) {
        this_node->PrintAllChildState();
        logger::critical("calling should always be an available action");
    }

    // NOTE(kwok): biased according to the specified favor
    switch (bias_favor) {
        case BIASED_CALLING: {
            continuation_distr_rnb[call_sibling_idx] *= BIASED_SCALE;
            break;
        }
        case BIASED_RAISING: {
            for (int i = call_sibling_idx + 1; i < a_max; i++) {
                continuation_distr_rnb[i] *= BIASED_SCALE;
            }
            break;
        }
        case BIASED_FOLDING: {
            if (call_sibling_idx == 1) {
                continuation_distr_rnb[0] *= BIASED_SCALE;
            }
            break;
        }
        case BIASED_NONE: {
            break; // NOTE(kwok): use the precomputed blueprint strategy as is
        }
    }

    // NOTE(kwok): re-normalization
    float local_distr_sum = 0.0;
    for (auto a = 0; a < a_max; a++) {
        local_distr_sum += continuation_distr_rnb[a];
    }

    float factor = 1.f / local_distr_sum;
    for (auto a = 0; a < a_max; a++) {
        continuation_distr_rnb[a] *= factor;
    }

    // NOTE(kwok): pick an action
    int sampled_a = RndXorShift<float>(continuation_distr_rnb, a_max, x, y, z, (1 << 16));
    if (sampled_a == -1) {
        logger::debug("💢blueprint zip problem");
        blueprint_->PrintNodeStrategy(this_node, b, STRATEGY_ZIPAVG);
        exit(1);
    }

#if 0
    logger::debug("select child leaf node %d", this_node->children[sampled_a]->GetLastActionCode());
#endif

    auto *sampled_child = this_node->children[sampled_a];
    double cfu = RolloutWalkLeafTreeWithBiasFavor(sampled_child, trainee, hands_info, bias_favors_for_all);
    return cfu;
}

/// Pluribus' way to update WAVG.
void ScalarCfrWorker::WavgUpdateSideWalk(Node *this_node, int trainee, sPrivateHandsInfo &hands_info,
                                         Strategy *target_strategy)
{
    if (this_node->IsTerminal()) {
        return;
    }

    if (this_node->IsLeafNode()) {
        return;
    }

    bool is_my_turn = this_node->GetActingPlayer() == trainee;
    auto a_max = this_node->GetAmax();

    if (is_my_turn) {
        // NOTE(kwok): The agent's turn.
        int r = this_node->GetRound();
        auto n = this_node->GetN();
        Bucket_t b = hands_info.buckets_[trainee][r];
        float distr_rnb[a_max];
        target_strategy->ComputeStrategy(this_node, b, distr_rnb, cfr_param_->strategy_cal_mode_);
        int sampled_a = RndXorShift<float>(distr_rnb, a_max, x, y, z, (1 << 16));
        if (sampled_a == -1) {
            logger::debug("🚨new strategy regret problem");
            target_strategy->PrintNodeStrategy(this_node, b, cfr_param_->strategy_cal_mode_);
        }
        auto rnba = target_strategy->ag_->kernel_->hash_rnba(r, n, b, sampled_a);
#ifdef DEBUG_EAGER_LOOKUP
        // TODO(kwok): 🦊
        target_strategy->eager_uint_wavg_[rnba] += 1;
#endif
        // TODO(kwok): 🐦
        // TODO(kwok): Test if this operation blocks.
        target_strategy->uint_wavg_->upsert(rnba, [](auto &n) { n++; }, 1);
        WavgUpdateSideWalk(this_node->children[sampled_a], trainee, hands_info, target_strategy);
    } else {
        // NOTE(kwok): The opponent's turn.
        for (auto a = 0; a < a_max; a++) {
            WavgUpdateSideWalk(this_node->children[a], trainee, hands_info, target_strategy);
        }
    }
}
