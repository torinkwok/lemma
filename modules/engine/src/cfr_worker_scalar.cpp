#include "cfr_worker.h"

double ScalarCfrWorker::Solve(Board_t board)
{
    auto ag = strategy_->ag_;
    auto active_players = AbstractGame::GetActivePlayerNum();
    auto private_hands_info = sPrivateHandsInfo(active_players, board, gen);

    // NOTE(kwok): assemble local root private hand beliefs based on the given board
    std::array<
            sPrivateHandBelief *,
            2 // FIXME(kwok): The number of players is not supposed to be fixed to 2.
    > local_root_beliefs{};
    for (int p = 0; p < active_players; p++) {
        local_root_beliefs[p] = new sPrivateHandBelief(&ag->root_hand_beliefs_for_all_[p]);
        local_root_beliefs[p]->NormalizeExcludeBoard(board);
    }

    double sum_cfus = 0.0;
    static const int n_priv_hand_samples = 1; // TODO(kwok): What about blueprint computing?
    for (int i = 0; i < n_priv_hand_samples; i++) {
        // NOTE(kwok): On each iteration, we start by sampling all of chance’s actions: the public chance
        // events visible to each player, as well as the private chance events that are visible to only a
        // subset of the players. In poker, this corresponds to randomly choosing the public cards revealed
        // to the players, and the private cards that each player is dealt.
        //
        // NOTE(kwok): Here we're sampling all the private chance events, i.e. the private hands of each player.
        // The public chance event, i.e. the public cards, has already been sampled by the outside invoker of
        // ScalaraCfrWorker::Solve.
        private_hands_info.SamplePrivateHandsForAll(ag, local_root_beliefs);
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
        // Note that at a terminal node z ∈ Z, it takes O(1) work to determine the sum_cfus for player i, u_i(z).

        // NOTE(kwok): Walk down the training tree alternatively.
        for (int trainee = 0; trainee < active_players; trainee++) {
            sum_cfus += WalkTree(trainee, ag->root_node_, private_hands_info);
        }
    }

    // conduct a sidewalk for updating WAVG
    if (cfr_param_->avg_side_update_ && cfr_param_->rm_avg_update == AVG_CLASSIC) {
        int n_sidewalk_iters = 50;
        for (int i = 0; i < n_sidewalk_iters; i++) {
            private_hands_info.SamplePrivateHandsForAll(ag, local_root_beliefs);
            for (int trainee_pos = 0; trainee_pos < active_players; trainee_pos++) {
                WavgUpdateSideWalk(trainee_pos, ag->root_node_, private_hands_info);
            }
        }
    }

    sum_cfus /= (n_priv_hand_samples * 2 * ag->GetBigBlind());

    for (auto b: local_root_beliefs) {
        delete b;
    }

    return sum_cfus;
}

/// NOTE(kwok): Where CFR iterations taking place.
double ScalarCfrWorker::WalkTree(int trainee, Node *this_node, sPrivateHandsInfo &hand_info)
{
    if (this_node->IsTerminal()) {
        return EvalTermNode(trainee, this_node, hand_info);
    }
    // NOTE(kwok): depth limited solving, till the local root of next round (rather than the local
    // root of this round)
    if (this_node->IsLeafNode()) {
        return EvalLeafRootNode(trainee, this_node, hand_info);
    }
    return EvalInterNode(trainee, this_node, hand_info);
}

double ScalarCfrWorker::EvalTermNode(int trainee, Node *this_node, sPrivateHandsInfo &hand_info)
{
    if (this_node->IsShowdown()) {
        int stake = this_node->GetStake(trainee);
        int payoff = hand_info.payoff_[trainee];
        stake *= payoff; // tie 0, win 1, lose -1
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
double ScalarCfrWorker::EvalLeafRootNode(int trainee, Node *leaf_root_node, sPrivateHandsInfo &hand_info)
{
    // NOTE(kwok): Should be the last round. It's the Pluribus way.
    auto r = leaf_root_node->GetRound() - 1;
    /*
     * If exists, read from the cache:
     * (PREFLOP) RAM cache from file
     * (FLOP) rollout/NN cache if needed
     */
    auto b0 = 0;
    auto b1 = 0;
    if (cfr_param_->depth_limited_cache_) {
        leaf_root_node->CreateValueCacheIfNull();
        if (r > HOLDEM_ROUND_FLOP) {
            logger::critical("We don't do DLS for post-flop");
        }
        // FIXME(kwok): The number of players is not supposed to be fixed to 2.
        b0 = hand_info.buckets_[0][r];
        b1 = hand_info.buckets_[1][r];
        if (leaf_root_node->value_cache_->Exists(b0, b1)) {
            auto player0_cfu = leaf_root_node->value_cache_->GetValue(b0, b1);
            // FIXME(kwok): The number of players is not supposed to be fixed to 2.
            if (trainee == 1) {
                player0_cfu *= -1;
            }
            return player0_cfu;
        }
    }

    /*
     * ROLLOUT for the current trainee rather than pairwise or alternatively
     */
    double player0_rollout_cfu = RolloutLeafRootNode(leaf_root_node, hand_info);

    /*
     * or use NN
     */

    // insert to cache
    if (cfr_param_->depth_limited_cache_) {
        leaf_root_node->value_cache_->SetValue(b0, b1, player0_rollout_cfu);
    }

    // FIXME(kwok): The number of players is not supposed to be fixed to 2.
    // only in terms of player 0
    return trainee == 0 ? player0_rollout_cfu : -player0_rollout_cfu;
}

double ScalarCfrWorker::EvalInterNode(int trainee, Node *this_node, sPrivateHandsInfo &hand_info)
{
    int acting_player = this_node->GetActingPlayer();
    bool is_trainee_turn = acting_player == trainee;
    auto a_max = this_node->GetAmax();
    auto r = this_node->GetRound();
    auto n = this_node->GetN();
    auto b = hand_info.buckets_[acting_player][r];
    auto rnb0 = strategy_->ag_->kernel_->hash_rnba(r, n, b, 0);
    // RNBA debug_copy = rnb0;
    // printf("debug copy = %llu", debug_copy);

    if (is_trainee_turn) {
        double children_cfus[a_max];
        bool prune_flag[a_max];
        for (auto a = 0; a < a_max; a++) {
            prune_flag[a] = false;
            auto next_node = this_node->children[a];
            // NOTE(kwok): Do pruning if the flag is set. Skip river nodes and nodes leading to terminal.
            if (iter_prune_flag && !next_node->IsTerminal() && next_node->GetRound() != HOLDEM_ROUND_RIVER) {
                // if (iter_prune_flag && !next_node->IsTerminal()) {
                // if (strategy_->int_regret_[rnb0 + a] <= cfr_param_->rollout_prune_thres) {

                // if (strategy_->int_regret_->try_emplace(rnb0 + a).first->second <= cfr_param_->rollout_prune_thres) {
                //     prune_flag[a] = true;
                //     continue;
                // }

                // TODO(kwok): ❓
                // TODO(kwok): Test if this operation blocks.
                strategy_->int_regret_->insert(rnb0 + a);
                strategy_->int_regret_->find_fn(rnb0 + a, [&](const auto &regret)
                                                {
                                                    // TODO(kwok): ❓
                                                    if (regret <= cfr_param_->rollout_prune_thres) {
                                                        prune_flag[a] = true;
                                                    }
                                                }
                );
                // TODO(kwok): ❓
                if (prune_flag[a]) {
                    continue;
                }
            }
            children_cfus[a] = WalkTree(trainee, next_node, hand_info);
        }

        // NOTE(kwok): only WEIGHTED_RESPONSE is supported
        if (cfr_param_->cfu_compute_acting_playing != WEIGHTED_RESPONSE) {
            logger::critical("scalar does not support best response eval");
        }

        double this_node_cfu = 0.0;

        float distr_rnb[a_max];
        // NOTE(kwok): query RNBA indexed strategy data to compute this_node_cfu
        strategy_->ComputeStrategy(this_node, b, distr_rnb, cfr_param_->strategy_cal_mode_);
        for (auto a = 0; a < a_max; a++) {
            if (prune_flag[a]) {
                // NOTE(kwok): Pruned nodes won't be taken into account when calculating the final CFU
                // of the current node.
                continue;
            }
            this_node_cfu += distr_rnb[a] * children_cfus[a];
        }

        for (auto a = 0; a < a_max; a++) {
            if (prune_flag[a]) {
                continue;
            }
            int diff = (int) round(children_cfus[a] - this_node_cfu); // NOTE(kwok): The regret value
            // double temp_reg = strategy_->int_regret_[rnb0 + a] + diff; // NOTE(kwok): Accumulate regret values
            // double temp_reg = strategy_->int_regret_->try_emplace(rnb0 + a).first->second +
            //                   diff; // NOTE(kwok): Accumulate regret values

            // double temp_reg = strategy_->int_regret_->try_emplace(rnb0 + a).first->second +
            //                   diff; // NOTE(kwok): Accumulate regret values
            // // clamp it
            // int new_reg = (int) std::fmax(temp_reg, cfr_param_->rm_floor);
            // // total regret should have a ceiling
            // if (new_reg > 2107483647) { // 2147483647 - 4 * 10^7
            //     logger::critical(
            //             "if the regret overflowed, think about the possibility of the regret %d not being converging. [temp_reg %f][diff %f][cfu_a %f][this_node_cfu %f]",
            //             new_reg,
            //             temp_reg,
            //             diff,
            //             children_cfus[a],
            //             this_node_cfu
            //     );
            // }
            // // strategy_->int_regret_[rnb0 + a] = new_reg;
            // strategy_->int_regret_->operator[](rnb0 + a) = new_reg;

            // TODO(kwok): ❓
            strategy_->int_regret_->insert(rnb0 + a);
            strategy_->int_regret_->update_fn(rnb0 + a, [&](auto &regret)
                                              {
                                                  double temp_reg = regret + diff; // NOTE(kwok): Accumulate regret values
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

        return this_node_cfu;
    } else {
        // NOTE(kwok): Non-trainee's turn. Sample an action using their own computed strategy.
        float distr_rnb[a_max];
        strategy_->ComputeStrategy(this_node, b, distr_rnb, cfr_param_->strategy_cal_mode_);
        int sampled_a = RndXorShift<float>(distr_rnb, a_max, x, y, z, (1 << 16));
        if (sampled_a == -1) {
            strategy_->PrintNodeStrategy(this_node, b, cfr_param_->strategy_cal_mode_);
            logger::critical("new strategy regret problem in main walk, opp side");
        }

        if (!cfr_param_->avg_side_update_ /* set by "side_walk" */ && cfr_param_->rm_avg_update == AVG_CLASSIC) {
            // update wavg. don't use it for blueprint training as it explodes quickly
            //    for (auto a = 0; a < a_max; a++) {
            //      strategy_->uint_wavg_[rnb0 + a] += distr_rnb[a] * 1000;
            //    }
            // strategy_->uint_wavg_[rnb0 + sampled_a] += 1;

            // strategy_->uint_wavg_->try_emplace(rnb0 + sampled_a).first->second += 1;

            // TODO(kwok): ❓
            strategy_->uint_wavg_->upsert(rnb0 + sampled_a, [&](auto &n) { n++; }, 1);
        }

        return WalkTree(trainee, this_node->children[sampled_a], hand_info);
    }
}

double ScalarCfrWorker::RolloutWalkLeafTreeWithBiasFavor(int trainee,
                                                         Node *this_node,
                                                         sPrivateHandsInfo &hand_info,
                                                         int *bias_favors_for_all)
{
#if 0
    this_node->PrintState("leaf node: ");
#endif
    if (this_node->IsTerminal()) {
        return EvalTermNode(trainee, this_node, hand_info);
    }
    return RolloutLeafInterNodeWithBiasFavor(trainee, this_node, hand_info, bias_favors_for_all);
}

double ScalarCfrWorker::RolloutLeafRootNode(Node *leaf_root_node, sPrivateHandsInfo &hand_info)
{
    // NOTE(kwok): Map leaf_root_node to a node from the blueprint. We only consider the blueprint during
    // the rollout.
    NodeMatchResult condition;
    // todo: change the strategy match node
    blueprint_->ag_->MapStateToNode(leaf_root_node->state_, condition);

    auto *matched_node = condition.matched_node_;

    // NOTE(kwok): rollout for `n_rollout_iters` times starting from the matched node till we hit terminals
    sPrivateHandsInfo subgame_priv_hands_info(hand_info.num_players, hand_info.external_sampled_board_, gen);
    // FIXME(kwok): The number of players is not supposed to be fixed to 2.
    subgame_priv_hands_info.internal_sampled_priv_hands_[0] = hand_info.internal_sampled_priv_hands_[0];
    subgame_priv_hands_info.internal_sampled_priv_hands_[1] = hand_info.internal_sampled_priv_hands_[1];

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
    //          2 trainees (a.k.a. traversers) FIXME(kwok): The number of players is not supposed to be fixed to 2.
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
                trainee_bias_favor_cfus[trainee_bias_favor] = RolloutWalkLeafTreeWithBiasFavor(
                        trainee, matched_node, subgame_priv_hands_info, bias_favors_for_all
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

double ScalarCfrWorker::RolloutLeafInterNodeWithBiasFavor(int trainee,
                                                          Node *this_node,
                                                          sPrivateHandsInfo &hand_info,
                                                          int *bias_favors_for_all)
{
    auto r = this_node->GetRound();
    auto acting_player = this_node->GetActingPlayer();
    auto b = hand_info.buckets_[acting_player][r];
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
            continuation_distr_rnb[call_sibling_idx] *= BIASED_SCALER;
            break;
        }
        case BIASED_RAISING: {
            for (int i = call_sibling_idx + 1; i < a_max; i++) {
                continuation_distr_rnb[i] *= BIASED_SCALER;
            }
            break;
        }
        case BIASED_FOLDING: {
            if (call_sibling_idx == 1) {
                continuation_distr_rnb[0] *= BIASED_SCALER;
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
    double cfu = RolloutWalkLeafTreeWithBiasFavor(trainee, sampled_child, hand_info, bias_favors_for_all);
    return cfu;
}

/// Pluribus' way to update WAVG.
void ScalarCfrWorker::WavgUpdateSideWalk(int trainee_pos, Node *this_node, sPrivateHandsInfo &hand_info)
{
    if (this_node->IsTerminal()) {
        return;
    }

    if (this_node->IsLeafNode()) {
        return;
    }

    bool is_my_turn = this_node->GetActingPlayer() == trainee_pos;
    auto a_max = this_node->GetAmax();

    if (is_my_turn) {
        // NOTE(kwok): The agent's turn.
        int r = this_node->GetRound();
        auto n = this_node->GetN();
        Bucket_t b = hand_info.buckets_[trainee_pos][r];
        float distr_rnb[a_max];
        strategy_->ComputeStrategy(this_node, b, distr_rnb, cfr_param_->strategy_cal_mode_);
        int sampled_a = RndXorShift<float>(distr_rnb, a_max, x, y, z, (1 << 16));
        if (sampled_a == -1) {
            logger::debug("🚨new strategy regret problem");
            strategy_->PrintNodeStrategy(this_node, b, cfr_param_->strategy_cal_mode_);
        }
        auto rnba = strategy_->ag_->kernel_->hash_rnba(r, n, b, sampled_a);
        // strategy_->uint_wavg_[rnba] += 1;
        // strategy_->uint_wavg_->try_emplace(rnba).first->second += 1;
        // TODO(kwok): ❓
        // TODO(kwok): Test if this operation blocks.
        strategy_->uint_wavg_->upsert(rnba, [](auto &n) { n++; }, 1);
        WavgUpdateSideWalk(trainee_pos, this_node->children[sampled_a], hand_info);
    } else {
        // NOTE(kwok): The opponent's turn.
        for (auto a = 0; a < a_max; a++) {
            WavgUpdateSideWalk(trainee_pos, this_node->children[a], hand_info);
        }
    }
}
