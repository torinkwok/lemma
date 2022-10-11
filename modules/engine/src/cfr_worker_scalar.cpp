#include "cfr_worker.h"

double ScalarCfrWorker::Solve(Board_t board)
{
    auto ag = strategy_->ag_;

    auto active_players = AbstractGame::GetActivePlayerNum();
    auto private_hands_info = sPrivateHandsInfo(active_players, board, gen);
    // generate root belief based on the board
    std::array<
            sPrivateHandBelief *,
            2 // FIXME(kwok): The number of players is not supposed to be fixed to 2.
    > local_root_beliefs{};
    for (int p = 0; p < active_players; p++) {
        // do preparations
        local_root_beliefs[p] = new sPrivateHandBelief(&ag->root_hand_belief_[p]);
        local_root_beliefs[p]->NormalizeExcludeBoard(board);
    }

    static const int REPS = 1000;
    double utility = 0.0;
    for (int loc_iter = 0; loc_iter < REPS; loc_iter++) {
        // NOTE(kwok): On each iteration, we start by sampling all of chance’s actions: the public chance
        // events visible to each player, as well as the private chance events that are visible to only a
        // subset of the players. In poker, this corresponds to randomly choosing the public cards revealed
        // to the players, and the private cards that each player is dealt.
        //
        // NOTE(kwok): Here we're sampling all the private chance events, i.e. the private hands of each player.
        // The public chance event, i.e. the public cards, has been sampled by the invoker of ScalaraCfrWorker::Solve.
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
        // Note that at a terminal node z ∈ Z, it takes O(1) work to determine the utility for player i, u_i(z).

        // NOTE(kwok): Walk down the training tree alternatively.
        for (int trainee = 0; trainee < active_players; trainee++) {
            utility += WalkTree(trainee, ag->root_node_, private_hands_info);
        }
    }

    // Do a sidewalk for updating WAVG
    if (cfr_param_->avg_side_update_ && cfr_param_->rm_avg_update == AVG_CLASSIC) {
        int SIDE_REPS = 50;
        for (int loc_iter = 0; loc_iter < SIDE_REPS; loc_iter++) {
            private_hands_info.SamplePrivateHandsForAll(ag, local_root_beliefs);
            for (int trainee_pos = 0; trainee_pos < active_players; trainee_pos++) {
                WavgUpdateSideWalk(trainee_pos, ag->root_node_, private_hands_info);
            }
        }
    }

    utility /= (REPS * 2 * ag->GetBigBlind());

    for (auto hb: local_root_beliefs) {
        delete hb;
    }

    return utility;
}

/// NOTE(kwok): Where CFR iterations taking place.
double ScalarCfrWorker::WalkTree(int trainee_pos, Node *this_node, sPrivateHandsInfo &hand_info)
{
    if (this_node->IsTerminal()) {
        return EvalTermNode(trainee_pos, this_node, hand_info);
    }
    // NOTE(kwok): depth limited solving, till the local root of next round (rather than the local
    // root of this round)
    if (this_node->IsLeafNode()) {
        return EvalRootLeafNode(trainee_pos, this_node, hand_info);
    }
    return EvalIntermediateChoiceNode(trainee_pos, this_node, hand_info);
}

double ScalarCfrWorker::EvalTermNode(int trainee_pos, Node *this_node, sPrivateHandsInfo &hand_info)
{
    if (this_node->IsShowdown()) {
        int stake = this_node->GetStake(trainee_pos);
        stake *= hand_info.payoff_[trainee_pos]; // tie 0, win 1, lose -1
        return (double) stake;
    } else {
        // fold
        int stake = this_node->GetStake(trainee_pos);
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
double ScalarCfrWorker::EvalRootLeafNode(int trainee_pos, Node *this_node, sPrivateHandsInfo &hand_info)
{
    // NOTE(kwok): Should be the last round. It's the Pluribus way.
    auto r = this_node->GetRound() - 1;
    /*
     * If exists, read from the cache:
     * (PREFLOP) RAM cache from file
     * (FLOP) rollout/NN cache if needed
     */
    auto b0 = 0;
    auto b1 = 0;
    if (cfr_param_->depth_limited_cache_) {
        this_node->CreateValueCacheIfNull();
        if (r > HOLDEM_ROUND_FLOP) {
            logger::critical("We don't do DLS for post-flop");
        }
        // FIXME(kwok): The number of players is not supposed to be fixed to 2.
        b0 = hand_info.buckets_[0][r];
        b1 = hand_info.buckets_[1][r];
        if (this_node->value_cache_->Exists(b0, b1)) {
            auto player0_cfu = this_node->value_cache_->GetValue(b0, b1);
            // FIXME(kwok): The number of players is not supposed to be fixed to 2.
            if (trainee_pos == 1) {
                player0_cfu *= -1;
            }
            return player0_cfu;
        }
    }

    /*
     * ROLLOUT for the current trainee rather than pairwise or alternatively
     */
    double trainee_cfu = LeafRootRollout(trainee_pos, this_node, hand_info);

    /*
     * or use NN
     */

    // insert to cache
    if (cfr_param_->depth_limited_cache_) {
        this_node->value_cache_->SetValue(b0, b1, trainee_cfu);
    }

    // FIXME(kwok): The number of players is not supposed to be fixed to 2.
    // only in terms of player 0
    return trainee_pos == 0 ? trainee_cfu : -trainee_cfu;
}

double ScalarCfrWorker::EvalIntermediateChoiceNode(int trainee, Node *this_node, sPrivateHandsInfo &hand_info)
{
    int acting_player = this_node->GetActingPlayer();
    bool is_trainee_turn = acting_player == trainee;
    auto a_max = this_node->GetAmax();
    auto r = this_node->GetRound();
    auto n = this_node->GetN();
    auto b = hand_info.buckets_[acting_player][r];
    auto rnb0 = strategy_->ag_->kernel_->hash_rnba(r, n, b, 0);

    if (is_trainee_turn) {
        double children_cfus[a_max];
        bool prune_flag[a_max];
        for (auto a = 0; a < a_max; a++) {
            prune_flag[a] = false;
            auto next_node = this_node->children[a];
            // NOTE(kwok): Do pruning if the flag is set. Skip river nodes and nodes leading to terminal.
            if (iter_prune_flag && !next_node->IsTerminal() && next_node->GetRound() != HOLDEM_ROUND_RIVER) {
                // if (iter_prune_flag && !next_node->IsTerminal()) {
                if (strategy_->int_regret_[rnb0 + a] <= cfr_param_->rollout_prune_thres) {
                    prune_flag[a] = true;
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
            double temp_reg = strategy_->int_regret_[rnb0 + a] + diff; // NOTE(kwok): Accumulate regret values
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
            strategy_->int_regret_[rnb0 + a] = new_reg;
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
            strategy_->uint_wavg_[rnb0 + sampled_a] += 1;
        }

        return WalkTree(trainee, this_node->children[sampled_a], hand_info);
    }
}

double ScalarCfrWorker::WalkLeafTree(int trainee_pos,
                                     Node *this_node,
                                     sPrivateHandsInfo &hand_info,
                                     int *c_strategy)
{
#if 0
    this_node->PrintState("leaf node: ");
#endif
    if (this_node->IsTerminal()) {
        return EvalTermNode(trainee_pos, this_node, hand_info);
    }
    return LeafChoiceRollout(trainee_pos, this_node, hand_info, c_strategy);
}

double ScalarCfrWorker::LeafRootRollout(int trainee, Node *this_node, sPrivateHandsInfo &hand_info)
{
    // NOTE(kwok): Map this_node to a node from the blueprint. We only consider the blueprint
    // during the rollout.
    NodeMatchResult condition;
    // todo: change the strategy match node
    blueprint_->ag_->MapStateToNode(this_node->state_, condition);

    auto *matched_node = condition.matched_node_;

    // NOTE(kwok): Rollout for `rollout_rep` times starting from the matched node till we hit terminals.
    sPrivateHandsInfo subgame_priv_hands_info(hand_info.num_players, hand_info.board_, gen);
    // FIXME(kwok): The number of players is not supposed to be fixed to 2.
    subgame_priv_hands_info.hand_[0] = hand_info.hand_[0];
    subgame_priv_hands_info.hand_[1] = hand_info.hand_[1];

    // NOTE(kwok): Fill the board according to the round we are currently at
    auto r = this_node->GetRound();
    int n_init_board_cards = r == HOLDEM_ROUND_PREFLOP ? 0 : 3;
    for (int c = n_init_board_cards; c < HOLDEM_MAX_BOARD; c++) {
        // fill the rest of the board cards array with placeholders
        subgame_priv_hands_info.board_.cards[c] = IMPOSSIBLE_CARD;
    }

    int rollout_rep = cfr_param_->depth_limited_rollout_reps_;

    /*
     * V1: Externally sampling
     */
    // Allocate regrets for each four strategy. Should the regrets be global?
    // FIXME(kwok): The number of players is not supposed to be fixed to 2.
    double all_players_regrets[2][MAX_META_STRATEGY];
    for (auto &player_regrets: all_players_regrets) {
        for (auto &regret: player_regrets) {
            regret = 0;
        }
    }

    double final_cfus[2]; // FIXME(kwok): The number of players is not supposed to be fixed to 2.

    /*
     * reps 3
     *      traverser 2 FIXME(kwok): The number of players is not supposed to be fixed to 2.
     *          strategy 4
     * doing 24 probing in total
     */
    for (int rollout_i = 0; rollout_i < rollout_rep; rollout_i++) {
        // TODO(kwok): Make this loop body a testable function.
        // NOTE(kwok): Sample a board, which must not crash with the current hands.
        HoldemDeck deck{subgame_priv_hands_info.board_};
        deck.Shuffle();
        int n_curr_board_cards = n_init_board_cards;
        int deck_cursor = 0;
        while (n_curr_board_cards <= HOLDEM_MAX_BOARD) {
            auto sampled_public_card = deck.cards_[deck_cursor++];
            if (VectorIdxCrashesWithCard(subgame_priv_hands_info.hand_[0], sampled_public_card)
                || VectorIdxCrashesWithCard(subgame_priv_hands_info.hand_[1], sampled_public_card)) {
                continue;
            }
            subgame_priv_hands_info.board_.cards[n_curr_board_cards++] = sampled_public_card;
        }

        // subgame_priv_hands_info.board_.Print();
        subgame_priv_hands_info.SetBucketAndPayoff(blueprint_->ag_);

        // FIXME(kwok): The number of players is not supposed to be fixed to 2.
        for (int p = 0; p < 2; p++) {
            // NOTE(kwok): pick an action for the opponent based on its strategy profile
            float opp_distr[MAX_META_STRATEGY];
            GetPolicy<double>(opp_distr, MAX_META_STRATEGY, all_players_regrets[1 - p]);
            auto opp_sampled_a = RndXorShift<float>(opp_distr, MAX_META_STRATEGY, x, y, z, (1 << 16));
            if (opp_sampled_a == -1) {
                logger::warn("🚨depth limit meta strategy regret problem");
                for (int g = 0; g < MAX_META_STRATEGY; g++) {
                    logger::warn("action %d with the regret of %f", g, all_players_regrets[1 - p][g]);
                }
            }

            double action_cfus[MAX_META_STRATEGY];

            // For each strategy of the current acting player
            for (int s = 0; s < MAX_META_STRATEGY; s++) {
                // FIXME(kwok): The number of players is not supposed to be fixed to 2.
                int c_strategy[2];
                c_strategy[p] = s;
                c_strategy[1 - p] = opp_sampled_a;
                // rollout with the strategy combination
                action_cfus[s] = WalkLeafTree(p, matched_node, subgame_priv_hands_info, c_strategy);
            }

            // NOTE(kwok): pick an action for the trainee p
            float trainee_distr[MAX_META_STRATEGY];
            GetPolicy<double>(trainee_distr, MAX_META_STRATEGY, all_players_regrets[p]);
            double cfu = 0.0;
            for (int s = 0; s < MAX_META_STRATEGY; s++) {
                cfu += trainee_distr[s] * action_cfus[s];
            }

            // If at the final iter
            if (rollout_i == rollout_rep - 1) {
                final_cfus[p] = cfu;
            } else {
                // update the regrets for the current player
                for (int s = 0; s < MAX_META_STRATEGY; s++) {
                    double diff = action_cfus[s] - cfu;
                    all_players_regrets[p][s] += diff;
                }
            }
        }
    }

    return final_cfus[0]; // only for player 0
}

double ScalarCfrWorker::LeafChoiceRollout(int trainee_pos,
                                          Node *this_node,
                                          sPrivateHandsInfo &hand_info,
                                          int *p_meta)
{
    auto r = this_node->GetRound();
    auto acting_player = this_node->GetActingPlayer();
    auto b = hand_info.buckets_[acting_player][r];
    int a_max = this_node->GetAmax();

    float distr_rnb[a_max];
    blueprint_->ComputeStrategy(this_node, b, distr_rnb, STRATEGY_ZIPAVG);

    // adjust according to the biased strategy
    auto continuation_strategy = p_meta[acting_player];
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

    switch (continuation_strategy) {
        case BIASED_CALLING: {
            distr_rnb[call_sibling_idx] *= BIASED_SCALER;
            break;
        }
        case BIASED_RAISING: {
            for (int i = call_sibling_idx + 1; i < a_max; i++) {
                distr_rnb[i] *= BIASED_SCALER;
            }
            break;
        }
        case BIASED_FOLDING: {
            if (call_sibling_idx == 1) {
                distr_rnb[0] *= BIASED_SCALER;
            }
            break;
        }
        case BIASED_NONE: {
            // do nothing
            break;
        }
    }

    // renomalize the distribution
    float local_distr_sum = 0.0;
    for (auto a = 0; a < a_max; a++) local_distr_sum += distr_rnb[a];
    float scaler = 1.f / local_distr_sum;
    for (auto a = 0; a < a_max; a++) distr_rnb[a] *= scaler;

    // pick an action
    int sampled_a = RndXorShift<float>(distr_rnb, a_max, x, y, z, (1 << 16));
    if (sampled_a == -1) {
        logger::debug("💢blueprint zip problem");
        blueprint_->PrintNodeStrategy(this_node, b, STRATEGY_ZIPAVG);
        exit(1);
    }

#if 0
    logger::debug("select child leaf node %d", this_node->children[sampled_a]->GetLastActionCode());
#endif

    double cfu = WalkLeafTree(trainee_pos, this_node->children[sampled_a], hand_info, p_meta);
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
        strategy_->uint_wavg_[rnba] += 1;
        WavgUpdateSideWalk(trainee_pos, this_node->children[sampled_a], hand_info);
    } else {
        // NOTE(kwok): The opponent's turn.
        for (auto a = 0; a < a_max; a++) {
            WavgUpdateSideWalk(trainee_pos, this_node->children[a], hand_info);
        }
    }
}
