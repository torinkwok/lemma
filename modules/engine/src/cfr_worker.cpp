#include "cfr_worker.h"
#include "node.h"

double VectorCfrWorker::Solve(Board_t board)
{
    auto ag = strategy_->ag_;
    auto starting_round = ag->root_state_.round;
    if (starting_round < 3 /* pre-flop, flop, and turn */) {
        delete priv_hand_kernel;
        priv_hand_kernel = new sPrivateHandKernel(board, starting_round);
        priv_hand_kernel->AbstractHandKernel(&ag->bucket_reader_);
    } else if (starting_round == 3 /* river */ && priv_hand_kernel == nullptr) {
        // cache the hand kernel if just solving round RIVER subgame
        priv_hand_kernel = new sPrivateHandKernel(board, starting_round);
        priv_hand_kernel->AbstractHandKernel(&ag->bucket_reader_);
    }

    // pruning with prob
    ConditionalPrune();

    // TODO(kwok): Cache the results for RIVER subgames.
    sPrivateHandBelief local_root_belief[2];
    auto active_players = ag->GetActivePlayerNum();
    for (int p = 0; p < active_players; p++) {
        // preparations
        local_root_belief[p].CopyValue(&ag->root_hand_beliefs_for_all_[p]);
        local_root_belief[p].NormalizeExcludeBoard(board);
    }

    double sum_cfus = 0.0;

    // NOTE(kwok): There are three different methods of sampling chance events that have slower iterations,
    // but do more work on each iteration:
    //
    //   0. Chance-Sampling
    //      * ⬇ (A SCALAR for us, A SCALAR for the opponent)
    //      * ⬆ A SCALAR (the sampled counterfactual value v ̃_i(σ, I) for player i)
    //
    //   1. Opponent-Public Chance Sampling
    //      * ⬇ (A VECTOR for us, A SCALAR for the opponent)
    //      * ⬆ A VECTOR(our counterfactual value for each of our private chance outcomes)
    //
    //   2. Self-Public Chance Sampling
    //      * ⬇ (A SCALAR for us, A VECTOR for the opponent)
    //      * ⬆ A SCALAR (the counterfactual value for our sampled outcome)
    //
    //   3. Public Chance Sampling
    //      * ⬇ (A VECTOR for us, A VECTOR for the opponent)
    //      * ⬆ A VECTOR (containing the counterfactual value for each of our n information set)
    //      * At the terminal nodes, we seemingly have an O(n^2) computation, as for each of our n information
    //        sets, we must consider all n of the opponent’s possible private outcomes in order to compute our
    //        utility for that information set. However, if the payoffs at terminal nodes are structured in some
    //        way, we can often reduce this to an O(n) evaluation that returns exactly the same value as the
    //        O(n^2) evaluation. Doing so gives PCS the advantage of both SPCS (accurate strategy updates) and
    //        OPCS (many strategy updates) for the same evaluation cost of either.
    switch (mode_) {
        case CFR_VECTOR_PAIRWISE_SOLVE: {
            Ranges starting_ranges{active_players};
            for (int p = 0; p < starting_ranges.num_player_; p++) {
                starting_ranges.beliefs_[p].CopyValue(&local_root_belief[p]);
                starting_ranges.beliefs_[p].Scale(REGRET_SCALER);
            }
            auto cfu = WalkTree_Pairwise(ag->root_node_, &starting_ranges);
            for (int p = 0; p < starting_ranges.num_player_; p++) {
                cfu->beliefs_[p].DotMultiply(&local_root_belief[p]);
            }
            sum_cfus = cfu->ValueSum();
            // logger::debug("%f | %f", local_root_belief[0].BeliefSum(), local_root_belief[1].BeliefSum());
            delete cfu;
            break;
        }
        case CFR_VECTOR_ALTERNATE_SOLVE: {
            // NOTE(kwok): walk down the vector training tree alternatively
            for (int trainee = 0; trainee < active_players; trainee++) {
                // FIXME(kwok): The number of players is not supposed to be fixed to 2.
                auto opp_local_root_belief = local_root_belief[1 - trainee];
                opp_local_root_belief.Scale(REGRET_SCALER);
                sPrivateHandBelief *cfu_p = WalkTree_Alternate(ag->root_node_, trainee, &opp_local_root_belief);
                cfu_p->DotMultiply(&local_root_belief[trainee]);  // illegal parts are 0, so it is fine.
                sum_cfus += cfu_p->BeliefSum();
                delete cfu_p;
            }
            break;
        }
        default: {
            logger::critical("invalid CFR mode %d", mode_);
        }
    }

    // FIXME(kwok): The number of players is not supposed to be fixed to 2.
    double avg_cfu = sum_cfus / (2 * ag->GetBigBlind());
    return avg_cfu;
}

// TODO(kwok): Adding back the pruned is necessary.
// We ignore those pruned hands anyway, hence no need to set -1. The program won't panic.
sPrivateHandBelief *VectorCfrWorker::WalkTree_Alternate(Node *this_node, int trainee, sPrivateHandBelief *opp_belief)
{
    if (opp_belief->AllZero()) {
        return new sPrivateHandBelief(0.0);
    }
    if (this_node->IsTerminal()) {
        auto cfu = new sPrivateHandBelief(0.0);
        priv_hand_kernel->hand_eval_kernel.FastTerminalEval(opp_belief->belief_,
                                                            cfu->belief_,
                                                            this_node->GetStake(trainee),
                                                            this_node->IsShowdown());
        return cfu;
    }
    return EvalChoiceNode_Alternate(this_node, trainee, opp_belief);
}

sPrivateHandBelief *
VectorCfrWorker::EvalChoiceNode_Alternate(Node *this_node, int trainee, sPrivateHandBelief *opp_belief)
{
    auto a_max = this_node->GetAmax();
    bool is_my_turn = trainee == this_node->GetActingPlayer();

    /*
     * ROLLOUT
     */
    // copy the reach_ranges to child ranges
    std::vector<sPrivateHandBelief *> child_reach_ranges;
    child_reach_ranges.reserve(a_max);
    if (is_my_turn) {
        // shallow copy
        for (int a = 0; a < a_max; a++) {
            child_reach_ranges.emplace_back(opp_belief);
        }
    } else {
        // deep copy and rollout
        for (int a = 0; a < a_max; a++) {
            child_reach_ranges.emplace_back(new sPrivateHandBelief(opp_belief));
        }
        RangeRollout(this_node, opp_belief, child_reach_ranges);
    }

    /*
     * WALK DOWN THE TREE RECURSIVELY and collect the CFUs of all children
     */

    std::vector<sPrivateHandBelief *> children_cfus;
    children_cfus.reserve(a_max);
    for (int a = 0; a < a_max; a++) {
        auto *child_node = this_node->children[a];
        auto *child_node_reach_range = child_reach_ranges[a];
        children_cfus.emplace_back(WalkTree_Alternate(child_node, trainee, child_node_reach_range));
    }

    /*
     * COMPUTE CFU of this node
     */

    // Precompute strategy only if it's my turn. Range rollout is fine, it is computed on the fly.
    float *all_belief_distr_1dim = nullptr;
    if (is_my_turn) {
        // unless we have pruning, we don't need to skip any one. it is a bit wasteful but makes it more accurate.
        all_belief_distr_1dim = new float[FULL_HAND_BELIEF_SIZE * a_max];
        for (auto &i: priv_hand_kernel->valid_priv_hand_vector_idxes) {
            int offset = a_max * i;
            auto b = priv_hand_kernel->GetBucket(this_node->GetRound(), i);
            strategy_->ComputeStrategy(this_node, b, all_belief_distr_1dim + offset, cfr_param_->strategy_cal_mode_);
        }
    }

    // NOTE(kwok): A utility of +1 is given for a win, and −1 for a loss.
    auto this_node_cfu = new sPrivateHandBelief(0.0);
    CFU_COMPUTE_MODE compute_mode = is_my_turn
                                    ? cfr_param_->cfu_compute_acting_playing
                                    : cfr_param_->cfu_compute_opponent;
    ComputeCfu(this_node, child_reach_ranges, children_cfus, this_node_cfu, compute_mode, all_belief_distr_1dim);

    /*
     * REGRET LEARNING
     */
    // only learning on the trainee's node
    if (cfr_param_->regret_learning_on) {
        if (is_my_turn) {
            CollectRegrets(this_node, children_cfus, this_node_cfu);
        }
    }

    // delete child pop up this_node_cfu
    for (int a = 0; a < a_max; a++) {
        delete children_cfus[a];
        if (!is_my_turn) {
            delete child_reach_ranges[a];
        }
    }

    if (is_my_turn) {
        delete[] all_belief_distr_1dim;
    }

    return this_node_cfu;
}

// ⬇ (A VECTOR for us, A VECTOR for the opponent)
// ⬆ A VECTOR (containing the counterfactual value for each of our n information set)
Ranges *VectorCfrWorker::WalkTree_Pairwise(Node *this_node, Ranges *reach_ranges)
{
    if (reach_ranges->ReturnReady(iter_prune_flag)) {
        return new Ranges(reach_ranges, "empty");
    }

    if (this_node->IsTerminal()) {
        // some range are 0. do some early return processing
        auto this_cfu = new Ranges(reach_ranges); // copy the masks as well.
        for (auto my_pos = 0; my_pos < reach_ranges->num_player_; my_pos++) {
            priv_hand_kernel->hand_eval_kernel.FastTerminalEval(
                    reach_ranges->beliefs_[
                            // FIXME(kwok): The number of players is not supposed to be fixed to 2.
                            1 - my_pos
                    ].belief_,
                    this_cfu->beliefs_[my_pos].belief_,
                    this_node->GetStake(my_pos),
                    this_node->IsShowdown()
            );
        }
#if DEV > 1
        // FIXME(kwok): The number of players is not supposed to be fixed to 2.
        for (int p = 0; p < 2; p++) {
            if (!this_cfu->beliefs_[p].TopoAligned(&reach_ranges->beliefs_[p])) {
                logger::error("misaligned belief %d topo %d != %d", p,
                              this_cfu->beliefs_[p].CountPrunedEntries(),
                              reach_ranges->beliefs_[p].CountPrunedEntries());
            }
        }
#endif
        return this_cfu;
    }

    return EvalChoiceNode_Pairwise(this_node, reach_ranges);
}

Ranges *VectorCfrWorker::EvalChoiceNode_Pairwise(Node *this_node, Ranges *reach_ranges)
{
    auto a_max = this_node->GetAmax();
    /*
     * ROLLOUT
     */
    //copy the reach_ranges to child ranges for both players
    std::vector<Ranges *> child_reach_ranges;
    child_reach_ranges.reserve(a_max);
    for (int a = 0; a < a_max; a++) {
        child_reach_ranges.emplace_back(new Ranges(reach_ranges));
    }

    //only need to roll out the acting_player
    int actor = this_node->GetActingPlayer();
    auto actor_child_beliefs = ExtractBeliefs(child_reach_ranges, actor);
    RangeRollout(this_node, &reach_ranges->beliefs_[actor], actor_child_beliefs);

    /*
     * WALK TREE RECURSE DOWN
     */
    std::vector<Ranges *> child_cfu;
    child_cfu.reserve(a_max);
    for (int a = 0; a < a_max; a++) {
        child_cfu.emplace_back(WalkTree_Pairwise(this_node->children[a], child_reach_ranges[a]));
    }

    /*
     * COMPUTE CFU
     */
    auto cfu = new Ranges(reach_ranges, "empty");
    //acting player
    auto actor_child_cfu = ExtractBeliefs(child_cfu, actor);
    //todo: need to add the memory thing here.
    ComputeCfu(this_node,
               actor_child_beliefs,
               actor_child_cfu,
               &cfu->beliefs_[actor],
               cfr_param_->cfu_compute_acting_playing, nullptr
    );
    //opponent
    int opp_pos = 1 - actor;
    ComputeCfu(this_node,
               ExtractBeliefs(child_reach_ranges, opp_pos),
               ExtractBeliefs(child_cfu, opp_pos),
               &cfu->beliefs_[opp_pos],
               cfr_param_->cfu_compute_opponent, nullptr
    );

    /*
     * REGRET LEARNING
     */
    //only on actor
    if (cfr_param_->regret_learning_on) {
        CollectRegrets(this_node, actor_child_cfu, &cfu->beliefs_[actor]);
    }

    //delete the util belief from child nodes
    for (int a = 0; a < a_max; a++) {
        // pop up
        delete child_cfu[a];
        //goes down
        delete child_reach_ranges[a];
    }
    return cfu;
}

/*
 * For both side,
 * - do belief_distr update, and do pruning
 * - do wavg update
 */
void VectorCfrWorker::RangeRollout(Node *this_node, sPrivateHandBelief *belief_distr,
                                   std::vector<sPrivateHandBelief *> &child_ranges)
{
    auto r = this_node->GetRound();
    auto a_max = this_node->GetAmax();
    auto n = this_node->GetN();
    auto frozen_b = this_node->frozen_b;

    for (auto &combo_index: priv_hand_kernel->valid_priv_hand_vector_idxes) {
        // greedily outer skip the 0 belief and board cards
        if (belief_distr->IsPruned(combo_index)) {
            continue;
        }

        // the same if zeroed
        if (belief_distr->IsZero(combo_index)) {
            continue;
        }

        auto b = priv_hand_kernel->GetBucket(r, combo_index);
        auto rnb0 = strategy_->ag_->kernel_->hash_rnba(r, n, b, 0);

        // NOTE(kwok): The strategy probabilities if the opponent's private hand was combo_index.
        float distr_rnb[a_max];
        strategy_->ComputeStrategy(this_node, b, distr_rnb, cfr_param_->strategy_cal_mode_);

        /*
         * update WAVG, if not frozen.
         */
        if (cfr_param_->rm_avg_update == AVG_CLASSIC && b != frozen_b) {
            double reach_i = belief_distr->belief_[combo_index] * pow(10, this_node->reach_adjustment[b]);

            // adjustment adaptively goes up
            if (reach_i > 0.0 && reach_i < 100) {
                pthread_mutex_lock(&this_node->mutex_);
                reach_i = belief_distr->belief_[combo_index] * pow(10, this_node->reach_adjustment[b]);
                while (true) {
                    if (reach_i > 100) {
                        break;
                    }
                    this_node->reach_adjustment[b] += 1;
                    reach_i *= 10;
                    for (int a = 0; a < a_max; a++) {
                        strategy_->ulong_wavg_[rnb0 + a] *= 10;
                        // this_node->ulong_wavg_[this_node->HashBa(b, a)] *= 10;
                    }
                }
                pthread_mutex_unlock(&this_node->mutex_);
            }

            // adjustment adaptively goes down
            if (reach_i > pow(10, 15)) {
                pthread_mutex_lock(&this_node->mutex_);
                // get the latest value
                reach_i = belief_distr->belief_[combo_index] * pow(10, this_node->reach_adjustment[b]);
                while (true) {
                    if (reach_i < pow(10, 13)) { // goes down two steps
                        break;
                    }
                    this_node->reach_adjustment[b] -= 2;
                    reach_i *= 0.01;
                    for (int a = 0; a < a_max; a++) {
                        strategy_->ulong_wavg_[rnb0 + a] *= 0.01;
                    }
                }
                pthread_mutex_unlock(&this_node->mutex_);
            }

            // final update
            reach_i = belief_distr->belief_[combo_index] * pow(10, this_node->reach_adjustment[b]);
            if (reach_i > pow(10, 13 + r)) {
                logger::critical(
                        "too large!! round %d is reach = %.16f | adjusted = %.16f | reach_adjustment[%d] = %d",
                        r,
                        belief_distr->belief_[combo_index],
                        reach_i,
                        b,
                        this_node->reach_adjustment[b]
                );
            }

            for (int a = 0; a < a_max; a++) {
                double new_wavg = strategy_->ulong_wavg_[rnb0 + a] + (reach_i * distr_rnb[a]);
                strategy_->ulong_wavg_[rnb0 + a] = (ULONG_WAVG) new_wavg;
            }
        }

        /*
         * belief_distr update + pruning
         */
        for (int a = 0; a < a_max; a++) {
            double action_prob = distr_rnb[a];

            // check pruning if action_prob = 0, skipping river node and terminal node
            if (iter_prune_flag
                && action_prob == 0.0
                && !this_node->children[a]->IsTerminal()
                && this_node->children[a]->GetRound() != HOLDEM_ROUND_RIVER) {
                auto regret = strategy_->double_regret_[rnb0 + a];
                if (regret <= cfr_param_->rollout_prune_thres) {
                    // prune it and continue
                    child_ranges[a]->Prune(combo_index);
                    continue;
                }
            }

            // zero all extremelly small values which <0.03
            if (action_prob < RANGE_ROLLOUT_PRUNE_THRESHOLD) {
                child_ranges[a]->Zero(combo_index);
            } else {
                // NOTE(kwok): child_ranges are all the same prior to this calculation
                // NOTE(kwok): A child node's reach probability equals:
                //
                //                   [how much we believe the opponent holding `combo_index`]
                //                                            ×
                // [how possible we will choose the action `a` in the case of the opponent holding `combo_index`]
                child_ranges[a]->belief_[combo_index] *= action_prob;
            }

            // safe-guarding
            auto new_v = child_ranges[a]->belief_[combo_index];
            if (new_v > 0.0 && new_v < pow(10, -14)) {
                logger::warn("🚨reach is too small %.16f [r %d]", new_v, this_node->GetRound());
            }
        }
    }
}

/**
 *
 * Reminder: if the next node is the root of a street, we need to handle the value properly. Skip all -1.
 *
 * @param out_belief
 * @param in_belief_map
 * @param hand_maps
 */
void VectorCfrWorker::ComputeCfu(Node *this_node,
                                 std::vector<sPrivateHandBelief *> child_reach_ranges,
                                 std::vector<sPrivateHandBelief *> child_cfus,
                                 sPrivateHandBelief *this_node_cfu,
                                 CFU_COMPUTE_MODE cfu_compute_mode,
                                 const float *all_belief_distr_1dim)
{
    int a_max = this_node->GetAmax();
    // auto r = this_node->GetRound();

    for (auto &i: priv_hand_kernel->valid_priv_hand_vector_idxes) {
        // outer pruning || lossless.
        if (this_node_cfu->IsPruned(i)) {
            continue;
        }

        // auto b = priv_hand_kernel->GetBucket(r, i);
        double final_value = 0.0;

        /*
         * compute the final value by cfu_compute_mode
         */
        switch (cfu_compute_mode) {
            case WEIGHTED_RESPONSE : {
                // multiply with strategy avg value.
                // whether the next node is the starting of a street does not matter.
                int offset = i * a_max;
                for (int a = 0; a < a_max; a++) {
                    // do it only when it is not pruned.
                    if (!child_reach_ranges[a]->IsPruned(i)) {
                        float weight = all_belief_distr_1dim[offset + a];
                        if (weight > 0.0) {
                            final_value += child_cfus[a]->belief_[i] * weight;
                        }
                    }
                }
                // safe guarding code
                if (final_value > pow(10, 15)) {
                    logger::critical("cfr too big %.16f at i = %d > 10^%d", final_value, i, 15);
                }
                // todo: if for computing the real this_node_cfu at each node, we need to weight it with reach
                break;
            }
            case SUM_RESPONSE : {
                for (int a = 0; a < a_max; a++) {
                    // do it only when it is not pruned. it should be pruned on the outlier,
                    if (!child_reach_ranges[a]->IsPruned(i)) {
                        final_value += child_cfus[a]->belief_[i];
                    }
                }
                break;
            }
            case BEST_RESPONSE : {
                // NOTE(kwok): The best response is a strategy for a player that is optimal against the opponent
                // strategy profile.
                final_value = -999999999;
                int offset = i * a_max;
                for (int a = 0; a < a_max; a++) {
                    if (all_belief_distr_1dim[offset + a] > 0) {
                        // Otherwise best response will become -1. Pruning will never be on in the best response cfu_compute_mode.
                        auto utility = child_cfus[a]->belief_[i];
                        if (utility != kBeliefPrunedFlag && utility > final_value) {
                            final_value = utility;
                        }
                    }
                }
                // todo: if computing the expl at this node.
                // weighted with weights, need to reweight with 1/ 1000, cuz it scales 1000 both on my reach and opp reach
                // final_value *= reach_ranges->ranges_[my_pos].belief_[i] / 1000.0;
                if (fabs(final_value + 999999999) < 0.001) {
                    //it means they are all pruned. probably the next node is the first node of a street. board crashes.
                    final_value = kBeliefPrunedFlag;
                }
                break;
            }
        }

        this_node_cfu->belief_[i] = final_value;
    }
}

void VectorCfrWorker::ConditionalPrune()
{
    if (cfr_param_->pruning_on) {
        iter_prune_flag = GenRndNumber(1, 100) <= cfr_param_->rollout_prune_prob * 100;
    }
}

void
VectorCfrWorker::CollectRegrets(Node *this_node,
                                std::vector<sPrivateHandBelief *> child_cfus,
                                sPrivateHandBelief *this_node_cfu)
{
    auto a_max = this_node->GetAmax();
    auto r = this_node->GetRound();
    auto n = this_node->GetN();
    auto frozen_b = this_node->frozen_b;

    // NOTE(kwok): Try all possible private hands for the opponent
    for (auto &combo_index: priv_hand_kernel->valid_priv_hand_vector_idxes) {
        // no need to update with the pruned hands, outer pruning
        if (this_node_cfu->IsPruned(combo_index)) {
            continue;
        }
        auto b = priv_hand_kernel->GetBucket(r, combo_index);
        // if frozen, no need to update regret
        if (b == frozen_b) {
            continue;
        }
        auto rnb0 = strategy_->ag_->kernel_->hash_rnba(r, n, b, 0);
        double cfu_combo_index = this_node_cfu->belief_[combo_index];
        // NOTE(kwok): Update the regrets, assuming the opponent's holding `combo_index`
        for (int a = 0; a < a_max; a++) {
            // already pruned, no need to update regret for this child node
            if (child_cfus[a]->IsPruned(combo_index)) {
                continue;
            }
            double cfu_child_a = child_cfus[a]->belief_[combo_index];
            double diff = cfu_child_a - cfu_combo_index; // the immediate regret
            double old_reg = strategy_->double_regret_[rnb0 + a];
            double new_reg = ClampRegret(old_reg, diff, cfr_param_->rm_floor);
            if (old_reg > pow(10, 15) || new_reg > pow(10, 15)) {
                logger::critical(
                        "[old reg %.16f] too big![new reg %.16f] [%.16f] [diff %.16f] [u_action %.16f] [this_node_cfu %.16f]",
                        old_reg,
                        new_reg,
                        cfr_param_->rm_floor,
                        diff,
                        cfu_child_a,
                        this_node_cfu
                );
            }
            strategy_->double_regret_[rnb0 + a] = new_reg;
        }
    }
}

std::vector<sPrivateHandBelief *> VectorCfrWorker::ExtractBeliefs(std::vector<Ranges *> &ranges, int pos)
{
    std::vector<sPrivateHandBelief *> beliefs;
    int size = ranges.size();
    beliefs.reserve(size);
    for (int a = 0; a < size; a++) {
        beliefs.emplace_back(&ranges[a]->beliefs_[pos]);
    }
    return beliefs;
}


