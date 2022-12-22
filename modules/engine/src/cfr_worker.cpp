#include "cfr_worker.h"
#include "node.h"

double VectorCfrWorker::Solve(Board_t board)
{
    auto *ag = strategy->ag_;
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

    // TODO(kwok): Cache the results for RIVER sub-games.
    sPrivateHandBelief local_root_belief[2]; // FIXME(kwok): The number of players is not supposed to be fixed to 2.
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
                starting_ranges.beliefs_[p].Scale(REGRET_SCALE_FACTOR);
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
            // walk down training tree alternatively, sampling all nodes
            for (int trainee = 0; trainee < active_players; trainee++) {
                // FIXME(kwok): The number of players is not supposed to be fixed to 2.
                auto opp_local_root_belief = local_root_belief[1 - trainee];
                opp_local_root_belief.Scale(REGRET_SCALE_FACTOR);
                sPrivateHandBelief *cfu_p =
                        WalkTree_Alternate(ag->root_node_, trainee, &opp_local_root_belief, strategy);
                cfu_p->DotMultiply(&local_root_belief[trainee]);
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
sPrivateHandBelief *VectorCfrWorker::WalkTree_Alternate(Node *this_node,
                                                        int trainee,
                                                        sPrivateHandBelief *opp_belief,
                                                        Strategy *target_strategy,
                                                        std::optional<CFU_COMPUTE_MODE> compute_node)
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
    bool is_trainee_turn = trainee == this_node->GetActingPlayer();
    if (!compute_node.has_value()) {
        compute_node = is_trainee_turn
                       ? cfr_param_->cfu_compute_acting_playing
                       : cfr_param_->cfu_compute_opponent;
    }
    return EvalChoiceNode_Alternate(this_node, trainee, opp_belief, target_strategy, compute_node.value());
}

sPrivateHandBelief *
VectorCfrWorker::EvalChoiceNode_Alternate(Node *this_node, int trainee, sPrivateHandBelief *opp_belief,
                                          Strategy *target_strategy, CFU_COMPUTE_MODE compute_mode)
{
    auto a_max = this_node->GetAmax();
    bool is_trainee_turn = trainee == this_node->GetActingPlayer();

    /* (0) calc reach range */

    // copy the reach_ranges to child ranges
    std::vector<sPrivateHandBelief *> child_reach_ranges;
    child_reach_ranges.reserve(a_max);
    if (is_trainee_turn) {
        for (int a = 0; a < a_max; a++) {
            child_reach_ranges.emplace_back(opp_belief); // shallow copy
        }
    } else {
        for (int a = 0; a < a_max; a++) {
            child_reach_ranges.emplace_back(new sPrivateHandBelief(opp_belief)); // deep copy
        }
        CalcReachRange(this_node, opp_belief, child_reach_ranges, target_strategy);
    }

    /* (1) collect CFUs from children */

    std::vector<sPrivateHandBelief *> children_cfus;
    children_cfus.reserve(a_max);
    for (int a = 0; a < a_max; a++) {
        auto *child_node = this_node->children[a];
        auto *child_node_reach_range = child_reach_ranges[a]; // will turn into `opp_belief` for child node
        children_cfus.emplace_back(
                // FIXME(kwok): Be careful that do not pass `compute_mode` as is
                WalkTree_Alternate(child_node, trainee, child_node_reach_range, target_strategy)
        );
    }

    /* (2) assemble CFU of the current node */

    // precompute strategy only if it is trainee's turn. range rollout is fine, which is computed on the fly
    float *all_belief_distr_1dim = nullptr;
    if (is_trainee_turn) {
        // unless we have pruning, we don't need to skip any one. it is a bit wasteful but makes it more accurate.
        all_belief_distr_1dim = new float[FULL_HAND_BELIEF_SIZE * a_max];
        for (auto &priv_hand_idx: priv_hand_kernel->valid_priv_hand_vector_idxes) {
            int offset = a_max * priv_hand_idx;
            auto b = priv_hand_kernel->GetBucket(this_node->GetRound(), priv_hand_idx);
            target_strategy->ComputeStrategy(this_node, b, all_belief_distr_1dim + offset,
                                             cfr_param_->strategy_cal_mode_
            );
        }
    } else {
        // FIXME(kwok): `all_belief_distr_1dim` may point to unallocated memory.
    }

    // NOTE(kwok): A utility of +1 is given for a win, and −1 for a loss.
    auto this_node_cfu = new sPrivateHandBelief(0.0);
    // FIXME(kwok): `all_belief_distr_1dim` may point to unallocated memory.
    ComputeCfu(this_node, child_reach_ranges, children_cfus, this_node_cfu, compute_mode, all_belief_distr_1dim);

    /* (3) learn from regrets */

    // only trainee learns from regrets
    if (cfr_param_->regret_learning_on) {
        if (is_trainee_turn) {
            CollectRegrets(this_node, children_cfus, this_node_cfu, target_strategy);
        }
    }
    // delete child pop up this_node_cfu
    for (int a = 0; a < a_max; a++) {
        delete children_cfus[a];
        if (!is_trainee_turn) {
            delete child_reach_ranges[a];
        }
    }
    if (is_trainee_turn) {
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
    CalcReachRange(this_node, &reach_ranges->beliefs_[actor], actor_child_beliefs, strategy);

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
        CollectRegrets(this_node, actor_child_cfu, &cfu->beliefs_[actor], strategy);
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

/// For both sides, conduct
///     * belief updating and pruning
///     * WAVG updating
void VectorCfrWorker::CalcReachRange(Node *this_node, sPrivateHandBelief *belief,
                                     std::vector<sPrivateHandBelief *> &child_ranges,
                                     Strategy *target_strategy)
{
    auto r = this_node->GetRound();
    auto a_max = this_node->GetAmax();
    auto sibling_i = this_node->GetN();
    auto frozen_b = this_node->frozen_b;
    for (auto &priv_hand_i: priv_hand_kernel->valid_priv_hand_vector_idxes) {
        if (belief->IsPruned(priv_hand_i)) {
            continue;
        }
        if (belief->IsZero(priv_hand_i)) {
            continue;
        }
        auto b = priv_hand_kernel->GetBucket(r, priv_hand_i);
        auto rnb0 = target_strategy->ag_->kernel_->hash_rnba(r, sibling_i, b, 0);
        float distr_rnb[a_max]; // strategy if private hand is priv_hand_i
        target_strategy->ComputeStrategy(this_node, b, distr_rnb, cfr_param_->strategy_cal_mode_);
        /* update WAVG, if not frozen */
        if (cfr_param_->rm_avg_update == AVG_CLASSIC && b != frozen_b) {
            double reach_i = belief->belief_[priv_hand_i] * pow(10, this_node->reach_adjustment[b]);
            // adjustment adaptively goes up
            if (reach_i > 0.0 && reach_i < 100) {
                pthread_mutex_lock(&this_node->mutex_);
                reach_i = belief->belief_[priv_hand_i] * pow(10, this_node->reach_adjustment[b]);
                while (true) {
                    if (reach_i > 100) {
                        break;
                    }
                    this_node->reach_adjustment[b] += 1;
                    reach_i *= 10;
                    for (int a = 0; a < a_max; a++) {
#ifdef DEBUG_EAGER_LOOKUP
                        // TODO(kwok): 🦊
                        target_strategy->eager_ulong_wavg_[rnb0 + a] *= 10;
#endif
                        // TODO(kwok): 🐦
                        target_strategy->ulong_wavg_->upsert(rnb0 + a, [](auto &n) { n *= 10; });
                        // this_node->ulong_wavg_[this_node->HashBa(b, a)] *= 10;
                    }
                }
                pthread_mutex_unlock(&this_node->mutex_);
            }
            // adjustment adaptively goes down
            if (reach_i > pow(10, 15)) {
                pthread_mutex_lock(&this_node->mutex_);
                // get the latest value
                reach_i = belief->belief_[priv_hand_i] * pow(10, this_node->reach_adjustment[b]);
                while (true) {
                    if (reach_i < pow(10, 13)) { // goes down two steps
                        break;
                    }
                    this_node->reach_adjustment[b] -= 2;
                    reach_i *= 0.01;
                    for (int a = 0; a < a_max; a++) {
#ifdef DEBUG_EAGER_LOOKUP
                        // TODO(kwok): 🦊
                        target_strategy->eager_ulong_wavg_[rnb0 + a] *= 0.01;
#endif
                        // TODO(kwok): 🐦
                        target_strategy->ulong_wavg_->upsert(rnb0 + a, [](auto &n) { n *= 0.01; });
#ifdef DEBUG_EAGER_LOOKUP
                        if (target_strategy->eager_ulong_wavg_[rnb0 + a] != target_strategy->ulong_wavg_->find(rnb0 + a)) {
                            logger::warn("🦊%ul vs. 🐦%ul",
                                             target_strategy->eager_ulong_wavg_[rnb0 + a],
                                             target_strategy->ulong_wavg_->find(rnb0 + a)
                            );
                        }
#endif
                    }
                }
                pthread_mutex_unlock(&this_node->mutex_);
            }
            // final update
            reach_i = belief->belief_[priv_hand_i] * pow(10, this_node->reach_adjustment[b]);
            if (reach_i > pow(10, 13 + r)) {
                logger::critical(
                        "too large!! round %d is reach = %.16f | adjusted = %.16f | reach_adjustment[%d] = %d",
                        r,
                        belief->belief_[priv_hand_i],
                        reach_i,
                        b,
                        this_node->reach_adjustment[b]
                );
            }
            for (int a = 0; a < a_max; a++) {
                double accu = reach_i * distr_rnb[a];
#ifdef DEBUG_EAGER_LOOKUP
                // TODO(kwok): 🦊
                double new_wavg = target_strategy->eager_ulong_wavg_[rnb0 + a] + accu;
                target_strategy->eager_ulong_wavg_[rnb0 + a] = (ULONG_WAVG) new_wavg;
#endif
                // TODO(kwok): 🐦
                target_strategy->ulong_wavg_->upsert(rnb0 + a, [&](auto &n) { n += accu; }, accu);
#ifdef DEBUG_EAGER_LOOKUP
                if (target_strategy->eager_ulong_wavg_[rnb0 + a] != target_strategy->ulong_wavg_->find(rnb0 + a)) {
                    logger::warn("🦊%lu vs. 🐦%lu vs. %lu",
                                 target_strategy->eager_ulong_wavg_[rnb0 + a],
                                 target_strategy->ulong_wavg_->find(rnb0 + a),
                                 (ULONG_WAVG) new_wavg
                    );
                }
#endif
            }
        }
        /* belief update + pruning */
        for (int a = 0; a < a_max; a++) {
            double action_prob = distr_rnb[a];
            // check pruning if action_prob = 0, skipping river node and terminal node
            if (iter_prune_flag
                && action_prob == 0.0
                && !this_node->children[a]->IsTerminal()
                && this_node->children[a]->GetRound() != HOLDEM_ROUND_RIVER) {
#ifdef DEBUG_EAGER_LOOKUP
                // TODO(kwok): 🦊
                bool is_fox_pruned = false;
                auto regret = target_strategy->eager_double_regret_[rnb0 + a];
                if (regret <= cfr_param_->rollout_prune_thres) {
                    // prune it and continue
                    child_ranges[a]->Prune(priv_hand_i);
                    is_fox_pruned = true;
                }
#endif
                // TODO(kwok): 🐦
                bool is_dove_pruned = false;
                target_strategy->double_regret_->insert(rnb0 + a);
                target_strategy->double_regret_->find_fn(
                        rnb0 + a,
                        [&](const auto &regret)
                        {
                            if (regret <= cfr_param_->rollout_prune_thres) {
                                // prune it and continue
                                child_ranges[a]->Prune(priv_hand_i);
                                is_dove_pruned = true;
                            }
                        }
                );

                if (is_dove_pruned
#ifdef DEBUG_EAGER_LOOKUP
                    && is_fox_pruned
#endif
                        ) {
                    continue;
                }
#ifdef DEBUG_EAGER_LOOKUP
                else if (is_dove_pruned ^ is_fox_pruned) {
                    logger::critical("🐦%d vs. 🦊%d", is_dove_pruned, is_fox_pruned);
                }
#endif
            }
            // zero all trivial values less than 0.03
            if (action_prob < RANGE_ROLLOUT_PRUNE_THRESHOLD) {
                child_ranges[a]->Zero(priv_hand_i);
            } else {
                // NOTE(kwok): child_ranges are all the same prior to this calculation
                // NOTE(kwok): A child node's reach probability equals:
                //
                //                   [how much we believe the opponent holding `priv_hand_i`]
                //                                            ×
                // [how possible we will choose the action `a` in the case of the opponent holding `priv_hand_i`]
                child_ranges[a]->belief_[priv_hand_i] *= action_prob;
            }
            // safe-guarding
            auto new_v = child_ranges[a]->belief_[priv_hand_i];
            if (new_v > 0.0 && new_v < pow(10, -14)) {
                logger::warn("🚨reach is too small %.16f [r %d]", new_v, this_node->GetRound());
            }
        }
    }
}

/// Reminder: If the next node is the root of a street, we get to handle the value properly. Skip all -1s.
void VectorCfrWorker::ComputeCfu(Node *this_node,
                                 std::vector<sPrivateHandBelief *> child_reach_ranges,
                                 std::vector<sPrivateHandBelief *> child_cfus,
                                 sPrivateHandBelief *this_node_cfu,
                                 CFU_COMPUTE_MODE cfu_compute_mode,
                                 const float *all_belief_distr_1dim) const
{
    int a_max = this_node->GetAmax();
    for (auto &priv_hand_idx: priv_hand_kernel->valid_priv_hand_vector_idxes) {
        // outer pruning || lossless
        if (this_node_cfu->IsPruned(priv_hand_idx)) {
            continue;
        }
        double expected_utility = 0.0;
        switch (cfu_compute_mode) {
            case WEIGHTED_RESPONSE : {
                // whether the next node is the root of a street does not matter
                int offset = priv_hand_idx * a_max;
                for (int a = 0; a < a_max; a++) {
                    // proceed only when candidate node has not been pruned
                    if (child_reach_ranges[a]->IsPruned(priv_hand_idx)) {
                        continue;
                    }
                    // FIXME(kwok): `all_belief_distr_1dim` may point to unallocated memory.
                    float weight = all_belief_distr_1dim[offset + a];
                    if (weight > 0.0) {
                        expected_utility += weight * child_cfus[a]->belief_[priv_hand_idx];
                    }
                }
                if (expected_utility > pow(10, 15)) {
                    logger::critical("expected utility of %.16f at priv_hand_idx=%d is considered too large (> 10^%d)",
                                     expected_utility, priv_hand_idx, 15
                    );
                }
                // todo: if for computing the real this_node_cfu at each node, we need to weight it with reach
                break;
            }
            case SUM_RESPONSE : {
                for (int a = 0; a < a_max; a++) {
                    // proceed only when candidate node has not been pruned. outliers should already
                    // been pruned
                    if (child_reach_ranges[a]->IsPruned(priv_hand_idx)) {
                        continue;
                    }
                    expected_utility += child_cfus[a]->belief_[priv_hand_idx];
                }
                break;
            }
            case BEST_RESPONSE : {
                // NOTE(kwok): The best response is a strategy for a player that is optimal against the opponent
                // strategy profile.
                expected_utility = -std::numeric_limits<double>::infinity();
                int offset = priv_hand_idx * a_max;
                for (int a = 0; a < a_max; a++) {
                    if (all_belief_distr_1dim[offset + a] > 0) {
                        // pruning will never happen in the BEST_RESPONSE mode.
                        auto utility = child_cfus[a]->belief_[priv_hand_idx];
                        if (utility != kBeliefPrunedFlag && utility > expected_utility) {
                            expected_utility = utility;
                        }
                    }
                }
                // todo: if computing the expl at this node.
                // weighted with weights, need to re-weight with 1/ 1000, cuz it scales 1000 both on my reach and opp reach
                // expected_utility *= reach_ranges->ranges_[my_pos].belief_[priv_hand_idx] / 1000.0;
                if (fabs(expected_utility + 999999999) < 0.001) {
                    // it means they are all pruned. probably the next node is the first node of a street. board crashes.
                    expected_utility = kBeliefPrunedFlag;
                }
                break;
            }
        }
        this_node_cfu->belief_[priv_hand_idx] = expected_utility;
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
                                std::vector<sPrivateHandBelief *> child_cfus, sPrivateHandBelief *this_node_cfu,
                                Strategy *target_strategy)
{
    auto a_max = this_node->GetAmax();
    auto r = this_node->GetRound();
    auto n = this_node->GetN();
    auto frozen_b = this_node->frozen_b;

    // iterate all possible private hands within the current infoset for `this_node`
    for (auto &priv_hand_i: priv_hand_kernel->valid_priv_hand_vector_idxes) {
        // no need to update pruned hands
        if (this_node_cfu->IsPruned(priv_hand_i)) {
            continue;
        }
        auto b_priv_hand_i = priv_hand_kernel->GetBucket(r, priv_hand_i);
        // if frozen, no need to update regret
        if (b_priv_hand_i == frozen_b) {
            continue;
        }
        auto rnb0_priv_hand_i = target_strategy->ag_->kernel_->hash_rnba(r, n, b_priv_hand_i, 0);
        double cfu_priv_hand_i = this_node_cfu->belief_[priv_hand_i];
        // NOTE(kwok): Update the regrets, assuming the opponent's holding `priv_hand_i`
        for (int a = 0; a < a_max; a++) {
            // already pruned, no need to update regret for this child node
            if (child_cfus[a]->IsPruned(priv_hand_i)) {
                continue;
            }
            double cfu_child_a = child_cfus[a]->belief_[priv_hand_i];
            double diff = cfu_child_a - cfu_priv_hand_i; // naive regret
#ifdef DEBUG_EAGER_LOOKUP
            // TODO(kwok): 🦊
            double old_reg = strategy_->eager_double_regret_[rnb0_priv_hand_i + a];
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
            strategy_->eager_double_regret_[rnb0_priv_hand_i + a] = new_reg;
#endif
            // TODO(kwok): 🐦
            target_strategy->double_regret_->insert(rnb0_priv_hand_i + a);
            target_strategy->double_regret_->update_fn(
                    rnb0_priv_hand_i + a,
                    [&](auto &regret)
                    {
                        double old_reg = regret;
                        double new_reg = ClampRegret(old_reg, diff, cfr_param_->rm_floor);
                        // the total regret should have a ceiling
                        if (old_reg > pow(10, 15) || new_reg > pow(10, 15)) {
                            logger::critical(
                                    "[old_reg %.16f] is considered too large: [new_reg %.16f] [cfr_param_->rm_floor %.16f] "
                                    "[diff %.16f] [cfu_child_a %.16f] [this_node_cfu %.16f]",
                                    old_reg, new_reg, cfr_param_->rm_floor,
                                    diff, cfu_child_a, this_node_cfu
                            );
                        }
                        regret = new_reg;
                    }
            );
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
