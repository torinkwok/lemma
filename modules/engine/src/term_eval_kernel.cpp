#include "term_eval_kernel.h"

void TermEvalKernel::Prepare(Board_t *board_ptr)
{
    board = *board_ptr;
    int rank_id = 0;

    std::set<int> unique_ranks;
    for (Card_t low = 0; low < HOLDEM_MAX_DECK - 1; low++) {
        for (Card_t high = low + 1; high < HOLDEM_MAX_DECK; high++) {
            auto priv_hand = PrivHand_t{high, low};
            if (board_ptr->PrivHandCrash(priv_hand)) {
                continue;
            }
            int rank = RankHand(high, low, board_ptr);
            auto v_idx = ToVectorIndex(high, low);
            // TODO(kwok): Let `sPrivHandRank` handle the calculation of the complete hand rank.
            sorted_infosets_by_rank[rank_id] = new sPrivHandRank{high, low, rank, v_idx};
            rank_id++;
            unique_ranks.insert(rank);
        }
    }

    if (rank_id != nCk_card(47, 2)) {
        logger::error("💢total number of hands is not correct = %d. We're expecting nCk_Card(47, 2)", rank_id);
    }

    // allocate memory on the heap
    n_unique_rank = unique_ranks.size();
    equal_ranks_first_idxs = new int[n_unique_rank];
    weaker_ranks_first_idxs = new int[n_unique_rank];

    Sort();
    PreStack();
}

void TermEvalKernel::Sort()
{
    std::sort(sorted_infosets_by_rank.begin(), sorted_infosets_by_rank.end(),
              [](const sPrivHandRank *lhs, const sPrivHandRank *rhs)
              {
                  return lhs->RankHighLowSort(rhs);
              }
    );

    // for (auto *a: sorted_infosets_by_rank) {
    //     a->Print();
    // }

    min_rank = sorted_infosets_by_rank[0]->rank;

    // NOTE(kwok): build an inverse lookup
    for (unsigned long rank_id = 0; rank_id < sorted_infosets_by_rank.size(); rank_id++) {
        auto h = sorted_infosets_by_rank[rank_id]->high_card;
        auto l = sorted_infosets_by_rank[rank_id]->low_card;
        rank_idxs_by_high_low[h][l] = rank_id;
        // Cardset c = emptyCardset();
        // addCardToCardset(&c, suitOfCard(list_[rank_id]->high_card, 4), rankOfCard(list_[rank_id]->high_card, 4));
        // addCardToCardset(&c, suitOfCard(list_[rank_id]->low_card, 4), rankOfCard(list_[rank_id]->low_card, 4));
        // printf("%s pos:%d\n", CardsToString(c.cards).c_str(), rank_id);
    }
}

void TermEvalKernel::PreStack()
{
    int cursor_rank = min_rank;
    int equal_ranks_first_i = 0;

    equal_ranks_first_idxs[equal_ranks_first_i] = 0;
    for (int i = 0; i < HOLDEM_MAX_HANDS_PERMUTATION_EXCLUDE_BOARD; i++) {
        auto priv_hand_rank = sorted_infosets_by_rank[i];
        auto rank = priv_hand_rank->rank;
        if (rank != cursor_rank) {
            cursor_rank = rank;
            equal_ranks_first_i++;
            equal_ranks_first_idxs[equal_ranks_first_i] = i;
        }
    }

    weaker_ranks_first_idxs[n_unique_rank - 1] = HOLDEM_MAX_HANDS_PERMUTATION_EXCLUDE_BOARD;
    for (int i = 0; i < n_unique_rank - 1; i++) {
        weaker_ranks_first_idxs[i] = equal_ranks_first_idxs[i + 1];
    }
}

// FIXME(kwok): The number of players is not supposed to be fixed to 2.
void TermEvalKernel::FastShowdownEval(const double *opp_full_belief,
                                      double *io_my_full_belief,
                                      int spent)
{
    // NOTE(kwok): transform the number of the opponents' beliefs from 1326 into 1081 accordingly
    double opp_belief_of_sorted_ranks[HOLDEM_MAX_HANDS_PERMUTATION_EXCLUDE_BOARD];
    for (auto i = 0; i < HOLDEM_MAX_HANDS_PERMUTATION_EXCLUDE_BOARD; i++) {
        opp_belief_of_sorted_ranks[i] = opp_full_belief[sorted_infosets_by_rank[i]->vector_idx];
    }

    double opp_nets_by_rank[n_unique_rank];
    double opp_nets_by_combo[HOLDEM_MAX_DECK * HOLDEM_MAX_DECK];
    for (auto &net: opp_nets_by_combo) net = 0;

    int skipping_ranks_of_combo[46 * HOLDEM_MAX_DECK]; // 46 = 52 - 5 - 1
    for (auto &rank_id: skipping_ranks_of_combo) rank_id = -1; // -1 as a cutoff flag

    StackOppShowdownProb(opp_belief_of_sorted_ranks, opp_nets_by_rank, opp_nets_by_combo, skipping_ranks_of_combo);

    int last_skipping_ranks_for_card[HOLDEM_MAX_DECK];
    for (auto &rank_id: last_skipping_ranks_for_card) rank_id = 0;

    // compute the drift
    double last_nets_for_card[HOLDEM_MAX_DECK];
    for (auto c = 0; c < HOLDEM_MAX_DECK; c++) {
        last_nets_for_card[c] = opp_nets_by_combo[ComboIdx(0, c)];
    }

    // NOTE(kwok): evaluate our full belief against the opponents'
    for (int rank_id = 0; rank_id < n_unique_rank; rank_id++) {
        double opp_rank_net = opp_nets_by_rank[rank_id];
        // NOTE(kwok): enumerate all the private hands ranked same
        for (int i = equal_ranks_first_idxs[rank_id]; i < weaker_ranks_first_idxs[rank_id]; i++) {
            auto v_idx = sorted_infosets_by_rank[i]->vector_idx;
            if (io_my_full_belief[v_idx] == kBeliefPrunedFlag) {
                continue;
            }
            double my_total_drift = 0.0;
            for (auto &c: sorted_infosets_by_rank[i]->GetHandPair()) {
                auto combo = ComboIdx(last_skipping_ranks_for_card[c], c);
                if (skipping_ranks_of_combo[combo] != rank_id) {
                    last_skipping_ranks_for_card[c]++;
                    last_nets_for_card[c] = opp_nets_by_combo[ComboIdx(last_skipping_ranks_for_card[c], c)];
                }
                my_total_drift += last_nets_for_card[c]; // no need to delete double count [high, low] as it must be the same value
            }
            double my_net = opp_rank_net - my_total_drift;
            io_my_full_belief[v_idx] = my_net * spent;
        }
    }
}

void TermEvalKernel::NaiveShowdownEval(const double *opp_belief,
                                       double *io_my_full_belief,
                                       int spent)
{
    auto begin = std::chrono::steady_clock::now();

    for (auto my_pos = 0; my_pos < HOLDEM_MAX_HANDS_PERMUTATION_EXCLUDE_BOARD; my_pos++) {
        auto v_idx = sorted_infosets_by_rank[my_pos]->vector_idx;
        // pruning
        if (io_my_full_belief[v_idx] == kBeliefPrunedFlag) {
            continue;
        }
        auto weighted_sum = 0.0;
        for (auto opp_pos = 0; opp_pos < HOLDEM_MAX_HANDS_PERMUTATION_EXCLUDE_BOARD; opp_pos++) {
            if (my_pos == opp_pos) {
                continue;
            }
            auto weight = opp_belief[sorted_infosets_by_rank[opp_pos]->vector_idx];
            // no hand belief or just too low
            if (weight == 0) {
                continue;
            }
            if (sorted_infosets_by_rank[my_pos]->RankEqual(sorted_infosets_by_rank[opp_pos])) {
                continue;
            }
            if (sorted_infosets_by_rank[my_pos]->CardCrash(sorted_infosets_by_rank[opp_pos])) {
                continue;
            }
            // compare
            weighted_sum += my_pos > opp_pos ? weight : -1.0 * weight;
        }
        io_my_full_belief[v_idx] = weighted_sum * spent;
    }

    auto end = std::chrono::steady_clock::now();
    auto total_time_ms = std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count();
    logger::trace("naive showdown || term node eval takes = %d  (microseconds)", total_time_ms);
}

/**
 * not board-card safe
 * @param opp_full_belief
 * @param io_my_full_belief
 * @param spent
 * @param win_loss_multiplier win = 1, lose = -1
 */
void TermEvalKernel::FastFoldEval(const double *opp_full_belief,
                                  double *io_my_full_belief,
                                  int spent)
{
    double folding_drift[HOLDEM_MAX_DECK];
    for (auto &i: folding_drift) {
        i = 0;
    }
    double base = 0.0;
    StackFoldingProb(opp_full_belief, folding_drift, base);
    for (auto my_pos = 0; my_pos < HOLDEM_MAX_HANDS_PERMUTATION_EXCLUDE_BOARD; my_pos++) {
        auto v_idx = sorted_infosets_by_rank[my_pos]->vector_idx;
        // pruning
        if (io_my_full_belief[v_idx] == kBeliefPrunedFlag) {
            continue;
        }
        auto high_low = FromVectorIndex(v_idx);
        auto total_drift =
                folding_drift[high_low.first]
                + folding_drift[high_low.second] - opp_full_belief[v_idx]; // - double count
        auto net = base - total_drift;
        io_my_full_belief[v_idx] = net * spent;
    }
}

void TermEvalKernel::NaiveFoldEval(const double *opp_full_belief,
                                   double *io_my_belief,
                                   int spent)
{
    auto begin = std::chrono::steady_clock::now();
    for (auto my_pos = 0; my_pos < FULL_HAND_BELIEF_SIZE; my_pos++) {
        // pruning
        if (io_my_belief[my_pos] == kBeliefPrunedFlag) {
            continue;
        }
        auto weighted_sum = 0.0;
        for (auto opp_pos = 0; opp_pos < FULL_HAND_BELIEF_SIZE; opp_pos++) {
            if (my_pos == opp_pos) {
                continue;
            }
            //check if i and j clashes
            if (VectorIdxClash(my_pos, opp_pos)) {
                continue;
            }
            auto weight = opp_full_belief[opp_pos];
            //no hand belief or just too low
            if (weight > 0.0) {
                weighted_sum += weight;
            }
        }
        io_my_belief[my_pos] = weighted_sum * spent;
    }
    auto end = std::chrono::steady_clock::now();
    auto total_time_ms = std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count();
    logger::trace("naive fold || term node eval takes = %d  (microseconds)", total_time_ms);

}

void TermEvalKernel::StackFoldingProb(const double *opp_full_belief,
                                      double *drift_by_card, double &sum)
{
    for (int i = 0; i < FULL_HAND_BELIEF_SIZE; i++) {
        double w = opp_full_belief[i];
        //skipping 0 and -1
        if (w > 0.0) {
            sum += w;
            auto high_low = FromVectorIndex(i);
            drift_by_card[high_low.first] += w;
            drift_by_card[high_low.second] += w;
        }
    }
}

/**
 * the net sum is calculated iteratively...
 * using a skipping linked list.
 * @param opp_belief_of_sorted_ranks
 * @param io_opp_nets_by_rank
 * @param io_opp_nets_by_combo
 * @param io_skipping_ranks_of_combo , default -1
 */
void TermEvalKernel::StackOppShowdownProb(const double *opp_belief_of_sorted_ranks,
                                          double *io_opp_nets_by_rank,
                                          double *io_opp_nets_by_combo,
                                          int *io_skipping_ranks_of_combo)
{
    double opp_belief_total_sum = 0.0; // used for the computation of `io_opp_nets_by_rank`
    double opp_belief_sums_by_rank[n_unique_rank]; // used for the computation of `io_opp_nets_by_rank`
    for (auto &i: opp_belief_sums_by_rank) {
        i = 0;
    }

    // temporary object for constructing io_skipping_ranks_of_combo
    // e.g. 0->4, 1->67, 2->89, 3->126... skip_list_idx -> rank
    int last_skipping_ranks_by_card[HOLDEM_MAX_DECK];
    for (auto &rank_id: last_skipping_ranks_by_card) {
        rank_id = -1;
    }

    // skipping linked list no more than 52 - 5 - 1 = 46 unique rank.  rank * card
    // in practise it is a lot less
    double opp_belief_sums_by_rank_card[46][HOLDEM_MAX_DECK];
    for (auto &sums: opp_belief_sums_by_rank_card) {
        for (double &s: sums) {
            s = 0;
        }
    }

    double opp_belief_sums_by_card[HOLDEM_MAX_DECK];
    for (auto &s: opp_belief_sums_by_card) {
        s = 0;
    }

    /* NOTE(kwok): The First Iteration */

    for (int rank_id = 0; rank_id < n_unique_rank; rank_id++) {
        for (int i = equal_ranks_first_idxs[rank_id]; i < weaker_ranks_first_idxs[rank_id]; i++) {
            double w = opp_belief_of_sorted_ranks[i];
            // WARNING: FIXME(kwok): the below code snippet makes showdown not able to handle 0 items in opp belief!
            // Maybe because it makes the skipping rank_list part not accurate, i.e. missing some steps.
            // if (w == 0) {
            //     continue;
            // }
            opp_belief_sums_by_rank[rank_id] += w;
            // NOTE(kwok): deal with nets in terms of individual card
            for (auto &c: sorted_infosets_by_rank[i]->GetHandPair()) {
                int combo = ComboIdx(last_skipping_ranks_by_card[c], c);
                if (combo < 0 || io_skipping_ranks_of_combo[combo] != rank_id) {
                    // NOTE(kwok): `combo` being negative means that `last_skipping_ranks_by_card[c]` has
                    // yet to be initialized (still being -1). Self-increase to zero it.
                    last_skipping_ranks_by_card[c]++;
                    // NOTE(kwok): `combo` remains negative here. `combo + HOLDEM_MAX_DECK` here is equivalent
                    // to `ComboIdx(last_skipping_ranks_by_card[c], c)`.
                    io_skipping_ranks_of_combo[combo + HOLDEM_MAX_DECK] = rank_id;
                }
                opp_belief_sums_by_rank_card[last_skipping_ranks_by_card[c]][c] += w;
                opp_belief_sums_by_card[c] += w;
            }
        }
        opp_belief_total_sum += opp_belief_sums_by_rank[rank_id];
    }

    /* NOTE(kwok): The Second Iteration. Dynamic Programming Structures. */

    // NOTE(kwok): compute the opponents' nets by rank
    for (int rank_id = 0; rank_id < n_unique_rank; rank_id++) {
        if (rank_id == 0) {
            // NOTE(kwok): If `rank_id` is 0, hands bucketed by it will always lose to hands with
            // other rank ids.
            io_opp_nets_by_rank[rank_id] = opp_belief_sums_by_rank[rank_id] - opp_belief_total_sum;
            continue;
        }
        io_opp_nets_by_rank[rank_id] =
                // NOTE(kwok): take advantage of the value previously computed for the last rank id
                io_opp_nets_by_rank[rank_id - 1]
                // This term is equivalent to the previous iteration's `opp_belief_sums_by_rank[rank_id]`
                // FIXME(kwok): Is this term repeatly included?
                + opp_belief_sums_by_rank[rank_id - 1]
                + opp_belief_sums_by_rank[rank_id];
    }

    // NOTE(kwok): compute the opponents' nets by combo
    for (auto c = 0; c < HOLDEM_MAX_DECK; c++) {
        // ignoring non-legit card
        if (board.CardCrash(c)) {
            continue;
        }
        int rank_id = 0; // the last skipping index
        while (true) { // TODO(kwok): Refactor this while-loop into a standard for-loop.
            if (rank_id == 0) {
                // NOTE(kwok): If `rank_id` is 0, hands with the same rank which contains `c` will always
                // lose to other hands containing the same `c`.
                io_opp_nets_by_combo[ComboIdx(rank_id, c)] =
                        opp_belief_sums_by_rank_card[rank_id][c] - opp_belief_sums_by_card[c];
                rank_id++;
                continue;
            }
            int combo = ComboIdx(rank_id, c);
            // cut off if -1
            if (io_skipping_ranks_of_combo[combo] == -1) {
                break;
            }
            io_opp_nets_by_combo[combo] =
                    // NOTE(kwok): Take advantage of the value previously computed for the last rank id.
                    // This `combo - HOLDEM_MAX_DECK` is equivalent to `ComboIdx(rank_id - 1, c)`.
                    io_opp_nets_by_combo[combo - HOLDEM_MAX_DECK /* combo with the last rank_id */]
                    // This term is equivalent to the previous iteration's `opp_belief_sums_by_rank[rank_id]`
                    // FIXME(kwok): Is this term repeatly included?
                    + opp_belief_sums_by_rank_card[rank_id - 1][c]
                    + opp_belief_sums_by_rank_card[rank_id][c];
            rank_id++;
        }
    }
}

inline int TermEvalKernel::ComboIdx(int rank_id, int card)
{
    return rank_id * HOLDEM_MAX_DECK + card;
}

// FIXME(kwok): The number of players is not supposed to be fixed to 2.
void
TermEvalKernel::FastTerminalEval(const double *opp_full_belief, double *io_my_full_belief, int spent, bool showdown)
{
    if (showdown) {
        FastShowdownEval(opp_full_belief, io_my_full_belief, spent);
    } else {
        FastFoldEval(opp_full_belief, io_my_full_belief, spent);
    }
}
