#include "term_eval_kernel.h"

void TermEvalKernel::Prepare(Board_t *board_ptr)
{
    board = *board_ptr;
    int rank_i = 0;

    std::set<int> unique_ranks;
    for (Card_t low = 0; low < HOLDEM_MAX_DECK - 1; low++) {
        for (Card_t high = low + 1; high < HOLDEM_MAX_DECK; high++) {
            auto priv_hand = PrivHand_t{high, low};
            if (board_ptr->PrivHandCrash(priv_hand)) {
                continue;
            }
            int rank = RankHand(high, low, board_ptr);
            auto vector_idx = ToVectorIndex(high, low);
            // TODO(kwok): Let `sPrivHandRank` handle the calculation of the complete hand rank.
            sorted_infoset_by_rank[rank_i] = new sPrivHandRank{high, low, rank, vector_idx};
            rank_i++;
            unique_ranks.insert(rank);
        }
    }

    if (rank_i != nCk_card(47, 2)) {
        logger::error("💢total number of hands is not correct = %d. We're expecting nCk_Card(47, 2)", rank_i);
    }

    // allocate memory on the heap
    n_unique_rank = unique_ranks.size();
    rank_first_equal_index = new int[n_unique_rank];
    rank_first_losing_index = new int[n_unique_rank];

    Sort();
    PreStack();
}

void TermEvalKernel::Sort()
{
    std::sort(sorted_infoset_by_rank.begin(), sorted_infoset_by_rank.end(),
              [](const sPrivHandRank *lhs, const sPrivHandRank *rhs)
              {
                  return lhs->RankHighLowSort(rhs);
              }
    );

    // for (auto *a: sorted_infoset_by_rank) {
    //     a->Print();
    // }

    min_rank = sorted_infoset_by_rank[0]->rank;

    // NOTE(kwok): build an inverse lookup
    for (unsigned long rank_i = 0; rank_i < sorted_infoset_by_rank.size(); rank_i++) {
        auto h = sorted_infoset_by_rank[rank_i]->high_card;
        auto l = sorted_infoset_by_rank[rank_i]->low_card;
        rank_indices_by_high_low[h][l] = rank_i;
        // Cardset c = emptyCardset();
        // addCardToCardset(&c, suitOfCard(list_[rank_i]->high_card, 4), rankOfCard(list_[rank_i]->high_card, 4));
        // addCardToCardset(&c, suitOfCard(list_[rank_i]->low_card, 4), rankOfCard(list_[rank_i]->low_card, 4));
        // printf("%s pos:%d\n", CardsToString(c.cards).c_str(), rank_i);
    }
}

void TermEvalKernel::PreStack()
{
    int last_seen_rank = min_rank;
    int rank_i = 0;
    rank_first_equal_index[rank_i] = 0;

    for (int i = 0; i < HOLDEM_MAX_HANDS_PERMUTATION_EXCLUDE_BOARD; i++) {
        auto priv_hand_rank = sorted_infoset_by_rank[i];
        auto rank = priv_hand_rank->rank;
        if (rank != last_seen_rank) {
            last_seen_rank = rank;
            rank_i++;
            // store new rank's value and index
            rank_first_equal_index[rank_i] = i;
        }
    }

    for (int j = 0; j < n_unique_rank - 1; j++) {
        rank_first_losing_index[j] = rank_first_equal_index[j + 1];
    }
    rank_first_losing_index[n_unique_rank - 1] = HOLDEM_MAX_HANDS_PERMUTATION_EXCLUDE_BOARD;
}

// FIXME(kwok): The number of players is not supposed to be fixed to 2.
void TermEvalKernel::FastShowdownEval(const double *opp_full_belief,
                                      double *io_my_full_belief,
                                      int spent)
{
    // NOTE(kwok): transform the number of the opponents' beliefs from 1326 into 1081 accordingly
    double sorted_opp_belief_of_ranks[HOLDEM_MAX_HANDS_PERMUTATION_EXCLUDE_BOARD];
    for (auto i = 0; i < HOLDEM_MAX_HANDS_PERMUTATION_EXCLUDE_BOARD; i++) {
        sorted_opp_belief_of_ranks[i] = opp_full_belief[sorted_infoset_by_rank[i]->vector_idx];
    }

    double nets_per_rank[n_unique_rank];
    double nets_per_card[HOLDEM_MAX_DECK * HOLDEM_MAX_DECK];
    for (auto &net: nets_per_card) net = 0;

    int card_skipping_ranks[46 * HOLDEM_MAX_DECK]; // 52 - 5 - 1
    for (auto &rank: card_skipping_ranks) rank = -1; // -1 as a cutoff flag

    StackShowdownProb(sorted_opp_belief_of_ranks, nets_per_rank, nets_per_card, card_skipping_ranks);

    int card_last_skip_idxes[HOLDEM_MAX_DECK];
    for (auto &i: card_last_skip_idxes) i = 0;

    // compute the drift
    double card_last_net[HOLDEM_MAX_DECK];
    for (auto c = 0; c < HOLDEM_MAX_DECK; c++) {
        card_last_net[c] = nets_per_card[ComboIdx(0, c)];
    }

    for (int rank_i = 0; rank_i < n_unique_rank; rank_i++) {
        double base = nets_per_rank[rank_i];
        for (int j = rank_first_equal_index[rank_i]; j < rank_first_losing_index[rank_i]; j++) {
            auto v_idx = sorted_infoset_by_rank[j]->vector_idx;
            if (io_my_full_belief[v_idx] == kBeliefPrunedFlag) {
                continue;
            }
            double total_drift = 0.0;
            for (auto &c: sorted_infoset_by_rank[j]->GetHandPair()) {
                auto idx = ComboIdx(card_last_skip_idxes[c], c);
                if (card_skipping_ranks[idx] != rank_i) {
                    card_last_skip_idxes[c]++;
                    card_last_net[c] = nets_per_card[ComboIdx(card_last_skip_idxes[c], c)];
                }
                total_drift += card_last_net[c]; // No need to delete double count [high, low] as it must be the same value
            }
            double net = base - total_drift;
            io_my_full_belief[v_idx] = net * spent;
        }
    }
}

void TermEvalKernel::NaiveShowdownEval(const double *opp_belief,
                                       double *io_my_full_belief,
                                       int spent)
{
    auto begin = std::chrono::steady_clock::now();

    for (auto my_pos = 0; my_pos < HOLDEM_MAX_HANDS_PERMUTATION_EXCLUDE_BOARD; my_pos++) {
        auto v_idx = sorted_infoset_by_rank[my_pos]->vector_idx;
        // pruning
        if (io_my_full_belief[v_idx] == kBeliefPrunedFlag) {
            continue;
        }
        auto weighted_sum = 0.0;
        for (auto opp_pos = 0; opp_pos < HOLDEM_MAX_HANDS_PERMUTATION_EXCLUDE_BOARD; opp_pos++) {
            if (my_pos == opp_pos) {
                continue;
            }
            auto weight = opp_belief[sorted_infoset_by_rank[opp_pos]->vector_idx];
            // no hand belief or just too low
            if (weight == 0) {
                continue;
            }
            if (sorted_infoset_by_rank[my_pos]->RankEqual(sorted_infoset_by_rank[opp_pos])) {
                continue;
            }
            if (sorted_infoset_by_rank[my_pos]->CardCrash(sorted_infoset_by_rank[opp_pos])) {
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
        auto v_idx = sorted_infoset_by_rank[my_pos]->vector_idx;
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
 * @param opp_full_belief
 * @param io_nets_per_rank
 * @param io_nets_per_card
 * @param io_card_skipping_ranks , default -1
 */
void TermEvalKernel::StackShowdownProb(const double *opp_full_belief,
                                       double *io_nets_per_rank,
                                       double *io_nets_per_card,
                                       int *io_card_skipping_ranks)
{
    double rank_net_accu_sums[n_unique_rank];
    for (auto &i: rank_net_accu_sums) {
        i = 0;
    }

    // temporary object for constructing io_card_skipping_ranks
    // e.g. 0->4, 1->67, 2->89, 3->126... skip_list_idx -> rank
    int card_last_skipping_list_idxes[HOLDEM_MAX_DECK];
    for (auto &i: card_last_skipping_list_idxes) {
        i = -1;
    }

    // skipping linked list no more than 52 - 5 - 1 = 46 unique rank.  rank * card
    // in practise it is a lot less
    double card_rank_sum[46][HOLDEM_MAX_DECK];
    for (auto &i: card_rank_sum) {
        for (double &j: i) {
            j = 0;
        }
    }

    double card_sum[HOLDEM_MAX_DECK];
    for (auto &i: card_sum) {
        i = 0;
    }

    // first list iteration, tally the sum of belief, by rank, by card
    double opp_belief_sum = 0.0;
    for (int rank_i = 0; rank_i < n_unique_rank; rank_i++) {
        for (int j = rank_first_equal_index[rank_i]; j < rank_first_losing_index[rank_i]; j++) {
            double w = opp_full_belief[j];
            // WARNING: FIXME(kwok): the below code snippet makes showdown not able to handle 0 items in opp belief!
            // Maybe because it makes the skipping rank_list part not accurate, i.e. missing some steps.
            // if (w == 0)
            //     continue;
            rank_net_accu_sums[rank_i] += w;
            // card
            for (auto &c: sorted_infoset_by_rank[j]->GetHandPair()) {
                int idx = ComboIdx(card_last_skipping_list_idxes[c], c);
                //if idx < 0. then it is not init at all, ++ to 0
                if (idx < 0 || io_card_skipping_ranks[idx] != rank_i) {
                    card_last_skipping_list_idxes[c]++;
                    // present this one to rank_i, equivalent to ComboIdx(card_last_skipping_list_idxes[c], c);, after ++
                    io_card_skipping_ranks[idx + HOLDEM_MAX_DECK] = rank_i;
                }
                card_rank_sum[card_last_skipping_list_idxes[c]][c] += w;
                card_sum[c] += w;
            }
        }
        opp_belief_sum += rank_net_accu_sums[rank_i];
    }

    // second iter, iteratively compute needed values by rank

    // computing the base  //add the rank-1 and rank local sum
    for (int rank_i = 0; rank_i < n_unique_rank; rank_i++) {
        if (rank_i == 0) {
            io_nets_per_rank[rank_i] = rank_net_accu_sums[rank_i] - opp_belief_sum;
            continue;
        }
        io_nets_per_rank[rank_i] =
                io_nets_per_rank[rank_i - 1] + rank_net_accu_sums[rank_i - 1] + rank_net_accu_sums[rank_i];
    }

    // compute the card drift, skipping linked list.
    for (auto c = 0; c < HOLDEM_MAX_DECK; c++) {
        // skipping non-legit card
        if (board.CardCrash(c)) {
            continue;
        }
        int skip_idx = 0;
        while (true) {
            if (skip_idx == 0) {
                io_nets_per_card[ComboIdx(skip_idx, c)] = card_rank_sum[skip_idx][c] - card_sum[c];
                skip_idx++;
                continue;
            }
            int combo_idx = ComboIdx(skip_idx, c);
            // cut off if -1
            if (io_card_skipping_ranks[combo_idx] == -1) {
                break;
            }
            io_nets_per_card[combo_idx] =
                    io_nets_per_card[combo_idx - HOLDEM_MAX_DECK] + card_rank_sum[skip_idx][c]
                    +
                    card_rank_sum[skip_idx - 1][c]; //card_net[combo_idx - HOLDEM_MAX_DECK]  equals combo(skip_idx-1, c)
            skip_idx++;
        }
    }
}

inline int TermEvalKernel::ComboIdx(int rank, int card)
{
    return rank * HOLDEM_MAX_DECK + card;
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
