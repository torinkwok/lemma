#ifndef BULLDOG_MODULES_ENGINE_SRC_TERM_EVAL_KERNEL_H_
#define BULLDOG_MODULES_ENGINE_SRC_TERM_EVAL_KERNEL_H_

#include "hand_belief.h"

static const int REGRET_SCALER = 1000; // turn into milli big blind. also counter the regret integer rounding problem

struct sPrivHandRank
{
    Card_t high_card;
    Card_t low_card;
    /// Computed and stored by external consumers, taking a external board into account.
    int rank;
    uint16_t vector_idx;

    [[nodiscard]] inline bool CheckCardCrash(Card_t check_card) const
    {
        if (low_card == check_card) {
            return true;
        } else {
            return high_card == check_card;
        }
    }

    inline std::array<Card_t, 2> GetHandPair()
    {
        return std::array<Card_t, 2>{high_card, low_card};
    }

    /// A utility function for sorting by rank.
    [[nodiscard]] inline bool RankLower(int that_rank) const
    {
        return rank < that_rank;
    }

    inline bool CardCrash(const sPrivHandRank *that) const
    {
        // NOTE(kwok): short-circuited check
        if (high_card < that->low_card) return false;
        if (low_card > that->high_card) return false;
        // NOTE(kwok): real check
        if (CheckCardCrash(that->high_card)) {
            return true;
        } else {
            return CheckCardCrash(that->low_card);
        }
    }

    inline bool RankEqual(const sPrivHandRank *that) const
    {
        return rank == that->rank;
    }

    /// Sort not only by rank but also by cards: rank > high > low
    inline bool RankHighLowSort(const sPrivHandRank *that) const
    {
        // don't even think further about it if the ranks are not equal
        if (rank < that->rank) return true;
        if (rank > that->rank) return false;
        // if the ranks are equal
        if (high_card > that->high_card) return true;
        if (high_card < that->high_card) return false;
        // if the ranks are equal, and if the high cards are equal as well
        if (low_card > that->low_card) return true;
        if (low_card < that->low_card) return false;
        return false;
    }

    void Print() const
    {
        logger::debug("high %d low %d | rank %d | idx = %d", high_card, low_card, rank, vector_idx);
    }
};

/**
 * v0 - naive O(N2) eval | 30000mcs
 * v1 - with 52C2 hands | 17000 mcs
 * v2 - with 47C2 hands | 13000 mcs
 * v3 - cache card exclusion value | 4000 mcs
 * v4 - using integer iterator not stl iteration | 3200 mcs
 * v5 - cache rank index and net value | 1200 mcs
 * v6 - passing double matrix reference | 1160 mcs
 * v7 - cache exclusion index | 300 mcs
 * v8 - A vs B. cache A drift value | 45mcs
 * v9 - building exclusion in stack up, resursively 13 mcs.
 *
 *  future:
 *  - use all primitive type better, by size.
 *  - use sparse rank (accoridng to the unique rank by my and opp.
 */


class TermEvalKernel
{
public:
    std::array<sPrivHandRank *, HOLDEM_MAX_HANDS_PERMUTATION_EXCLUDE_BOARD> showdown_sorted_hand_ranks;

    Board_t board;
    int min_rank = 0;
    int n_unique_rank = 0;
    int *rank_first_equal_index;  //rank starting
    int *rank_first_losing_index; //next rank starting
    uint16_t high_low_pos[52][52];

    //prepare related
    void Prepare(Board_t *board_ptr);

    inline void Sort();;

    void PreStack();

    //showdown
    void FastShowdownEval(double *opp_full_belief, double *my_full_belief, int spent);

    void NaiveShowdownEval(double *opp_belief, double *my_full_belief, int spent);

    void StackShowdownProb(double *opp_belief, double *rank_net_win_prob, double *card_rank_net,
                           int *card_skipping_rank_list);

    //fold eval
    void FastFoldEval(double *opp_full_belief, double *my_full_belief, int spent);

    void NaiveFoldEval(double *opp_belief, double *my_belief, int spent);

    void StackFoldingProb(double *opp_belief, double *drift_by_card, double &sum);

    void FastTerminalEval(double *opp_belief, double *my_full_belief, int spent, bool showdown);

    virtual ~TermEvalKernel()
    {
        for (auto a: showdown_sorted_hand_ranks) delete a;
        delete[] rank_first_losing_index;
        delete[] rank_first_equal_index;
    }

    static int ComboIdx(int rank, int card);
};

#endif //BULLDOG_MODULES_ENGINE_SRC_TERM_EVAL_KERNEL_H_
