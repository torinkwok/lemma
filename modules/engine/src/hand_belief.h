#ifndef BULLDOG_MODULES_ENGINE_SRC_HAND_BELIEF_H_
#define BULLDOG_MODULES_ENGINE_SRC_HAND_BELIEF_H_

#include <map>
#include <vector>
#include <set>
#include <bulldog/logger.hpp>
#include <bulldog/card_util.h>
#include <limits>
#include <cmath>
#include "engine_util.h"

/// NOTE(kwok):
/// - A private hand belief value is a hand expectation of the opponent's private hand, throughout the
///   whole game tree from the root to terminals.
///
/// - Whenever the opponents take an action, their the belief distributions are updated via Bayes’s rule.

const double BELIEF_VECTOR_1081_DEFAULT = 1.0 / 1081.0; // 1081 = [(52 - 5)(52 - 5 - 1)] / 2
const double BELIEF_VECTOR_1326_DEFAULT = 1.0 / 1326.0; // 1326 = (52 * 51) / 2
const int FULL_HAND_BELIEF_SIZE = HOLDEM_MAX_HANDS_PERMUTATION;

static const double kBeliefPrunedFlag = -1;

struct sPrivateHandBelief
{
    sPrivateHandBelief()
    {
        SetAll(BELIEF_VECTOR_1326_DEFAULT);
    }

    explicit sPrivateHandBelief(double v)
    {
        SetAll(v);
    }

    explicit sPrivateHandBelief(sPrivateHandBelief *that)
    {
        CopyValue(that);
    }

    double belief_[FULL_HAND_BELIEF_SIZE]{};

    VectorIndex SampleHand(unsigned long &x, unsigned long &y, unsigned long &z)
    {
        return RndXorShift<double>(belief_, FULL_HAND_BELIEF_SIZE, x, y, z, (1 << 16));
    }

    void CleanCrashHands(VectorIndex v)
    {
        for (int i = 0; i < FULL_HAND_BELIEF_SIZE; i++) {
            if (VectorIdxClash(v, i)) {
                Zero(i);
            }
        }
    }

    void Zero(int idx)
    {
        if (IsPruned(idx)) {
            return;
        }
        belief_[idx] = 0;
    }

    bool IsZero(int idx)
    {
        // NOTE(kwok): In vectorized worker, rollout might make the belief become extremelly small.
        // Using epsilon with 10^-10 might not be safe enough. Comparing to 0 should be fine.
        // return fabs(belief_[idx] - 0.0) < DOUBLE_EPSILON;
        return belief_[idx] == 0;
    }

    void Prune(int idx)
    {
        belief_[idx] = kBeliefPrunedFlag;
    }

    bool IsPruned(int idx)
    {
        return fabs(belief_[idx] - kBeliefPrunedFlag) < DOUBLE_EPSILON;
    }

    /// Typically used in ranges propogations. -1 also consdiered as zero in this case.
    bool AllZero()
    {
        for (double &i: belief_) {
            if (i > 0.0) {
                return false;
            }
        }
        return true;
    }

    /// Typically used in range propogation, assuming all > 0 except -1.
    bool AllPruned()
    {
        for (double &i: belief_) {
            if (i > kBeliefPrunedFlag) {
                return false;
            }
        }
        return true;
    }

    // TODO: Test it.
    /// Typically used in range propogation, assuming all > 0 except -1.
    void Normalize()
    {
        double sum = BeliefSum();
        if (sum == 0) {
            logger::warn("hand belief sum is zero, hack to 1.0/1326");
            for (auto i = 0; i < FULL_HAND_BELIEF_SIZE; i++) {
                if (!IsPruned(i)) {
                    belief_[i] = 1.0 / FULL_HAND_BELIEF_SIZE;
                }
            }
        } else {
            double scaler = 1.0 / sum;
            for (auto i = 0; i < FULL_HAND_BELIEF_SIZE; i++) {
                if (!IsPruned(i)) {
                    belief_[i] *= scaler;
                }
            }
        }
    }

    void PrintNonZeroBelief();

    void ExcludeBoard(Board_t &board)
    {
        for (auto i = 0; i < FULL_HAND_BELIEF_SIZE; i++) {
            // if crashed with the board, prune the belief and continue
            auto high_low = FromVectorIndex(i);
            if (board.CardCrash(high_low.first) || board.CardCrash(high_low.second)) {
                Prune(i);
                continue;
            }
        }
    }

    void NormalizeExcludeBoard(Board_t &board)
    {
        ExcludeBoard(board);
        Normalize();
    }

    double BeliefSum()
    {
        double sum = 0;
        for (int i = 0; i < FULL_HAND_BELIEF_SIZE; i++) {
            if (!IsPruned(i)) { // skip the explicitly pruned only
                sum += belief_[i];
            }
        }
        return sum;
    }

    bool HandEquals(sPrivateHandBelief *that)
    {
        for (auto i = 0; i < FULL_HAND_BELIEF_SIZE; i++) {
            if (!IsPruned(i)) {
                auto a = belief_[i];
                auto b = that->belief_[i];
                if (std::fabs(a - b) >= 0.0001) {
                    logger::debug("🚨🚨🚨 this = %.15f, while that = %.15f || idx = %d",
                                  a,
                                  b,
                                  i
                    );
                    return false;
                }
            }
        }
        return true;
    }

    /// Set the probabilities for all hands, overriding the default -1
    void SetAll(double nv)
    {
        for (double &i: belief_) {
            i = nv;
        }
    }

    void DotMultiply(sPrivateHandBelief *that)
    {
        // NOTE(kwok): They must be topologically same
        for (auto i = 0; i < FULL_HAND_BELIEF_SIZE; i++) {
            if (!IsPruned(i) && !that->IsPruned(i)) {
                belief_[i] *= that->belief_[i];
            }
        }
    }

    int CountPrunedEntries()
    {
        int count = 0;
        for (auto i = 0; i < FULL_HAND_BELIEF_SIZE; i++) {
            if (IsPruned(i)) {
                count++;
            }
        }
        return count;
    }

    /// Copy everything including the pruned.
    void CopyValue(sPrivateHandBelief *that)
    {
        for (auto i = 0; i < FULL_HAND_BELIEF_SIZE; i++) {
            belief_[i] = that->belief_[i];
        }
    }

    void Scale(double factor)
    {
        for (auto i = 0; i < FULL_HAND_BELIEF_SIZE; i++) {
            if (!IsPruned(i)) { // skip 0 and -1
                belief_[i] *= factor;
            }
        }
    }

    /// NOTE(kwok): Delete all entries < 10^-6. An average range is about 1 / 1326 ~= 8 * 10^-4.
    /// Typically used in range estimation. expect > 0 except -1.
    void Purify();

    /// NOTE(kwok): Typically used in SGS range estimations.
    /// lIt should be postive only. Skipping both 0 and 1 is fine.
    int NonZeroBeliefCount();

    /// NOTE(kwok): All the pruned should be topologically alligned.
    bool TopoAligned(sPrivateHandBelief *that);

    void SetAllUnmaskedHands(double v);
};

struct Ranges
{
    explicit Ranges(int players)
    {
        num_player_ = players;
        beliefs_ = new sPrivateHandBelief[players];
    }

    explicit Ranges(Ranges *that)
    {
        num_player_ = that->num_player_;
        beliefs_ = new sPrivateHandBelief[num_player_];
        Copy(that);
    }

    // empty only for now
    Ranges(Ranges *that, const std::string &directive)
    {
        if (directive == "empty") {
            num_player_ = that->num_player_;
            beliefs_ = new sPrivateHandBelief[num_player_];
            Copy(that);
            for (int p = 0; p < num_player_; p++) {
                beliefs_[p].SetAllUnmaskedHands(0.0);
            }
        } else {
            logger::critical("unsupported range constructr mode %s", directive);
        }
    }

    void Copy(Ranges *that_range) const
    {
        for (int p = 0; p < num_player_; p++) {
            beliefs_[p].CopyValue(&that_range->beliefs_[p]);
        }
    }

    virtual ~Ranges()
    {
        delete[] beliefs_;
    }

    [[nodiscard]] double ValueSum() const
    {
        double v = 0.0;
        for (int p = 0; p < num_player_; p++) {
            v += beliefs_[p].BeliefSum();
        }
        return v;
    }

    int num_player_;
    sPrivateHandBelief *beliefs_ = nullptr;

    [[nodiscard]] bool ReturnReady(bool check_pruning) const
    {
        if (check_pruning) {
            // return if any side is all pruned. cuz it does not want to go down while it will be all
            // zero even pruning is not on
            bool ready = beliefs_[0].AllPruned() || beliefs_[1].AllPruned();
            if (ready) {
                return true;
            }
        }
        // then check the all zero condition.
        return beliefs_[0].AllZero() && beliefs_[1].AllZero();
    }
};

#endif //BULLDOG_MODULES_ENGINE_SRC_HAND_BELIEF_H_
