#ifndef AUTODIDACT_MODULES_ENGINE_SRC_CFR_WORKER_H_
#define AUTODIDACT_MODULES_ENGINE_SRC_CFR_WORKER_H_

#include "cfr_progress.h"
#include <pthread.h>
#include <vector>
#include "strategy.h"
#include "profiling.h"
#include "hand_kernel.hpp"

const double REGRET_EPSILON = 0.000000001;

const int BIASED_CALLING = 0;
const int BIASED_RAISING = 1;
const int BIASED_FOLDING = 2;
const int BIASED_NONE = 3;
const int MAX_META_STRATEGY = 4;

const int BIASED_SCALE = 5;

class CfrWorker
{
public:
    bool iter_prune_flag = false;

    CfrWorker(Strategy *blueprint,
              Strategy *strategy,
              sCfrParam *cfr_param,
              std::vector<Board_t> &my_flops,
              unsigned long long seed)
            : blueprint_(blueprint), strategy(strategy), cfr_param_(cfr_param),
              my_flops_(my_flops)
    {
        gen.seed(seed);
        std::uniform_int_distribution<int> distr(std::numeric_limits<int>::min(), std::numeric_limits<int>::max());
        x = std::abs(distr(gen));
        y = std::abs(distr(gen));
        z = std::abs(distr(gen));
    }

    // no need to delete anything. these pointers are owned by other classes.
    virtual ~CfrWorker() = default;

    Strategy *blueprint_;
    Strategy *strategy;
    sCfrParam *cfr_param_;
    std::vector<Board_t> &my_flops_;
    std::mt19937 gen;
    unsigned long x, y, z;

    CFR_MODE mode_;

    // return exploitability mbb/g
    virtual double Solve(Board_t board, bool calc_bru_explo, double *out_bru_explo) = 0;

    static double ClampRegret(double old_reg, double diff, double floor)
    {
        if (old_reg < floor) {
            logger::warn("old regret %.16f < floor %.16f", old_reg, floor);
        }
        double temp_reg = old_reg + diff; // to prevent the +diff makes reg overflow.
        if (temp_reg < floor) {
            temp_reg = floor;
        }
        return temp_reg;
    }

    void SetWalkingMode(CFR_MODE mode)
    {
        mode_ = mode;
    }
};

struct sPrivateHandsInfo
{
    unsigned long x, y, z;
    int num_players;

    /// NOTE(kwok): The board is specified by invokers. Typically sampled by a external public chance
    /// events sampler, such as `SampleSequentialFullBoard()`.
    Board_t external_sampled_board_;

    // FIXME(kwok): The number of players is not supposed to be fixed to 2.
    VectorIndex internal_sampled_priv_hands_[2];

    // FIXME(kwok): The number of players is not supposed to be fixed to 2.
    int payoff_[2];

    // FIXME(kwok): The number of players is not supposed to be fixed to 2.
    Bucket_t buckets_[2][4];

    sPrivateHandsInfo(int num_players, Board_t board, std::mt19937 &ran_gen)
            : num_players(num_players), external_sampled_board_(board)
    {
        std::uniform_int_distribution<int> distr(std::numeric_limits<int>::min(), std::numeric_limits<int>::max());
        x = std::abs(distr(ran_gen));
        y = std::abs(distr(ran_gen));
        z = std::abs(distr(ran_gen));
    }

    virtual ~sPrivateHandsInfo() = default;

    void SetBucketAndPayoff(AbstractGame *ag)
    {
        // set buckets by round
        int rank[2]; // FIXME(kwok): The number of players is not supposed to be fixed to 2.
        for (int player_pos = 0; player_pos < num_players; player_pos++) {
            auto high_low_pair = FromVectorIndex(internal_sampled_priv_hands_[player_pos]);
            rank[player_pos] = RankHand(high_low_pair.first, high_low_pair.second, &external_sampled_board_);
#if DEV > 1
            if (rank[player_pos] < 0) {
                logger::critical("rank < 0 | %d | [high %d] [low %d]", rank[player_pos], high_low_pair.first,
                                 high_low_pair.second
                );
                external_sampled_board_.Print();
            }
#endif
            for (int r = ag->root_node_->GetRound(); r <= ag->GetMaxRound(); r++) {
                auto bucket = ag->bucket_reader_.GetBucket_HighLowPair_Board_Round(high_low_pair.first,
                                                                                   high_low_pair.second,
                                                                                   &external_sampled_board_, r
                );
                // NOTE(kwok): discriminate buckets at different rounds
                buckets_[player_pos][r] = bucket;
            }
        }
        // FIXME(kwok): The number of players is not supposed to be fixed to 2.
        // set winning flag
        payoff_[0] = rank[0] > rank[1]
                     ? 1
                     : rank[0] == rank[1] ? 0 : -1;
        payoff_[1] = payoff_[0] * -1;
    }

    /// Assuming the root belief are already safe.
    void SamplePrivateHandsForAll(AbstractGame *ag, std::array<sPrivateHandBelief *, 2> &root_hand_belief)
    {
        // sample a private hand pair for player 0
        internal_sampled_priv_hands_[0] = root_hand_belief[0]->SampleHand(x, y, z);

        // sample a private hand pair for player 1
        // FIXME(kwok): The number of players is not supposed to be fixed to 2.
        while (true) {
            VectorIndex vidx_1 = root_hand_belief[1]->SampleHand(x, y, z);
            // and the new sampled pair must not crash with hand 0
            if (!VectorIdxClash(internal_sampled_priv_hands_[0], vidx_1)) {
                // we have got a legit private hand pair
                internal_sampled_priv_hands_[1] = vidx_1;
                break;
            }
        }

        // logger::debug("sampled hands | %s | %s", VectorIdxToString(vidx[0]), VectorIdxToString(vidx[1]));

#if DEV > 1
        // ensure that nothing crashes with board
        for (unsigned short i: internal_sampled_priv_hands_) {
            auto high_low = FromVectorIndex(i);
            if (external_sampled_board_.CardCrash(high_low.first) ||
                external_sampled_board_.CardCrash(high_low.second)) {
                logger::critical("error in hand sampling");
            }
        }
#endif

        SetBucketAndPayoff(ag);
    }
};

const int ACTION_DIMENSION = 0;
const int CARD_DIMENSION = 1;

class ScalarCfrWorker : public CfrWorker
{
public:
    ScalarCfrWorker(Strategy *blueprint,
                    Strategy *strategy,
                    sCfrParam *cfr_param,
                    std::vector<Board_t> &my_flops,
                    unsigned long long seed)
            :
            CfrWorker(blueprint, strategy, cfr_param, my_flops, seed)
    {
    }

    double Solve(Board_t board, bool calc_bru_explo, double *out_bru_explo) override;

    double WalkTree(int trainee, Node *this_node, sPrivateHandsInfo &hand_info,
                    std::optional<CFU_COMPUTE_MODE> trainee_cfu_compute_mode_hint,
                    std::optional<STRATEGY_TYPE> trainee_strategy_type_hint, bool learn);

    static double EvalTermNode(int trainee, Node *this_node, sPrivateHandsInfo &hand_info);

    double EvalInterNode(int trainee, Node *this_node, sPrivateHandsInfo &hand_info,
                         std::optional<CFU_COMPUTE_MODE> trainee_cfu_compute_mode_hint,
                         std::optional<STRATEGY_TYPE> trainee_strategy_type_hint, bool learn);

    // Depth-Limited Solving
    double EvalLeafRootNode(int trainee, Node *leaf_root_node, sPrivateHandsInfo &hand_info);

    void ComputeCfu(Node *this_node, const double *children_cfus, double &out_this_node_cfu,
                    CFU_COMPUTE_MODE cfu_compute_mode, const float *distr_rnb, const bool *prune_flag) const;

    void CollectChildBRUs(Node *this_node, const double *children_brus, const double &this_node_bru,
                          Strategy *target_strategy, const bool *prune_flag, sPrivateHandsInfo &hand_info);

    void CollectRegrets(Node *this_node, const double *children_cfus, const double &this_node_cfu,
                        Strategy *target_strategy, const bool *prune_flag, sPrivateHandsInfo &hand_info);

    void WavgUpdateSideWalk(int trainee_pos, Node *this_node, sPrivateHandsInfo &hand_info);

    double RolloutWalkLeafTreeWithBiasFavor(int trainee, Node *this_node, sPrivateHandsInfo &hand_info,
                                            int *bias_favors_for_all);

    double RolloutLeafRootNode(Node *leaf_root_node, sPrivateHandsInfo &hand_info);

    double RolloutLeafInterNodeWithBiasFavor(int trainee, Node *this_node, sPrivateHandsInfo &hand_info,
                                             int *bias_favors_for_all);
};

class VectorCfrWorker : public CfrWorker
{
public:
    sPrivateHandKernel *priv_hand_kernel = nullptr;

    VectorCfrWorker(Strategy *blueprint,
                    Strategy *strategy,
                    sCfrParam *cfr_param,
                    std::vector<Board_t> &my_flops, unsigned long long seed)
            :
            CfrWorker(blueprint, strategy, cfr_param, my_flops, seed)
    {
    }

    ~VectorCfrWorker() override
    {
        delete priv_hand_kernel;
    };

    double Solve(Board_t board, bool calc_bru_explo, double *out_bru_explo) override;

    Ranges *WalkTree_Pairwise(Node *this_node, Ranges *reach_ranges);

    Ranges *EvalChoiceNode_Pairwise(Node *this_node, Ranges *reach_ranges);

    /*
     * cfr helper methods
     */
    void ComputeCfu(Node *this_node,
                    std::vector<sPrivateHandBelief *> child_reach_ranges,
                    std::vector<sPrivateHandBelief *> child_cfus,
                    sPrivateHandBelief *this_node_cfu,
                    CFU_COMPUTE_MODE cfu_compute_mode,
                    const float *all_belief_distr_1dim) const;

    void
    CalcReachRange(Node *this_node, sPrivateHandBelief *belief,
                   std::vector<sPrivateHandBelief *> &child_ranges, Strategy *target_strategy,
                   STRATEGY_TYPE strategy_type);

    void ConditionalPrune();

    void
    CollectChildBRUs(Node *this_node, std::vector<sPrivateHandBelief *> child_brus,
                     sPrivateHandBelief *this_node_bru, Strategy *target_strategy) const;

    void
    CollectRegrets(Node *this_node, std::vector<sPrivateHandBelief *> child_cfus,
                   sPrivateHandBelief *this_node_cfu, Strategy *target_strategy);

    std::vector<sPrivateHandBelief *> ExtractBeliefs(std::vector<Ranges *> &ranges, int pos);

    sPrivateHandBelief *WalkTree_Alternate(Node *this_node, int trainee, sPrivateHandBelief *opp_belief,
                                           Strategy *target_strategy,
                                           std::optional<CFU_COMPUTE_MODE> trainee_cfu_compute_node_hint,
                                           std::optional<STRATEGY_TYPE> trainee_strategy_type_hint, bool learn);

    sPrivateHandBelief *EvalChoiceNode_Alternate(Node *this_node, int trainee, sPrivateHandBelief *opp_belief,
                                                 Strategy *target_strategy,
                                                 std::optional<CFU_COMPUTE_MODE> trainee_cfu_compute_mode_hint,
                                                 std::optional<STRATEGY_TYPE> trainee_strategy_type_hint, bool learn);
};

#endif //AUTODIDACT_MODULES_ENGINE_SRC_CFR_WORKER_H_
