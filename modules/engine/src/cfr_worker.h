#ifndef BULLDOG_MODULES_ENGINE_SRC_CFR_WORKER_H_
#define BULLDOG_MODULES_ENGINE_SRC_CFR_WORKER_H_

#include "cfr_state.h"
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

const int BIASED_SCALER = 5;

class CfrWorker
{
public:
    bool iter_prune_flag = false;

    CfrWorker(Strategy *blueprint,
              Strategy *strategy,
              sCfrParam *cfr_param,
              std::vector<Board_t> &my_flops,
              unsigned long long seed)
            : blueprint_(
            blueprint), strategy_(strategy), cfr_param_(cfr_param), my_flops_(my_flops)
    {
        gen.seed(seed);
        std::uniform_int_distribution<int> distr(std::numeric_limits<int>::min(), std::numeric_limits<int>::max());
        x = std::abs(distr(gen));
        y = std::abs(distr(gen));
        z = std::abs(distr(gen));
    }

    virtual ~CfrWorker()
    {
        //no need to delete anything. these pointers are owned by other classes.
    }

    Strategy *blueprint_;
    Strategy *strategy_;
    sCfrParam *cfr_param_;
    std::vector<Board_t> &my_flops_;
    std::mt19937 gen;
    unsigned long x, y, z;

    CFR_MODE mode_;

    //return expl mbb/g
    virtual double Solve(Board_t board) = 0;

    static double ClampRegret(double old_reg, double diff, double floor)
    {
        if (old_reg < floor)
            logger::warn("old regret %.16f < floor %.16f", old_reg, floor);

        double temp_reg = old_reg + diff; //to prevent the +diff makes reg overflow.
        if (temp_reg < floor)
            temp_reg = floor;
        return temp_reg;
    }

    void SetWalkingMode(CFR_MODE mode)
    {
        mode_ = mode;
    }
};

struct HandInfo
{
    HandInfo(int num_players, Board_t board, std::mt19937 &ran_gen) : num_players(num_players), board_(board)
    {
        std::uniform_int_distribution<int> distr(std::numeric_limits<int>::min(), std::numeric_limits<int>::max());
        x = std::abs(distr(ran_gen));
        y = std::abs(distr(ran_gen));
        z = std::abs(distr(ran_gen));
    }

    virtual ~HandInfo()
    {
    }

    void SetBucketAndPayoff(AbstractGame *ag)
    {
        // set buckets by round
        int rank[2];
        for (int player_pos = 0; player_pos < num_players; player_pos++) {
            auto high_low_pair = FromVectorIndex(hand_[player_pos]);
            rank[player_pos] = RankHand(high_low_pair.first, high_low_pair.second, &board_);
#if DEV > 1
            if (rank[player_pos] < 0) {
                logger::critical("rank < 0 | %d | [high %d] [low %d]", rank[player_pos], high_low_pair.first,
                                 high_low_pair.second);
                board_.Print();
            }
#endif
            for (int r = ag->root_node_->GetRound(); r <= ag->GetMaxRound(); r++) {
                auto bucket = ag->bucket_reader_.GetBucket_HighLowPair_Board_Round(high_low_pair.first,
                                                                                   high_low_pair.second,
                                                                                   &board_, r);
                // NOTE(kwok): Differentiate buckets at different rounds
                buckets_[player_pos][r] = bucket;
            }
        }
        // set winning flag
        payoff_[0] = rank[0] > rank[1]
                     ? 1
                     : rank[0] == rank[1] ? 0 : -1;
        payoff_[1] = payoff_[0] * -1;
    }

    // assuming the root belief are already safe
    void Sample(AbstractGame *ag, std::array<sPrivateHandBelief *, 2> &root_hand_belief)
    {
        // sample a private hand pair for player 0
        hand_[0] = root_hand_belief[0]->SampleHand(x, y, z);

        // sample a private hand pair for player 1
        while (true) {
            VectorIndex vidx_1 = root_hand_belief[1]->SampleHand(x, y, z);
            //and it should not crash with hand 0
            if (!VectorIdxClash(hand_[0], vidx_1)) {
                // we have gotten a legal hand
                hand_[1] = vidx_1;
                break;
            }
        }
        // logger::debug("sampled hands | %s | %s", VectorIdxToString(vidx[0]), VectorIdxToString(vidx[1]));
#if DEV > 1
        //none crash with board
        for (unsigned short i: hand_) {
            auto high_low = FromVectorIndex(i);
            if (board_.CardCrash(high_low.first) || board_.CardCrash(high_low.second)) {
                logger::critical("error in hand sampling");
            }
        }
#endif
        SetBucketAndPayoff(ag);
    }

    unsigned long x, y, z;
    int num_players;
    Board_t board_;

    //by sample
    VectorIndex hand_[2];
    //
    int payoff_[2];
    Bucket_t buckets_[2][4];
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
                    unsigned long long seed) :
            CfrWorker(blueprint, strategy, cfr_param, my_flops, seed)
    {
    }

    double Solve(Board_t board) override;

    double WalkTree(int trainee_pos, Node *this_node, HandInfo &hand_info);

    static double EvalTermNode(int trainee_pos, Node *this_node, HandInfo &hand_info);

    // Depth-Limited Solving
    double EvalRootLeafNode(int trainee_pos, Node *this_node, HandInfo &hand_info);

    double EvalIntermediateChoiceNode(int trainee_pos, Node *this_node, HandInfo &hand_info);

    void WavgUpdateSideWalk(int trainee_pos, Node *this_node, HandInfo &hand_info);

    double LeafRootRollout(int trainee_pos, Node *this_node, HandInfo &hand_info);

    double WalkLeafTree(int trainee_pos, Node *this_node, HandInfo &hand_info, int *c_strategy);

    double LeafChoiceRollout(int trainee_pos, Node *this_node, HandInfo &hand_info, int *p_meta);
};

class VectorCfrWorker : public CfrWorker
{
public:
    sHandKernel *hand_kernel = nullptr;

    VectorCfrWorker(Strategy *blueprint,
                    Strategy *strategy,
                    sCfrParam *cfr_param,
                    std::vector<Board_t> &my_flops, unsigned long long seed) :
            CfrWorker(blueprint, strategy, cfr_param, my_flops, seed)
    {
    }

    virtual ~VectorCfrWorker()
    {
        delete hand_kernel;
    };

    double Solve(Board_t board) override;

    Ranges *WalkTree_Pairwise(Node *this_node, Ranges *reach_ranges);

    Ranges *EvalChoiceNode_Pairwise(Node *this_node, Ranges *reach_ranges);

    /*
     * cfr helper methods
     */
    void ComputeCfu(Node *this_node,
                    std::vector<sPrivateHandBelief *> child_reach_ranges,
                    std::vector<sPrivateHandBelief *> child_cfu,
                    sPrivateHandBelief *cfu,
                    CFU_COMPUTE_MODE mode,
                    const float *p_double);

    void RangeRollout(Node *this_node, sPrivateHandBelief *belief_distr, std::vector<sPrivateHandBelief *> &child_ranges);

    void ConditionalPrune();

    void RegretLearning(Node *this_node, std::vector<sPrivateHandBelief *> child_cfu, sPrivateHandBelief *cfu);

    std::vector<sPrivateHandBelief *> ExtractBeliefs(std::vector<Ranges *> &ranges, int pos);

    sPrivateHandBelief *WalkTree_Alternate(Node *this_node, int actor, sPrivateHandBelief *opp_belief);

    sPrivateHandBelief *EvalChoiceNode_Alternate(Node *this_node, int trainee, sPrivateHandBelief *opp_belief);
};

#endif //BULLDOG_MODULES_ENGINE_SRC_CFR_WORKER_H_
