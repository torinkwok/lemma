#ifndef AUTODIDACT_MODULES_ENGINE_SRC_HAND_KERNEL_HPP_
#define AUTODIDACT_MODULES_ENGINE_SRC_HAND_KERNEL_HPP_

#include "bucket_reader.hpp"
#include "hand_belief.h"
#include "term_eval_kernel.h"
#include <set>

struct sPrivateHandKernel
{
    Bucket_t bucket_by_round_vector_idx[HOLDEM_MAX_ROUNDS][FULL_HAND_BELIEF_SIZE];
    TermEvalKernel hand_eval_kernel;
    Board_t external_sampled_board;
    int starting_round;
    int valid_priv_hand_vector_idxes[HOLDEM_MAX_HANDS_PERMUTATION_EXCLUDE_BOARD];

    sPrivateHandKernel(const Board_t &board, int starting_round)
            : external_sampled_board(board), starting_round(starting_round)
    {
        int cursor = 0;
        for (int i = 0; i < FULL_HAND_BELIEF_SIZE; i++) {
            auto high_low_pair = FromVectorIndex(i);
            if (external_sampled_board.CardCrash(high_low_pair.first) ||
                external_sampled_board.CardCrash(high_low_pair.second)) {
                continue;
            }
            valid_priv_hand_vector_idxes[cursor] = i;
            cursor++;
        }
    }

    Bucket_t GetBucket(int r, int i)
    {
        auto b = bucket_by_round_vector_idx[r][i];
        if (b == INVALID_BUCKET) {
            logger::critical("should not have b = -1. should be filtered by pruned checking. at %d", i);
        }
        return b;
    }

    void AbstractHandKernel(BucketReader *bucket_reader)
    {
        for (int round = starting_round; round < HOLDEM_MAX_ROUNDS; round++) {
            for (auto vector_idx = 0; vector_idx < FULL_HAND_BELIEF_SIZE; vector_idx++) {
                auto high_low = FromVectorIndex(vector_idx);
                // skipping crash, set as -1
                if (external_sampled_board.CardCrashTillRound(high_low.first, round)
                    || external_sampled_board.CardCrashTillRound(high_low.second, round)) {
                    bucket_by_round_vector_idx[round][vector_idx] = INVALID_BUCKET;
                    continue;
                }
                auto bucket = bucket_reader->GetBucket_HighLowPair_Board_Round(
                        high_low.first, high_low.second, &external_sampled_board, round
                );
                bucket_by_round_vector_idx[round][vector_idx] = bucket;
            }
        }

#if DEV > 1
        // count -1 sum
        for (int round = starting_round; round < HOLDEM_MAX_ROUNDS; round++) {
            int sum_board = HoldemSumBoardMap[round];
            int actual_count = 0;
            int supposed_count = nCk_card(52, 2) - nCk_card(52 - sum_board, 2);
            for (auto i = 0; i < FULL_HAND_BELIEF_SIZE; i++) {
                if (bucket_by_round_vector_idx[round][i] == INVALID_BUCKET) {
                    actual_count++;
                }
            }
            if (actual_count != supposed_count) {
                logger::critical("impossible bucket count not adding up %d != %d", actual_count, supposed_count);
            }
        }
#endif

        hand_eval_kernel.Prepare(&external_sampled_board);
    };
};

#endif //AUTODIDACT_MODULES_ENGINE_SRC_HAND_KERNEL_HPP_
