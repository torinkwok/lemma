#ifndef AUTODIDACT_CONSTANT_H
#define AUTODIDACT_CONSTANT_H

#include <cstdint>
#include <algorithm>
#include <libcuckoo/cuckoohash_map.hh>

using namespace libcuckoo;

extern "C" {
#include "autodidact/game.h"
}

using ULONG_WAVG = uint64_t;
using UINT_WAVG = uint32_t;
using DOUBLE_REGRET = double;
using INT_REGRET = int;
using ZIPAVG = uint8_t;

/// strategy purification, etc.
const double DEFAULT_BAYESIAN_TRANSITION_FILTER = 0.03;
const double RANGE_ROLLOUT_PRUNE_THRESHOLD = 0.01;

const double DOUBLE_EPSILON = 0.00000000001;

void CheckAvgSum(float *avg, int size);

bool IsAvgUniform(float *avg, int size);

void NormalizePolicy(float *avg, int size);

template<typename T>
int GetPolicy(float *inout_distr, int size, cuckoohash_map<size_t, T> *regrets, size_t offset = 0)
{
    // bool integral = std::is_integral<T>::value;
    T positive_v[size];
    T sum_pos_v = 0;
    for (int a = 0; a < size; a++) {
        // TODO(kwok): 🐦
        regrets->template insert(offset + a);
        regrets->template find_fn(
                offset + a,
                [&](const auto &regret)
                {
                    positive_v[a] = std::max<T>(0, regret);
                    sum_pos_v += positive_v[a];
                }
        );
    }
    // normalization
    if (sum_pos_v > 0) {
        for (int a = 0; a < size; a++) {
            inout_distr[a] = (float) positive_v[a] / sum_pos_v;
            // necessary?
            if (inout_distr[a] < 0) {
                return 1;
            }
        }
    } else {
        for (int a = 0; a < size; a++) {
            inout_distr[a] = (float) 1.0 / (float) size;
        }
    }
    return 0;
}

template<typename T>
int GetPolicy(float *inout_distr, int size, T *regrets, size_t offset = 0)
{
    // bool integral = std::is_integral<T>::value;
    T positive_v[size];
    T sum_pos_v = 0;
    for (int a = 0; a < size; a++) {
        positive_v[a] = std::max<T>(0, regrets[offset + a]);
        //    T v = regrets[a];
        //    positive_v[a] = v > 0 ? v : 0;;
        //      new_pos_reg[a] = regret_[rnba] > 0.0 ? regret_[rnba] : 0.0;
        // in multithreading setting this may have problem.
        sum_pos_v += positive_v[a];
    }
    // normalization
    if (sum_pos_v > 0) {
        for (int a = 0; a < size; a++) {
            inout_distr[a] = (float) positive_v[a] / sum_pos_v;
            // necessary?
            if (inout_distr[a] < 0) {
                return 1;
            }
        }
    } else {
        for (int a = 0; a < size; a++) {
            inout_distr[a] = (float) 1.0 / (float) size;
        }
    }
    return 0;
}

#endif //AUTODIDACT_CONSTANT_H
