#include "hand_belief.h"
bool sPrivateHandBelief::TopoAligned(sPrivateHandBelief *that) {
  for (auto i = 0; i < FULL_HAND_BELIEF_SIZE; i++) {
    if (IsPruned(i))
      if (!that->IsPruned(i)) {
        logger::error("topo not aligned at %d | this %f | that %f", i, belief_[i], that->belief_[i]);
        return false;
      }

    if (that->IsPruned(i))
      if (!IsPruned(i)) {
        logger::error("topo not aligned at %d | this %f | that %f", i, belief_[i], that->belief_[i]);
        return false;
      }
  }
  return true;
}

void sPrivateHandBelief::SetAllUnmaskedHands(double v) {
  for (auto i = 0; i < FULL_HAND_BELIEF_SIZE; i++)
    if (!IsPruned(i))
      belief_[i] = v;
}
int sPrivateHandBelief::NonZeroBeliefCount() {
  int count = 0;
  for (double i : belief_)
    if (i > 0)
      count++;
  return count;
}

//0.00005 is really really low. purify that to speed up the eval. that might also be the case why opponent is out of range often.
void sPrivateHandBelief::Purify() {
  for (auto i = 0; i < FULL_HAND_BELIEF_SIZE; i++)
    if (belief_[i] < 5 * pow(10, -5))
      Zero(i);
  Normalize();
}
void sPrivateHandBelief::PrintNonZeroBelief() {
  if (NonZeroBeliefCount() <= 100) {
    std::string output = "hand belief | ";
    for (auto i = 0; i < FULL_HAND_BELIEF_SIZE; i++) {
      if (belief_[i] > 0) {
        output += VectorIdxToString(i) + ",";
        output += std::to_string(belief_[i]) + " | ";
      }
    }
    logger::debug(output);
  }
}