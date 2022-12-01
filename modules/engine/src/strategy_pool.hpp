#ifndef AUTODIDACT_MODULES_ENGINE_SRC_STRATEGY_POOL_HPP_
#define AUTODIDACT_MODULES_ENGINE_SRC_STRATEGY_POOL_HPP_

#include "strategy.h"

struct StrategyPool
{
    virtual ~StrategyPool()
    {
        for (auto s: pool_) {
            delete s;
        }
    }

    std::vector<Strategy *> pool_;
    bool sorted = false;

    Strategy *FindStrategy(MatchState *match_state, Game *game)
    {
        if (pool_.empty()) {
            logger::critical("no blueprint found. error");
        }

        if (pool_.size() == 1) {
            logger::debug("picked blueprint with depth %d",
                          GameDefaultStackDepth(&pool_.at(0)->ag_->game_));
            return pool_.at(0);
        }

        // Find by stack depth. Choose the ceiling one

        // Lazily sort the strategies by ascending stack depth if needed.
        if (!sorted) {
            std::sort(pool_.begin(), pool_.end(),
                      [](const Strategy *lhs, const Strategy *rhs)
                      {
                          auto lhs_depth = GameDefaultStackDepth(&lhs->ag_->game_);
                          auto rhs_depth = GameDefaultStackDepth(&rhs->ag_->game_);
                          if (lhs_depth == rhs_depth) {
                              logger::critical("we have blueprints of the same stack depth. not supported");
                          }
                          return lhs_depth < rhs_depth;
                      });
            sorted = true;
        }

        int match_state_stack_depth = StateStackDepth(&match_state->state, game);

        /*
         * find the blueprint with >= stack. the reason is,
         * the higher than real stack raise can be captured and translated into all in
         */
        auto it = std::find_if(pool_.begin(), pool_.end(),
                               [match_state_stack_depth](const Strategy *strategy)
                               {
                                   int depth = GameDefaultStackDepth(&strategy->ag_->game_);
                                   return depth >= match_state_stack_depth;
                               });
        Strategy *rtn_strategy;
        if (it == pool_.end()) {
            // None matched, return the hightest one.
            logger::debug("no blueprint with >= stackdepth found. use the highest one");
            rtn_strategy = pool_.back();
        } else {
            rtn_strategy = (*it);
        }

        logger::debug("picked blueprint with depth %d", GameDefaultStackDepth(&rtn_strategy->ag_->game_));
        return rtn_strategy;
    }

    void AddStrategy(Strategy *blueprint)
    {
        auto it = std::find(pool_.begin(), pool_.end(), blueprint);
        if (it == pool_.end()) {
            pool_.emplace_back(blueprint);
            sorted = false; // NOTE(kwok): Mark the order as broken. FIXME(kwok): Is this necessary?
        } else {
            logger::critical("loading duplicated blueprints");
        }
    }
};

#endif //AUTODIDACT_MODULES_ENGINE_SRC_STRATEGY_POOL_HPP_
