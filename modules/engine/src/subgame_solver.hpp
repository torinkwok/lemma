#ifndef BULLDOG_MODULES_ENGINE_SRC_SUBGAME_SOLVER_HPP_
#define BULLDOG_MODULES_ENGINE_SRC_SUBGAME_SOLVER_HPP_

#include "strategy.h"
#include "ag_builder.hpp"
#include "cfr.h"
#include "cfr_state.h"

enum SUBGAME_BUILT_CODE {
    UNSUPPORTED_SUBGAME,
    SOLVE_ON_NEW_ROUND_HERO_FIRST,
    SOLVE_ON_NEW_ROUND_OPPO_FIRST,
    SKIP_RESOLVING_SUBGAME,
    RESOLVING_NONE_PRIOR_ACTION,
    RESOLVING_WITH_PRIOR_ACTION
};

inline std::map<int, std::string> SubgameBuiltCodeMap{
        {SOLVE_ON_NEW_ROUND_HERO_FIRST, "SOLVE_ON_NEW_ROUND_HERO_FIRST"},
        {SOLVE_ON_NEW_ROUND_OPPO_FIRST, "SOLVE_ON_NEW_ROUND_OPPO_FIRST"},
        {SKIP_RESOLVING_SUBGAME,        "SKIP_RESOLVING_SUBGAME"},
        {RESOLVING_NONE_PRIOR_ACTION,   "RESOLVING_NONE_PRIOR_ACTION"},
        {RESOLVING_WITH_PRIOR_ACTION,   "RESOLVING_WITH_PRIOR_ACTION"}
};

struct SubgameSolver {
    virtual ~SubgameSolver() {
        delete cfr_;
        delete ag_builder_;
        delete convergence_state_;
        delete action_chooser_;
    }

    SubgameSolver() {
        convergence_state_ = new sCFRState();
        action_chooser_ = new ActionChooser();
    }

    std::string name_ = "default_sgs";

    // Outer triggers:  default to be matching by ronud, unless distances are specified
    int active_round = -1;

    // Only use in bigpot solver for now
    double active_offtree_min = 9999999;
    int active_bet_seq_min = 999999; //default not using.
    int active_sumpot_min = 9999999; //

    // Inner triggers:
    double resolve_offtree_min = 9999999;
    int resolve_sumpot_min = 9999999; //default to active it.
    int resolve_last_root_sumpot_min = 1000000000; //default not stepping back

    // Actors
    AGBuilder *ag_builder_ = nullptr;
    CFR *cfr_ = nullptr;
    sCFRState *convergence_state_;
    ActionChooser *action_chooser_;
    STRATEGY_TYPE strategy_type = STRATEGY_REG; //default at reg

    /*
     * water fall detection
     *  - first check if pushable
     *  - then pot
     *  - then bet_seq
     *  - then off_tree
     */
    bool CheckTriggerCondition(NodeMatchResult &condition) const {
        if (condition.matched_node_->GetRound() == active_round) {
            if (condition.matched_node_->GetSumPot() >= active_sumpot_min) {
                logger::info("    [SGS %s] : triggered by sumpot %d", name_, condition.matched_node_->GetSumPot());
                return true;
            } else {
                logger::info("        [SGS %s] : node_sum_pot=%d, sgs_conf_active_sumpot_min=%d", name_, condition.matched_node_->GetSumPot(), active_sumpot_min);
            }
            if (condition.bet_similarity_dist_ >= active_bet_seq_min) {
                logger::info("    [SGS %s] : triggered by bet sequence pattern distance %d", name_, condition.bet_similarity_dist_);
                return true;
            } else {
                logger::info("        [SGS %s] : node_bet_similarity_dist_=%d, sgs_conf_active_bet_seq_min=%d", name_, condition.bet_similarity_dist_, active_bet_seq_min);
            }
            if (condition.off_tree_dist_ >= active_offtree_min) {
                logger::info("    [SGS %s] : triggered by offtree distance %f", name_, condition.off_tree_dist_);
                return true;
            } else {
                logger::info("        [SGS %s] : node_off_tree_dist_=%f, sgs_active_offtree_min=%f", name_, condition.off_tree_dist_, active_offtree_min);
            }
        } else {
            logger::info("        [SGS %s] : node_round=%d, sgs_active_round=%d", name_, condition.matched_node_->GetRound(), active_round);
        }
        logger::info("    [SGS %s] : not triggered.", name_);
        return false;
    }

    /*
     * How to build the subgame defines the subgame solver.
     *
     * We follow the design of Pluribus implementation (unsafe, but always solve from the street root)
     *  - to force the previous unseen action into the abstraction
     *  - to force the last_strategy of already happened actions of mine frozen
     *
     * Specifically,
     *  - CheckNewRound()
     *  - resolving if needed
     *  - kth_action = 1, do step back, which is equivalent to the street root
     *  - > 1, step back to the root.
     */
    int BuildSubgame(AbstractGame *ag_out,
                     Strategy *last_strategy,
                     NodeMatchResult &match_result,
                     MatchState *ref_match_state) const {
        State &ref_state = ref_match_state->state;
        auto round = ref_state.round;

        // Never resolve back to pre-flop.
        if (round == HOLDEM_ROUND_PREFLOP) {
            // TODO(kwok): If we don't panic here, an `EXC_BAD_ACCESS` exception would be thrown by `Node::GetRound()` anyway.
            // TODO(kwok): Recover from it elegantly.
            logger::critical("    [SGS %s] : according to 10.1126/science.aay2400 (https://www.science.org/doi/10.1126/science.aay2400), we don't do sub-game solve in pre-flop", name_);
            return UNSUPPORTED_SUBGAME;
        }

        /*
         * CHECK NEW ROUND where the acting player has taken none actions
         *
         * In 2-player game, `action_kth` can only be either 0 or 1.
         */

        auto action_kth = ref_state.numActions[round];

        // New round: Lemma acts first.
        if (action_kth == 0) {
            logger::debug("    [SGS %s] : built subgame [step back 0] for new round for [round = %d] [action_kth = %d]",
                          name_, round,
                          action_kth);
            ag_builder_->Build(ag_out, &ref_state, nullptr, cfr_->cfr_param_.depth_limited);
            return SOLVE_ON_NEW_ROUND_HERO_FIRST;
        }

        // New round: Some of opponents act first.
        // FIXME(kwok): The number of players is not supposed to be fixed to 2.
        // if (action_kth < ag->game_.numPlayers) {
        if (action_kth == 1) {
            if (!this->BuildResolvingSubgame_(ag_out, ref_match_state, 1)) {
                return SKIP_RESOLVING_SUBGAME;
            }
            return SOLVE_ON_NEW_ROUND_OPPO_FIRST;
        }

        // Check which trigger activated sub-game resolving: off-tree or pot-limit?
        if (match_result.off_tree_dist_ >= resolve_offtree_min) {
            logger::debug("    [SGS %s] : subgame resolving activated by off_tree_dist", name_);
        } else if ((ref_state.spent[0] + ref_state.spent[1]) >= resolve_sumpot_min) {
            logger::debug("    [SGS %s] : subgame resolving activated by min_pot", name_);
        } else {
            logger::debug("    [SGS %s] : skip resolving subgame. resolving requirement not met", name_);
            return SKIP_RESOLVING_SUBGAME;
        }

        /* Build resolving sub-game. */

        auto this_round = ref_match_state->state.round;
        auto last_root_pot = last_strategy->ag_->root_state_.spent[0] + last_strategy->ag_->root_state_.spent[0];
        // `last_strategy` belongs to the same round && sumpot_min requirements are met.
        bool needs_step_back_to_last_root = last_strategy->ag_->root_node_->GetRound() == this_round
                                            && last_root_pot >= resolve_last_root_sumpot_min;

        // Firstly we need to figure out how many steps to take back.
        uint8_t nsteps_to_reverse;
        if (needs_step_back_to_last_root) {
            // If `last_startegy` belongs to the same round, resolve to the root of `last_strategy`.
            nsteps_to_reverse =
                    action_kth - last_strategy->ag_->root_state_.numActions[this_round]; // n-steps to the last root.
            logger::debug("    [SGS %s] : needs to step back to the last root. steps to reverse = ", name_,
                          nsteps_to_reverse);
        } else {
            logger::debug("    [SGS %s] : skipping resolve to street root cuz %d < %d", name_, last_root_pot,
                          resolve_last_root_sumpot_min);
            // If `last_strategy` belongs to the prior round, resolve to the root of the current street.
            // This is not very likely to happen though.
            nsteps_to_reverse = 1;
        }

        // TODO(kwok): To uncomment this line of code.
        // logger::debug("    [SGS %s] : resolving takes [step back %d]", name_, nsteps_to_reverse);
        if (!this->BuildResolvingSubgame_(ag_out, ref_match_state, nsteps_to_reverse)) {
            return SKIP_RESOLVING_SUBGAME;
        }

        ag_out->root_node_->PrintState("    new subgame root : ");
        return (nsteps_to_reverse > 1) ? RESOLVING_WITH_PRIOR_ACTION : RESOLVING_NONE_PRIOR_ACTION;
    }

    void ConfigWithJson(const char *config_file, BucketPool *bucket_pool) {
        std::filesystem::path dir(BULLDOG_DIR_CFG_ENG);
        std::filesystem::path filename(config_file);
        std::ifstream sgs_file(dir / filename);
        web::json::value sgs_conf;

        if (sgs_file.good()) {
            std::stringstream buffer;
            buffer << sgs_file.rdbuf();
            sgs_conf = web::json::value::parse(buffer);
        }

        logger::require_critical(sgs_conf.has_field("builder_file"), "please fill in build_file || sgs conf");
        logger::require_critical(sgs_conf.has_field("cfr_file"), "please fill in cfr_file || sgs conf");
        logger::require_critical(sgs_conf.has_field("trigger"), "please fill in trigger || sgs conf");

        auto sgs_conf_str = std::string(config_file);
        name_ = sgs_conf_str.substr(0, sgs_conf_str.length() - 5);

        if (sgs_conf.has_field("action_chooser")) {
            auto action_conf = sgs_conf.at("action_chooser");
            action_chooser_->ConfFromJson(action_conf);
        }

#if 0
        //make convergence self enclosed in the cfr
        auto conv_config = sgs_conf.at("convergence");
        convergence_state_->iteration =
            conv_config.has_field("iter") ? conv_config.at("iter").as_integer() : -1;
        convergence_state_->time_milli_seconds =
            conv_config.has_field("time_ms") ? conv_config.at("time_ms").as_double() : -1;
        convergence_state_->exploitability =
            conv_config.has_field("expl") ? conv_config.at("expl").as_double() : -1;
        convergence_state_->expl_std =
            conv_config.has_field("expl_std") ? conv_config.at("expl_std").as_double() : -1;
#endif

        cfr_ = new CFR(sgs_conf.at("cfr_file").as_string().c_str());
        convergence_state_->iteration = cfr_->cfr_param_.iteration;
        cfr_->BuildCMDPipeline();
        cfr_->profiling_writer_.prefix_ = cfr_->cfr_param_.name + "_" + std::to_string(convergence_state_->iteration);

        // Setting playing strategy.
        if (sgs_conf.has_field("playing_strategy")) {
            strategy_type = StrategyMap[sgs_conf.at("playing_strategy").as_string()];
        }

        std::filesystem::path file(sgs_conf.at("builder_file").as_string());
        // `AGBuilder::Build` gets invoked within `Build*Subgame` methods.
        ag_builder_ = new AGBuilder((dir / file), bucket_pool);

        web::json::value sgs_trigger = sgs_conf.at("trigger");

        /* Required options: */ {
            logger::require_warn(sgs_trigger.has_field("active_round"),
                                 "must specify round for sgs",
                                 nullptr);
            active_round = sgs_trigger.at("active_round").as_integer();
        }

        /* Optional options: */ {
            /* Outer triggers: */ {
                if (sgs_trigger.has_field("active_offtree_min")) {
                    active_offtree_min = sgs_trigger.at("active_offtree_min").as_double();
                }
                if (sgs_trigger.has_field("active_bet_seq_min")) {
                    active_bet_seq_min = sgs_trigger.at("active_bet_seq_min").as_integer();
                }
                if (sgs_trigger.has_field("active_sumpot_min")) {
                    active_sumpot_min = sgs_trigger.at("active_sumpot_min").as_integer();
                }
            }

            /* Inner triggers: */ {
                if (sgs_trigger.has_field("resolve_offtree_min")) {
                    resolve_offtree_min = sgs_trigger.at("resolve_offtree_min").as_double();
                }
                if (sgs_trigger.has_field("resolve_sumpot_min")) {
                    resolve_sumpot_min = sgs_trigger.at("resolve_sumpot_min").as_integer();
                }
                if (sgs_trigger.has_field("resolve_last_root_sumpot_min")) {
                    resolve_last_root_sumpot_min = sgs_trigger.at("resolve_last_root_sumpot_min").as_integer();
                    // if (resolve_last_root_sumpot_min > active_sumpot_min) {
                    //     logger::error_and_exit(
                    //             "resolve_last_root_sumpot_min %d must <= active_sumpot_min %d. "
                    //             "you may solve in middle and then try to resolve it from root again",
                    //             resolve_last_root_sumpot_min, active_sumpot_min);
                    // }
                }
            }
        }
    }

private:
    /*
     * Build a sub-game with `steps_to_reverse` steps back.
     */
    bool BuildResolvingSubgame_(AbstractGame *ag_out, MatchState *ref_match_state, int steps_to_reverse) const {
        State &ref_state = ref_match_state->state;
        auto *step_back_state = new State;

        if (StepBackAction(ag_builder_->game_, &ref_state, step_back_state, steps_to_reverse) == -1) {
            logger::warn(
                    "    [SGS %s] : stepping too many steps. you may have an empty state or invalid state. return false",
                    name_);
            return false;
        }

        if (ref_state.round != step_back_state->round) {
            logger::error("    [SGS %s] : can not step back to the last round", name_);
            return false;
        }

        ag_builder_->Build(ag_out, step_back_state, &ref_state, cfr_->cfr_param_.depth_limited);
        delete step_back_state;
        logger::debug("    [SGS %s] : built subgame [step back %d] for r = %d", name_, steps_to_reverse,
                      ref_state.round);
        return true;
    }
};

#endif //BULLDOG_MODULES_ENGINE_SRC_SUBGAME_SOLVER_HPP_
