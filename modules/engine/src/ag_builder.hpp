#ifndef AUTODIDACT_MODULES_ENGINE_SRC_AG_BUILDER_HPP_
#define AUTODIDACT_MODULES_ENGINE_SRC_AG_BUILDER_HPP_

#include "card_abs.hpp"
#include "action_abs.h"
#include "abstract_game.h"

class AGBuilder
{
public:
    Game *game_;

    virtual ~AGBuilder()
    {
        delete action_abs_;
        delete card_abs_;
        free(game_);
    }

    /// Prepare for the build of an abstract game.
    AGBuilder(const std::string &config_filepath, BucketPool *pool)
    {
        std::ifstream file(config_filepath);
        web::json::value config;
        if (file.good()) {
            std::stringstream buffer;
            buffer << file.rdbuf();
            config = web::json::value::parse(buffer);
            file.close();
        } else {
            logger::critical("%s is not good", config_filepath);
        }
        //storing raw config for final output
        raw_ = config;

        /* Create Game */
        std::filesystem::path game_dir(AUTODIDACT_DIR_CFG_GAME);
        FILE *game_file = fopen((game_dir / config.at("ag_builder").at("game_file").as_string()).c_str(), "rb");
        if (game_file == nullptr) {
            logger::critical("failed to open game file %s", game_file);
        }
        game_ = readGame(game_file);
        if (game_ == nullptr) {
            logger::critical("read game file failed%s", game_file);
        }
        /* initialize card abs */
        auto action_config = config.at("ag_builder").at("action_abs");

        card_abs_ = new CardAbs(config.at("ag_builder").at("card_abs"), pool);
        action_abs_ = new CompositeActionAbs(action_config);
        logger::info("ag_builder configured.");
    }

    /*
     * need to set the ag.root_reach_prob separately, unless it is an empty game.
     */
    void Build(AbstractGame *ag_out,
               State *root_state = nullptr, State *forced_state = nullptr,
               bool depth_limited = false)
    {
        if (root_state == nullptr) {
            // Building an empty game.
            State state;
            initState(game_, 0, &state);
            ag_out->root_state_ = state;
        } else {
            // Please set the root reach prob outside.
            ag_out->root_state_ = *root_state;
        }

        char line[1024];
        printState(game_, &ag_out->root_state_, 1024, line);

        logger::info("    tree root state = %s", line);
        {
            ag_out->raw_ = raw_;
            ag_out->game_ = *game_;
        }
        logger::info("ag_builder -> building betting tree...");
        {
            ag_out->depth_limited_ = depth_limited;
            ag_out->root_node_ = action_abs_->BuildBettingTree(&ag_out->game_, ag_out->root_state_, forced_state,
                                                               depth_limited);
        }
        Bucket_t bucket_count_by_round[4]{0, 0, 0, 0};
        logger::info("ag_builder -> building bucket reader...");
        {
            card_abs_->BuildReader(&ag_out->game_, &ag_out->root_state_, &ag_out->bucket_reader_);
            ag_out->bucket_reader_.GetBucketCounts(bucket_count_by_round);
        }
        logger::info("ag_builder -> building kernel...");
        {
            ag_out->BuildKernelFromRootNode(bucket_count_by_round);
            ag_out->Print();
        }
    }

private:
    web::json::value raw_;
    CompositeActionAbs *action_abs_;
    CardAbs *card_abs_;
};

#endif //AUTODIDACT_MODULES_ENGINE_SRC_AG_BUILDER_HPP_
