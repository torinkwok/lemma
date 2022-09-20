//
// Created by Torin Tung Kwok on 2022/9/17.
//

#include "bulldog/engine.h"
#import <filesystem>
#include "cxxopts.hpp"
#include "cpprest/json.h"
#include "../src/slumbot_connector.hpp"

int main(int argc, char *argv[]) {
    // Bucket bucket;
    // bucket.LoadClassicFromFlexbuffers();
    // auto all_cards_set = CardsetFromString("4s5s3hJhJc");
    // bucket.Get(&all_cards_set, NULL);

    cxxopts::Options options("Game Agent", "Manages connection with Poker Sites and executes Engine actions");
    options.add_options()
            ("g,game", "game file", cxxopts::value<std::string>(), "xxx.game")
            ("q,engine_params",
             "order specifc comma seperated list of parameters for engine",
             cxxopts::value<std::vector<std::string>>())
            ("c,connector", "connector for the platform we're connecting to", cxxopts::value<int>(), "0:acpc 1:slumbot")
            ("s,sessions", "number of outer sessions to run", cxxopts::value<int>()->default_value("1"))
            ("p,connector_params",
             " order specific comma seperated list of parameters for connector",
             cxxopts::value<std::vector<std::string>>())
            ("l,log_level", "log level", cxxopts::value<std::string>()->default_value("info"),
             "[debug, info, warn, err]")
            ("o,log_output", "console/file?", cxxopts::value<std::string>(), "default logs output to console")
            ("h,help", "print usage information");

    auto result = options.parse(argc, argv);

    // ACPC-flavor game initialization.
    Game *game = nullptr;
    std::filesystem::path dir(BULLDOG_DIR_CFG_GAME);
    std::filesystem::path filename(result["game"].as<std::string>());
    FILE *file = fopen((dir / filename).c_str(), "r");
    if (file == nullptr) {
        logger::critical(" [AGENT] : Failed to find file %s", (dir / filename));
    }

    game = readGame(file);
    if (game == nullptr) {
        logger::critical(" [AGENT] : Failed to read content of game file %s", (dir / filename));
    }

    auto connector = new SlumbotConnector(result["connector_params"].as<std::vector<std::string>>());

    { // All-in
        // auto json_str = R"({"action":"b200c/kb200c/kk/","board":["Ks","Kc","Th","Ts","Td"],"client_pos":0,"hole_cards":["Qh","7h"],"old_action":"b200c/kb200c/"})";
        auto json_str = R"({"action":"cb300c/b300b900c/cb2400b9600/","board":["Jc", "7h", "2c"],"client_pos":1,"hole_cards":["Qc","7c"],"old_action":"b200c/kb200c/"})";
        auto json = web::json::value::parse(json_str);
        connector->previous_act_result_json_ = json;
        from_json(connector->previous_act_result_json_, connector->slumbot_match_state_);

        MatchState match_state;
        connector->parse(game, &match_state);

        char line[MAX_LINE_LEN];
        printMatchState(game, &match_state, MAX_LINE_LEN, line);
        logger::info(" [AGENT] : %s", line);

        Action action{a_raise, 20000}; // All-in
        auto build_result = connector->build(game, &action, &match_state.state);
        assert(build_result == EXIT_SUCCESS);
        assert("b16400" == connector->action_str_);
    }

    { // Raise to 400
        auto json_str = R"({"action":"ck/kk/b100","board":["8d","4d","2s","9h"],"client_pos":1,"hole_cards":["As","8h"],"old_action":"ck/k"})";
        auto json = web::json::value::parse(json_str);
        connector->previous_act_result_json_ = json;
        from_json(connector->previous_act_result_json_, connector->slumbot_match_state_);

        MatchState match_state;
        connector->parse(game, &match_state);

        char line[MAX_LINE_LEN];
        printMatchState(game, &match_state, MAX_LINE_LEN, line);
        logger::info(" [AGENT] : %s", line);

        Action action{a_raise, 400};
        auto build_result = connector->build(game, &action, &match_state.state);
        assert(build_result == EXIT_SUCCESS);
        assert("b300" == connector->action_str_);
    }

    { // Raise to 1200
        auto json_str = R"({"action":"b300c/kk/b300","board":["Ks","Jc","6s","Jh"],"client_pos":1,"hole_cards":["9d","9c"],"old_action":"b300c/kk"})";
        auto json = web::json::value::parse(json_str);
        connector->previous_act_result_json_ = json;
        from_json(connector->previous_act_result_json_, connector->slumbot_match_state_);

        MatchState match_state;
        connector->parse(game, &match_state);

        char line[MAX_LINE_LEN];
        printMatchState(game, &match_state, MAX_LINE_LEN, line);
        logger::info(" [AGENT] : %s", line);

        Action action{a_raise, 1200};
        auto build_result = connector->build(game, &action, &match_state.state);
        assert(build_result == EXIT_SUCCESS);
        assert("b900" == connector->action_str_);
    }

    { // Raise to 2400
        auto json_str = R"({"action":"b300c/kk/b300b600b1200","board":["Ks","Jc","6s","Jh"],"client_pos":1,"hole_cards":["9d","9c"],"old_action":"b300c/kk"})";
        auto json = web::json::value::parse(json_str);
        connector->previous_act_result_json_ = json;
        from_json(connector->previous_act_result_json_, connector->slumbot_match_state_);

        MatchState match_state;
        connector->parse(game, &match_state);

        char line[MAX_LINE_LEN];
        printMatchState(game, &match_state, MAX_LINE_LEN, line);
        logger::info(" [AGENT] : %s", line);

        Action action{a_raise, 2400};
        auto build_result = connector->build(game, &action, &match_state.state);
        assert(build_result == EXIT_SUCCESS);
        assert("b1500" == connector->action_str_);
    }

    return 0;
}
