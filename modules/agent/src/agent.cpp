/* C / C++ includes */
#include <bulldog/engine.h>
#include "acpc_connector.hpp"
#include "slumbot_connector.hpp"
#include <cxxopts.hpp>
#include <filesystem>
#include "constants.hpp"
#include <cstdlib>
#include <cstdio>

int main(int argc, char *argv[])
{
    // argument parser
    // ./agent engine=acpc params=localhost,51000
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
             "[debug, info, warn, err]"
            )
            ("o,log_output", "console/file?", cxxopts::value<std::string>(), "default logs output to console")
            ("proxy", "proxy server", cxxopts::value<std::string>(), "proxy server")
            ("slumbot_session", "Slumbot session key", cxxopts::value<std::string>(), "")
            ("h,help", "print usage information");

    auto result = options.parse(argc, argv);

    if (result.count("help") || (result.arguments()).empty()) {
        std::cout << options.help() << std::endl;
        exit(0);
    }

    auto slumbot_session_key = result.count("slumbot_session") ? result["slumbot_session"].as<std::string>() : "";

    bool check_must_opts = result.count("engine")
                           || result.count("engine_params")
                           || result.count("connector")
                           || result.count("connector_params")
                           || result.count("game");
    if (!check_must_opts) {
        std::cout << "fill the must-have options" << std::endl;
        std::cout << options.help() << std::endl;
        exit(EXIT_FAILURE);
    }

    // Proxy
    bool proxy_provided = result.count("proxy") > 0;
    auto proxy_uri_str = proxy_provided ? result["proxy"].as<std::string>() : "";
    web::uri validated_proxy_uri{};
    bool is_proxy_validated;
    if (!proxy_uri_str.empty() && web::uri::validate(proxy_uri_str)) {
        validated_proxy_uri = web::uri(proxy_uri_str);
        is_proxy_validated = true;
    } else {
        is_proxy_validated = false;
        if (proxy_provided) {
            logger::warn("🚨%s is not a valid proxy URI. go without one.");
        }
    }

    // Log level
    std::string log_level = result["log_level"].as<std::string>();
    if (result.count("log_output")) {
        std::filesystem::path dir(BULLDOG_DIR_LOG);
        std::filesystem::path filename(result["log_output"].as<std::string>());
        logger::init_logger(dir / filename, log_level);
    } else {
        logger::init_logger(log_level);
    }

    // ACPC-flavor game initialization.
    Game *game = nullptr;
    std::filesystem::path dir(BULLDOG_DIR_CFG_GAME);
    std::filesystem::path filename(result["game"].as<std::string>());
    FILE *file = fopen((dir / filename).c_str(), "r");
    if (file == nullptr) {
        logger::critical(" [AGENT]: Failed to find file %s", (dir / filename));
    }
    game = readGame(file);
    if (game == nullptr) {
        logger::critical(" [AGENT]: Failed to read content of game file %s", (dir / filename));
    }

    // Engine initialization.
    auto &e_params = result["engine_params"].as<std::vector<std::string>>();
    const char *solver_filename = e_params[0].c_str();
    auto *engine = new Engine(solver_filename, game);
    TableContext table_context;
    table_context.session_game = *game;
    table_context.table_name_ = "slumbot_table";
    engine->SetTableContext(table_context);

    int sum_session_total = 0;
    int total_hands_played = 0;

    // Make stack aware. Need to do it after the engine is inited.
    // Or it will throw blueprint game not compatible error.
    // FIXME(kwok): What the hell is this?
    //  game->use_state_stack = 1;

    // Run Sessions.
    for (int i = 1; i <= result["sessions"].as<int>(); i++) {
        logger::info(" [AGENT]: RUNNING SESSION %d\n", i);

        // Connector initialization.
        BaseConnector *connector = nullptr;
        switch (result["connector"].as<int>()) {
            case acpc: {
                connector = new AcpcConnector(result["connector_params"].as<std::vector<std::string>>());
                if (connector->connect() != EXIT_SUCCESS) {
                    logger::critical(" [AGENT]: Failed to connect to ACPC dealer");
                }
                break;
            }
            case slumbot: {
                // TODO(kwok): Make this logic more elegant.
                if (is_proxy_validated) {
                    auto http_config = web::http::client::http_client_config();
                    http_config.set_proxy(web::web_proxy(validated_proxy_uri));
                    connector = new SlumbotConnector(result["connector_params"].as<std::vector<std::string>>(),
                                                     http_config
                    );
                } else {
                    connector = new SlumbotConnector(result["connector_params"].as<std::vector<std::string>>());
                }
                if (slumbot_session_key.empty()) {
                    if (!connector->connect()) {
                        logger::critical(" [AGENT]: failed to login on Slumbot");
                    }
                } else {
                    ((SlumbotConnector *) connector)->connectWithSession(slumbot_session_key);
                }
                break;
            }
            default: {
                logger::critical(" [AGENT]: Connector Not Implemented");
            }
        }

        int session_total = 0;

        /*should get the raw feed from remote server, parse it into typed objects game_feed,
         * then parse the typed obects into the poker_engine
         * */
        while (connector->get()) {
            /* Read the incoming match state */
            MatchState match_state;
            // TODO: Handle `-1` properly.
            if (connector->parse(game, &match_state) == EXIT_FAILURE) {
                logger::error(" [AGENT]: failed to parse into game state");
                break;
            }

            /* Ignoring game-over message */ // FIXME(kwok): Check if the `finished` state has been set correctly.
            if (stateFinished(&match_state.state)) {
                // TODO: 2p hack.
                std::string players[2];
                // FIXME(kwok): The number of players is not supposed to be fixed to 2.
                players[match_state.viewingPlayer] = "Lemma"; // TODO(kwok): Add engine number
                players[match_state.viewingPlayer ? 0 : 1] = "Slumbot";

                char line[MAX_LINE_LEN];
                printState(game, &match_state.state, MAX_LINE_LEN, line);
                std::string outcome = "";
                std::string player_names = "";

                for (int p = 0; p < game->numPlayers; p++) {
                    outcome += std::to_string((int) valueOfState(game, &match_state.state, p));
                    player_names += players[p];
                    if (p != (game->numPlayers - 1)) {
                        outcome += "|";
                        player_names += "|";
                    }
                }

                auto full_match_str = std::string(line) + ":" + outcome + ":" + player_names;
                logger::info(" [AGENT]: %s", full_match_str);

                total_hands_played++;

                // NOTE(kwok): It's me.
                double net = valueOfState(game, &match_state.state, match_state.viewingPlayer);
                session_total += net;
                logger::info(" [AGENT]: net win of this hand: %.3f | %s", net, full_match_str);
                logger::info(" [AGENT]: [hand %d] session total: %d", match_state.state.handId, session_total);

                engine->EvalShowdown(match_state);
                engine->RefreshEngineState();
                continue;
            }

            /* Ignore states that we are not acting in */
            // FIXME(kwok): Is this guardian code necessary?
            if (currentPlayer(game, &match_state.state) != match_state.viewingPlayer) {
                logger::debug(" [AGENT]: 🚨ignore state, not acting player");
                continue;
            }

            /* Pick an action to play. */
            char line[MAX_LINE_LEN];
            printMatchState(game, &match_state, MAX_LINE_LEN, line);
            logger::info(" [AGENT]: %s", line);
            Action action;
            //slumbot use normalized session.
            engine->GetActionBySession(match_state, action);
            //      logger::info(" [AGENT]: action returned from engine: %c%d", actionChars[action.type], action.size);

#if 0
            //detect if the action should be fixed
            //no longger needed, done within the engine interface
      //      if (!isValidAction(game, &match_state.state, 0, &action)) {
      //        logger::info(" [AGENT]: invalid action from engine: %c%d", actionChars[action.type], action.size);
      //        //try to fix it.
      //        if (!isValidAction(game, &match_state.state, 1, &action)) {
      //          logger::warn(" [AGENT]: unable to fix action");
      //          //normally it is a r20000 invalid. return call as a hack
      //          action.type = a_call;
      //        }
      //        logger::info(" [AGENT]: invalid action has been fixed to: %c%d", actionChars[action.type], action.size);
      //      }
#endif

            /* Craft a request. */
            // FIXME(kwok): If `engine->GetActionBySession` failed, `action` remains as an undefined value.
            if (connector->build(game, &action, &match_state.state) == EXIT_FAILURE) {
                logger::error(" [AGENT]: failed to build action");
                break;
            }

            /* Send it! */
            if (connector->send() == EXIT_FAILURE) {
                logger::error(" [AGENT]: failed to send action");
                break;
            }

            logger::info("\n\n\n");
        }

        logger::info(" [AGENT]: Session total = %d ", session_total);
        sum_session_total += session_total;
        logger::info(" [AGENT]: Running total = %d, Total hands played = %d", sum_session_total, total_hands_played);
        delete connector;
    }

    delete engine;
    free(game);
}
