#include "slumbot_connector.hpp"
#include <bulldog/game.h>

#include <cstdlib>
#include <numeric>
#include <cpprest/http_client.h>
#include <cpprest/filestream.h>
#include <bulldog/logger.hpp>
#include <sys/time.h>
#include <boost/algorithm/string.hpp>

//cpprest namespaces
//using namespace utility;                    // Common utilities like string conversions
//using namespace web;                        // Common features like URIs.
//using namespace web::http;                  // Common HTTP functionality
//using namespace web::http::client;          // HTTP client features
//using namespace concurrency::streams;       // Asynchronous streams

std::string acpcified_actions(std::string slumbot_actions) {
    std::string acpcified_actions = slumbot_actions;
    std::replace(acpcified_actions.begin(), acpcified_actions.end(), 'b', 'r');
    std::replace(acpcified_actions.begin(), acpcified_actions.end(), 'k', 'c');
    std::vector<std::string> streets;
    boost::split(streets, acpcified_actions, boost::is_any_of("/"));
    int max_bet = 0;
    for (size_t i = 0; i < streets.size(); i++) {
        std::string street_actions = streets[i];
        std::vector<std::string> betstrs;
        boost::split(betstrs, street_actions, boost::is_any_of("r"));
        int max_street_bet = max_bet;
        for (size_t j = 0; j < betstrs.size(); j++) {
            try {
                std::string betstr = betstrs[j];
                bool flag_call = false;
                bool flag_fold = false;
                if (betstr.length() > 1 && betstr[betstr.length() - 1] == 'c') {
                    flag_call = true;
                    std::remove(betstr.begin(), betstr.end(), 'c');
                } else if (betstr.length() > 1 and betstr[betstr.length() - 1] == 'f') {
                    flag_fold = true;
                    std::remove(betstr.begin(), betstr.end(), 'f');
                }
                int bet = std::stoi(betstr);
                bet += max_bet;
                max_street_bet = std::max(max_street_bet, bet);
                betstrs[j] = std::to_string(bet);
                if (flag_call) {
                    betstrs[j] += 'c';
                } else if (flag_fold) {
                    betstrs[j] += 'f';
                }
                betstrs[j] = 'r' + betstrs[j];
            } catch (std::invalid_argument &err) {
                continue;
            }
        }
        max_bet = max_street_bet;
        if (max_bet == 0) {
            max_bet = 100;
        }
        streets[i] = std::accumulate(betstrs.begin(),
                                     betstrs.end(),
                                     std::string{});
    }
    return boost::algorithm::join(streets, "/");
}

void from_json(const web::json::value &j, sSlumbotMatchState *s) {
    // {'old_action': 'b400c/kk/b200b400c/k', 'action': 'b400c/kk/b200b400c/kk',
    // 'client_pos': 1, 'hole_cards': ['7c', '2s'], 'board': ['Th', '9c', '4s',
    // '8c', '9d'], 'bot_hole_cards': ['Qc', 'Tc'], 'winnings': -800, 'won_pot':
    // -1600, 'session_num_hands': 2, 'baseline_winnings': -750, 'session_total':
    // -1200, 'session_baseline_total': -750}

    // p1_
    s->p1_ = j.has_field("client_pos") ? j.at("client_pos").as_integer() : 0;

    // holes_
    s->holes_ = "";
    if (j.has_field("hole_cards")) {
        auto holes = j.at("hole_cards").as_array();
        for (web::json::value h: holes) {
            s->holes_ += h.as_string();
        }
    }

    // board_
    s->board_ = "";
    if (j.has_field("board")) {
        auto slumbot_board = j.at("board").as_array();
        if (slumbot_board.size() == 0) {
            s->board_ = "";
        } else if (slumbot_board.size() == 3) {
            s->board_ = slumbot_board[0].as_string() + slumbot_board[1].as_string() + slumbot_board[2].as_string();
        } else if (slumbot_board.size() == 4) {
            s->board_ =
                    slumbot_board[0].as_string() + slumbot_board[1].as_string() + slumbot_board[2].as_string() + "/" +
                    slumbot_board[3].as_string();
        } else if (slumbot_board.size() == 5) {
            s->board_ =
                    slumbot_board[0].as_string() + slumbot_board[1].as_string() + slumbot_board[2].as_string() + "/" +
                    slumbot_board[3].as_string() + "/" + slumbot_board[4].as_string();
        } else {
            logger::error("undefined board size");
        }
    }

    // action_
    s->action_ = j.has_field("action") ? acpcified_actions(j.at("action").as_string()) : "";

    // ourb_ & oppb_
    s->ourb_ = s->oppb_ = 0xbe5;
    if (s->action_[s->action_.length() - 1] != 'c') {
        s->ourb_ *= -1; // so that ourb_ != oppb_
    }

    s->ps_ = j.has_field("ps") ? j.at("ps").as_integer() : 0; // x
    s->minb_ = j.has_field("minb") ? j.at("minb").as_integer() : 0;
    s->maxb_ = j.has_field("maxb") ? j.at("maxb").as_integer() : 0;

    //game ends state

    s->oppholes_ = "";
    if (j.has_field("bot_hole_cards")) {
        auto oppholes = j.at("bot_hole_cards").as_array();
        for (web::json::value h: oppholes) {
            s->oppholes_ += h.as_string();
        }
    }

    s->sd_ = j.has_field("sd") ? j.at("sd").as_integer() : 0; // x
    s->outcome_ = j.has_field("outcome") ? j.at("outcome").as_integer() : -1; // x
    s->stotal_ = j.has_field("stotal") ? j.at("stotal").as_integer() : 0;
    s->shands_ = j.has_field("shands") ? j.at("shands").as_integer() : 0; // x

    //account state values
    s->ltotal_ = j.has_field("ltotal") ? j.at("ltotal").as_integer() : 0;
    s->lconf_ = j.has_field("lconf") ? j.at("lconf").as_integer() : 0;
    s->lbtotal_ = j.has_field("lbtotal") ? j.at("lbtotal").as_integer() : 0;
    s->lbconf_ = j.has_field("lbconf") ? j.at("lbconf").as_integer() : 0;
    s->lhands_ = j.has_field("lhands") ? j.at("lhands").as_integer() : 0;
    s->sdtotal_ = j.has_field("sdtotal") ? j.at("sdtotal").as_integer() : 0;
    s->sdconf_ = j.has_field("sdconf") ? j.at("sdconf").as_integer() : 0;
    s->sdhands_ = j.has_field("sdhands") ? j.at("sdhands").as_integer() : 0;
    s->blbsdtotal_ = j.has_field("blbsdtotal") ? j.at("blbsdtotal").as_integer() : 0;
    s->blbsdconf_ = j.has_field("blbsdconf") ? j.at("blbsdconf").as_integer() : 0;
    s->blbsdhands_ = j.has_field("blbsdhands") ? j.at("blbsdhands").as_integer() : 0;
    s->clbsdtotal_ = j.has_field("clbsdtotal") ? j.at("clbsdtotal").as_integer() : 0;
    s->clbsdconf_ = j.has_field("lbsdconf") ? j.at("lbsdconf").as_integer() : 0;
    s->clbsdhands_ = j.has_field("clbsdhands") ? j.at("clbsdhands").as_integer() : 0;

    s->hip_ = j.has_field("hip") ? j.at("hip").as_integer() : 0; //always initialize to 0
};

SlumbotConnector::SlumbotConnector(const std::vector<std::string> &params) {
    username_ = params[0].c_str();
    password_ = params[1].c_str();
    iter_ = (unsigned int) std::strtoul(params[2].c_str(), nullptr, 0);
    iter_rec_ = (unsigned int) std::strtoul(params[2].c_str(), nullptr, 0);
    slumbot_match_state_ = new sSlumbotMatchState();
}

SlumbotConnector::SlumbotConnector(const std::vector<std::string> &params,
                                   web::http::client::http_client_config http_config) : SlumbotConnector(params) {
    _http_client_config = std::move(http_config);
}

SlumbotConnector::~SlumbotConnector() {
    delete slumbot_match_state_;
}

int SlumbotConnector::connect() {
    web::http::http_request loginRequest(web::http::methods::POST);
    loginRequest.headers().add(U("Content-Type"), U("application/json"));
    web::json::value loginRequestJsonBody;
    loginRequestJsonBody[U("username")] = web::json::value::string(U("lemma"));
    loginRequestJsonBody[U("password")] = web::json::value::string(
            U("ckp3t4kkbccHZFmosBZVsGibxz6MnaQ4Heof3uu3nkXtLwn7GVoMDhNrj6qe8ZCU"));
    loginRequest.set_body(loginRequestJsonBody);
    // FIXME(kwok): Encapsulate REST talks better.
    auto loginRequestJson = web::http::client::http_client(U("https://slumbot.com/api/login"), _http_client_config)
            .request(loginRequest)
            .then([](const web::http::http_response &response) {
                if (response.status_code() != 200) {
                    logger::error("(kwok) login returned " + std::to_string(response.status_code()));
                }
                return response.extract_json(true);
            })
            .then([this](const web::json::value &jsonObject) {
                logger::debug("(kwok) login returned:" + jsonObject.serialize());
                this->token_ = jsonObject.at(U("token")).as_string();
            });
    try { loginRequestJson.wait(); }
    catch (const std::exception &e) {
        logger::error(e.what());
        return EXIT_FAILURE;
    }
    return EXIT_FAILURE;
}

int SlumbotConnector::send() {
    web::http::http_request actRequest(web::http::methods::POST);
    actRequest.headers().add(U("Content-Type"), U("application/json"));
    web::json::value actRequestJsonBody;
    actRequestJsonBody[U("token")] = web::json::value::string(this->token_);
    actRequestJsonBody[U("incr")] = web::json::value::string(this->action_str_);
    actRequest.set_body(actRequestJsonBody);
    auto actRequestFuture = web::http::client::http_client(U("https://slumbot.com/api/act"), _http_client_config)
            .request(actRequest)
            .then([](const web::http::http_response &response) {
                if (response.status_code() != 200) {
                    logger::error("(kwok) slumbot_connector::send action returned Error" +
                                  std::to_string(response.status_code()));
                }
                return response.extract_json(true);
            })
            .then([this](const web::json::value &jsonObject) {
                logger::debug("(kwok) slumbot_connector::send action returned:" + jsonObject.serialize());
                this->previous_act_result_json_ = jsonObject;
                from_json(this->previous_act_result_json_, this->slumbot_match_state_);
            });
    try { actRequestFuture.wait(); }
    catch (const std::exception &e) {
        logger::error(e.what());
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}

bool SlumbotConnector::get() {
    // A hand has completed.
    if (this->has_showed_down()) {
        this->reset_hand();
        web::http::http_request newHandRequest(web::http::methods::POST);
        newHandRequest.headers().add(U("Content-Type"), U("application/json"));
        web::json::value newHandRequestJsonBody;
        newHandRequestJsonBody[U("token")] = web::json::value::string(this->token_);
        newHandRequest.set_body(newHandRequestJsonBody);
        // FIXME(kwok): Encapsulate REST talks better.
        auto newHandRequestFuture = web::http::client::http_client(U("https://slumbot.com/api/new_hand"),
                                                                   _http_client_config)
                .request(newHandRequest)
                .then([](const web::http::http_response &response) {
                    if (response.status_code() != 200) {
                        // FIXME(kwok): Check the type of error reported here.
                        throw std::runtime_error("(kwok) slumbot_connector::get next_hand returned Error");
                    }
                    return response.extract_json(true);
                })
                .then([this](const web::json::value &jsonObject) {
                    logger::trace("(kwok) slumbot_connector::get next_hand returned:" + jsonObject.serialize());
                    if (jsonObject.has_field(U("token"))) {
                        this->token_ = jsonObject.at(U("token")).as_string();
                    }
                    this->previous_act_result_json_ = jsonObject;
                    from_json(this->previous_act_result_json_, this->slumbot_match_state_);
                });
        try { newHandRequestFuture.wait(); }
        catch (const std::exception &e) {
            logger::error(e.what());
            return false;
        }
        return true;
        // In the mist of a hand.
    } else {
        return true;
    }
}

int SlumbotConnector::parse(const Game *game, MatchState *state) {
    if (slumbot_match_state_->p1_ >= game->numPlayers) {
        logger::critical("viewing player recieved from slumbot is not compatible with game");
    }
    initState(game, iter_rec_ - iter_, &state->state);

    // In Slumbot, `p1_ = 1` is the small blind. In ACPC, `p1_ = 0` is the small
    // blind.
    //
    // TODO(wolo): Should we move customizations such as below to a connector file?
    state->viewingPlayer = slumbot_match_state_->p1_ == 1 ? 0 : 1;
    std::string holes_str;
    if (state->viewingPlayer == 0) {
        holes_str += slumbot_match_state_->holes_;
        holes_str += "|";
        holes_str += slumbot_match_state_->oppholes_;
    } else {
        holes_str += slumbot_match_state_->oppholes_;
        holes_str += "|";
        holes_str += slumbot_match_state_->holes_;
    }
    std::string board_str;
    if (slumbot_match_state_->board_.empty()) {
        board_str = "";
    } else {
        board_str += "/";
        board_str += slumbot_match_state_->board_;
    }

    std::string ref_acpc_msg =
            "MATCHSTATE:" + std::to_string(state->viewingPlayer) + ":" + std::to_string(state->state.handId) +
            ":" + ":" + slumbot_match_state_->action_ + ":" + holes_str + board_str;
    logger::debug("slumbot_connector state in ref_acpc_msg: " + ref_acpc_msg);

    // Use readMatchStatePlus you'd like to supply custom ReadBettingFunction.
    //
    // NOTE(kwok): The `ReadBettingFunction` previously passed in is wrongly
    // `bsbrReadBetting`, which distorted the parsed match state.
    int c_acpc_msg_len = readMatchStatePlus(ref_acpc_msg.c_str(), game, state, bsbgReadBetting);
    char c_acpc_msg[MAX_LINE_LEN]; // TODO(kwok): Rename the `MAX_LINE_LEN` constant.
    printMatchState(game, state, MAX_LINE_LEN, c_acpc_msg);
    logger::debug("slumbot_connector state parsed: " + std::string(c_acpc_msg));

    if (c_acpc_msg_len < 0) {
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}

// Assemble Slumbot responses.
//
// TODO(wolo): Maybe wrap this ACPC standard to Slumbnot standard in a function?
// But all you need is the action after all.
int SlumbotConnector::build(const Game *game, Action *action, State *state) {
    if (action->type == a_call && state->round > 0 && state->numActions[state->round] == 0) {
        // (kwok) A CALL action leading a non-PRE-FLOP betting round is illegal.
        // Replace it with a CHECK. This occurs from time to time with random
        // strategy mode.
        action_str_ = "k";
        // FIXME(kwok): logger::error it when we're not in the random strategy mode.
        logger::warn("slumbot_connector built action string " + action_str_ +
                     " (in place of [c" + std::to_string(action->size) + "])");
        return EXIT_SUCCESS;
    }

    if (action->type == a_call && slumbot_match_state_->oppb_ == slumbot_match_state_->ourb_) {
        action_str_ = "k";
    } else if (action->type == a_raise) {
        action_str_ = "b";
    } else {
        action_str_ = actionChars[action->type];
    }

    if (game->bettingType == noLimitBetting && action->type == a_raise) {
        // FIXME(kwok): A hack only suitable for heads-up.
        unsigned short viewingPlayer = slumbot_match_state_->p1_ == 1 ? 0 : 1;
        unsigned short opponentPlayer = 1 - viewingPlayer;

        int32_t lower_bound_by_game = state->spent[opponentPlayer]; // NOTE(kwok): Total spent by the opponent.
        int32_t upper_bound_by_game = state->stackPlayer[viewingPlayer]; // NOTE(kwok): The depth of my stack.

        int32_t action_size_by_game = action->size;
        int32_t action_size_by_round = actionTranslate_bsbg2bsbr(action, state, game);

        // int32_t lower_bound_by_round = -1;
        // if ((lower_bound_by_round =
        //              state->spent[opponentPlayer] - state->sum_round_spent[state->round - 1][opponentPlayer]) < 0) {
        //     lower_bound_by_round = state->spent[opponentPlayer];
        // }
        //
        // bool is_all_in = action->size == state->stackPlayer[viewingPlayer];
        // bool has_opponent_bet_this_round = false;
        // if (state->round == 0) {
        //     has_opponent_bet_this_round = 0 != state->sum_round_spent[state->round][opponentPlayer];
        // } else {
        //     has_opponent_bet_this_round = state->sum_round_spent[state->round][opponentPlayer] >
        //                                   state->sum_round_spent[state->round - 1][opponentPlayer];
        // }
        // NOTE(kwok): If
        //      (1) the opponent has yet to bet for this round, or
        //      (2) this bet is an all-in,
        // `action_size_by_round` does not have to be twice as large.
        // if (!is_all_in && has_opponent_bet_this_round && action_size_by_round < lower_bound_by_round * 2) {
        //     logger::error(
        //             "💢the raise (by round) size must be at least twice as large as the opponent's last bet (by round), " +
        //             std::to_string(action_size_by_round) + "recieved");
        // }
        // NOTE(kwok): In no-limit betting strings, the raise action includes a size, which indicates the total number of
        // chips the player will have put into the pot after raising (i.e. the value they are raising to, not the value they
        // are raising by.)

        if (action_size_by_game < lower_bound_by_game || action_size_by_game > upper_bound_by_game) {
            logger::error("💢engine action size (by game): " + std::to_string(action_size_by_game) +
                          ", raise size (by game) must be inclusively within the range of [" +
                          std::to_string(lower_bound_by_game) + ", " +
                          std::to_string(upper_bound_by_game) + "]");
        }

        if (action_size_by_round < 0) {
            logger::error("💢raise (by round) must not be less than 0, "
                          + std::to_string(action_size_by_round)
                          + "recieved");
        }

        action_str_ += std::to_string(action_size_by_round);
    }

    logger::trace("slumbot_connector built action string " + action_str_);
    return EXIT_SUCCESS;
}
