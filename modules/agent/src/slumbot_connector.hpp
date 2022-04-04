//
// Created by Carmen C on 21/2/2020.
//

#ifndef BULLDOG_MODULES_CONNECTOR_INCLUDE_SLUMBOT_CONNECTOR_HPP_
#define BULLDOG_MODULES_CONNECTOR_INCLUDE_SLUMBOT_CONNECTOR_HPP_

#include "base_connector.hpp"
#include <vector>
#include <string>
#include <cpprest/http_client.h>

/* Calculations
 * The full computation would be to divide by 100.0 (the big blind) and
    multiply by 100.0 (to get bb/100 rather than bb/h).  We optimize by
    cancelling those two terms out.
 * bb100 = ltotal / lhands
 * baseline_bb100 = lbtotal / lhands
 */
struct sSlumbotMatchState {
  unsigned short int p1_; //1 if you are player 1, else 0
  std::string holes_;
  std::string board_;
  std::string action_;
  unsigned long int ps_; //pot size
  unsigned int ourb_;
  unsigned int oppb_;
  unsigned int minb_;
  unsigned int maxb_;
  std::string oppholes_;
  int sd_; //after showdown
  int outcome_;
  int stotal_; //session total
  int shands_; //session hands
  int ltotal_;
  int lconf_;
  int lbtotal_; //lifetime baseline total
  int lbconf_; //lifetime baseline conf
  int lhands_; //number of lifetime hands
  int sdtotal_; //showdown total
  int sdconf_; //showdown conf
  int sdhands_; //number of showdown hands
  int blbsdtotal_;
  int blbsdconf_;
  int blbsdhands_;
  int clbsdtotal_;
  int clbsdconf_;
  int clbsdhands_;
  int ai_; //action sequence count
  unsigned short int hip_; //hand in progress
  std::string msg;
  std::string token;
};

class SlumbotConnector : public BaseConnector {
 public:
  SlumbotConnector(const std::vector<std::string> &params);
  ~SlumbotConnector() override;
  int connect() override;
  int send() override;
  int parse(const Game *game, MatchState *match_state) override;
  int build(const Game *game, Action *action, State *state) override;
  bool get() override;

  std::string action_str_;
 private:
  const char *username_{};
  const char *password_{};
//  const int port_ = 80;
  const char *base_url_ = "http://www.slumbot.com";
  const char *cgi_uri_ = "cgi-bin/cgi_middleman";
  long int sid_{};
  unsigned int iter_{};
  unsigned int iter_rec_{};
  sSlumbotMatchState *slumbot_match_state_{};
  bool outcome_flag_ = false;
};

class SlumbotLiteConnector : public BaseConnector {
  public:
    SlumbotLiteConnector(const std::vector<std::string> &params) {
    }

    ~SlumbotLiteConnector() override {
    }

    int connect() override {
      web::http::http_request loginRequest(web::http::methods::POST);
      loginRequest.headers().add(U("Content-Type"), U("application/json"));
      web::json::value loginRequestJsonBody;
      loginRequestJsonBody[U("username")] = web::json::value::string(U("lemma"));
      loginRequestJsonBody[U("password")] = web::json::value::string(U("ckp3t4kkbccHZFmosBZVsGibxz6MnaQ4Heof3uu3nkXtLwn7GVoMDhNrj6qe8ZCU"));
      loginRequest.set_body(loginRequestJsonBody);
      // FIXME:
      auto loginRequestJson = web::http::client::http_client(U("https://slumbot.com/api/login"))
        .request(loginRequest)
        .then([](const web::http::http_response& response) {
          if (response.status_code() != 200) {
            logger::error("(kwok) login returned " + std::to_string(response.status_code()));
          }
          return response.extract_json(true);
        })
        .then([this](const web::json::value& jsonObject) {
          logger::debug("(kwok) login returned:" + jsonObject.serialize());
          this->token_ = jsonObject.at(U("token")).as_string().c_str();
        });
      try {loginRequestJson.wait();}
      catch (const std::exception &e) {logger::error(e.what()); return EXIT_FAILURE;}
      return EXIT_FAILURE;
    }

    bool get() override {
      // A hand has completed.
      if (this->has_showed_down()) {
        web::http::http_request newHandRequest(web::http::methods::POST);
        newHandRequest.headers().add(U("Content-Type"), U("application/json"));
        web::json::value newHandRequestJsonBody;
        newHandRequestJsonBody[U("token")] = web::json::value::string(this->token_);
        newHandRequest.set_body(newHandRequestJsonBody);
        // FIXME:
        auto newHandRequestFuture = web::http::client::http_client(U("https://slumbot.com/api/new_hand"))
          .request(newHandRequest)
          .then([](const web::http::http_response& response) {
            if (response.status_code() != 200) {
              // FIXME:
              throw std::runtime_error("(kwok) slumbot_connector::get next_hand returned Error");
            }
            return response.extract_json(true);
          })
          .then([this](const web::json::value& jsonObject) {
            logger::trace("(kwok) slumbot_connector::get next_hand returned:" + jsonObject.serialize());
            this->reset_hand();
            if (jsonObject.has_field(U("token"))) {
              this->token_ = jsonObject.at(U("token")).as_string().c_str();
            }
          });
        try {newHandRequestFuture.wait();}
        catch (const std::exception &e) {logger::error(e.what()); return false;}
        return true;
      // In the mist of a hand.
      } else {
        return true;
      }
    }

    int send() override {
      web::http::http_request actRequest(web::http::methods::POST);
      actRequest.headers().add(U("Content-Type"), U("application/json"));
      web::json::value actRequestJsonBody;
      actRequestJsonBody[U("token")] = web::json::value::string(this->token_);
      actRequestJsonBody[U("incr")] = web::json::value::string(this->action_str_);
      actRequest.set_body(actRequestJsonBody);
      auto actRequestFuture = web::http::client::http_client(U("https://slumbot.com/api/act"))
        .request(actRequest)
        .then([](const web::http::http_response& response) {
          if (response.status_code() != 200) {
            logger::error("(kwok) slumbot_connector::send action returned Error" + std::to_string(response.status_code()));
          }
          return response.extract_json(true);
        })
        .then([this](const web::json::value& jsonObject) {
          logger::debug("(kwok) slumbot_connector::send action returned:" + jsonObject.serialize());
          this->previous_act_result_json_ = jsonObject;
        });
      try {actRequestFuture.wait();}
      catch (const std::exception &e) {logger::error(e.what()); return EXIT_FAILURE;}
      return EXIT_FAILURE;
    }

    int parse(const Game *game, MatchState *match_state) override {
      std::string actions(this->previous_act_result_json_.at(U("action")).as_string().c_str());
    }

    int build(const Game *game, Action *action, State *state) override {
      this->action_str_ = "";
      if (action->type == a_call /* && slumbot_match_state_->oppb_ == slumbot_match_state_->ourb_ */) {
        this->action_str_ += 'k';
      } else if (action->type == a_raise) {
        this->action_str_ += 'b';
      } else {
        this->action_str_ += actionChars[action->type];
      }
      if (game->bettingType == noLimitBetting && action->type == a_raise) {
        long int action_size = actionTranslate_bsbg2bsbr(action, state, game);
        // if (action_size < (slumbot_match_state_->minb_) || action_size > (slumbot_match_state_->maxb_)) {
        //   logger::critical("engine action size:" + std::to_string(action_size) + ", raise size must be in the range of "
        //                   + std::to_string(slumbot_match_state_->minb_) + " - "
        //                   + std::to_string(slumbot_match_state_->maxb_));
        // }
        if (action_size < 0) {
          logger::critical("raise cannot be less than 0, " + std::to_string(action_size) + "recieved");
        }
        this->action_str_ += std::to_string(action_size);
      }

      logger::trace("(kwok) slumbot_connector built action string " + this->action_str_);
      return EXIT_SUCCESS;
    }

    std::string action_str_;

  private:
    const char *username_{};
    const char *password_{};
    const char *base_url_ = "https://slumbot.com";
    const char *path_ = "/api";

    // bool showed_down_ = false;
    const char* token_{};
    web::json::value previous_act_result_json_;

    bool has_showed_down() {
      return this->previous_act_result_json_.has_field(U("winnings"));
    }

    void reset_hand() {
      this->previous_act_result_json_ = web::json::value();
    }
};

#endif //BULLDOG_MODULES_CONNECTOR_INCLUDE_SLUMBOT_CONNECTOR_HPP_
