#ifndef AUTODIDACT_MODULES_CONNECTOR_INCLUDE_SLUMBOT_CONNECTOR_HPP_
#define AUTODIDACT_MODULES_CONNECTOR_INCLUDE_SLUMBOT_CONNECTOR_HPP_

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
struct sSlumbotMatchState
{
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

class SlumbotConnector : public BaseConnector
{
public:
    SlumbotConnector(const std::vector<std::string> &params);

    SlumbotConnector(const std::vector<std::string> &params, web::http::client::http_client_config http_config);

    ~SlumbotConnector() override;

    int connect() override;

    int connectWithSession(const std::string &session_key);

    int send() override;

    int parse(const Game *game, MatchState *match_state) override;

    int build(const Game *game, Action *action, State *state) override;

    bool get() override;

    std::string action_str_;

    // FIXME(kwok): The two members should've been private.
    sSlumbotMatchState *slumbot_match_state_{};
    web::json::value previous_act_result_json_;
private:
    const char *username_{};
    const char *password_{};
    //  const int port_ = 80;
    const char *base_url_ = "http://www.slumbot.com";
    const char *cgi_uri_ = "cgi-bin/cgi_middleman";
    long int sid_{};
    unsigned int iter_{};
    unsigned int iter_rec_{};

    bool outcome_flag_ = false;
    std::string token_;

    web::http::client::http_client_config _http_client_config{};

    bool has_showed_down()
    {
        return this->previous_act_result_json_.is_null() || this->previous_act_result_json_.has_field(U("winnings"));
    }

    void reset_hand()
    {
        this->previous_act_result_json_ = web::json::value();
    }
};

void from_json(const web::json::value &, sSlumbotMatchState *);

#endif //AUTODIDACT_MODULES_CONNECTOR_INCLUDE_SLUMBOT_CONNECTOR_HPP_
