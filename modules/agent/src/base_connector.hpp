#ifndef AUTODIDACT_MODULES_CONNECTOR_INCLUDE_BASE_CONNECTOR_HPP_
#define AUTODIDACT_MODULES_CONNECTOR_INCLUDE_BASE_CONNECTOR_HPP_

extern "C" {
#include "autodidact/game.h"
};

class BaseConnector
{
public:
    virtual ~BaseConnector() {};

    virtual int connect() = 0;

    virtual int parse(const Game *game, MatchState *match_state) = 0;

    //update game_feed and return response length
    virtual int build(const Game *game, Action *action, State *state) = 0;

    virtual int send() = 0;

    virtual bool get() = 0;
};

//bet size by round to bet size by game
static inline long int actionTranslate_bsbr2bsbg(Action *action, State *state, const Game *game)
{
    if (state->round > 0) {
        //                      cb300c/b300b900c/cb2400b9600
        //                      \----/\--------/\----------/
        //                         v       v         v
        //         $300 added to the pot   |         |
        //                                 v         |
        //            another $900 added to the pot  |
        //                                           v
        //                                    To be determined
        return action->size + state->sum_round_spent[state->round - 1][currentPlayer(game, state)];
    }
    return action->size;
}

//bet size by game to bet size by round
static inline long int actionTranslate_bsbg2bsbr(Action *action, State *state, const Game *game)
{
    if (state->round > 0) {
        // MATCHSTATE:0:0:20000|20000:cr300c/r600r1200c/cr3600r10800:Qc7c|/Jc7h2c/
        //                            \----/\---------/\-----------/
        //                              v        v           v
        //      Total $300 added to the pot      |           |
        //                                       v           |
        //                   Total $1200 added to the pot    |
        //                                                   v
        //                                         Total to be determined
        return action->size - state->sum_round_spent[state->round - 1][currentPlayer(game, state)];
    }
    return action->size;
}

static inline int bsbrReadAction(const char *string, const Game *game, Action *action, State *state)
{
    int c, r;

    //acpc game can automatically read 'k' and 'b' to 'c' and 'r' respectively
    if (string[0] == 'k' || string[0] == 'c') {
        action->type = a_call;
    } else if (string[0] == 'b' || string[0] == 'r') {
        action->type = a_raise;
    } else {
        action->type = a_fold;
    }

    if (action->type < 0) {
        return -1;
    }
    c = 1;

    if (action->type == a_raise && game->bettingType == noLimitBetting) {
        /* no-limit bet/raise needs to read a size */

        if (sscanf(&string[c], "%" SCNd32"%n", &action->size, &r) < 1) {
            return -1;
        }
        action->size = actionTranslate_bsbr2bsbg(action, state, game);
        c += r;
    } else {
        /* size is zero for anything but a no-limit raise */

        action->size = 0;
    }

    return c;
}

//Read Betting for Bet-Sum by Round, e.g. Slumbot / PokerTH
static inline int bsbrReadBetting(const char *string, const Game *game, State *state)
{
    int c, r;
    Action action;

    c = 0;
    while (1) {

        if (string[c] == 0) {
            break;
        }

        if (string[c] == ':') {
            ++c;
            break;
        }

        /* ignore / character */
        if (string[c] == '/') {
            ++c;
            continue;
        }

        r = bsbrReadAction(&string[c], game, &action, state);
        if (r < 0) {
            return -1;
        }

        if (!isValidAction(game, state, 0, &action)) {
            return -1;
        }

        doAction(game, &action, state);
        c += r;
    }

    return c;
}

//Read Betting for Bet-Sum by Game, e.g. ACPC
static inline int bsbgReadBetting(const char *string, const Game *game, State *state)
{
    int c, r;
    Action action;

    c = 0;
    while (1) {

        if (string[c] == 0) {
            break;
        }

        if (string[c] == ':') {
            ++c;
            break;
        }

        /* ignore / character */
        if (string[c] == '/') {
            ++c;
            continue;
        }

        r = readAction(&string[c], game, &action);
        if (r < 0) {
            return -1;
        }

        if (!isValidAction(game, state, 0, &action)) {
            return -1;
        }

        doAction(game, &action, state);
        c += r;
    }

    return c;
}

#endif //AUTODIDACT_MODULES_CONNECTOR_INCLUDE_BASE_CONNECTOR_HPP_
