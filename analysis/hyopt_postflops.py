import time
import numpy as np
import json
from hyperopt import hp, fmin, tpe

from common import run_agent


def set_hy_params(hy_params):
    print(f'⚙️{json.dumps(hy_params)}\n')
    for pot in ['0pot', 'bigpot', 'midpot']:
        with open(f'config/engine/delta/cfr/cfrs_r2_upoker-p_{pot}.json', 'w') as fh:
            json.dump(hy_params, fh, indent='    ')


def obj_func(args):
    set_hy_params(args)
    res_gen = run_agent(log=False)
    count = 0
    buffer = []
    max_iter = 1400
    try:
        for tup in res_gen:
            if count >= max_iter:
                res_gen.send('term')
            if count >= max_iter - 50:
                buffer.append(tup[-1])
            count += 1
    except StopIteration:
        avg = np.average(buffer)
        print(f'🎯avg = {avg}\n')
        return avg


def obj_func_mock(args):
    set_hy_params(args)
    time.sleep(3)
    return np.random.randint(100)


if __name__ == '__main__':
    space = {
        'cfr': {
            'num_threads': 4,
            'algo': 'vector_alt',
            'regret_matching': {
                'floor': hp.uniformint('floor', -20000.0, 0.0),
                'avg_update_on': hp.choice('avg_update_on', [-1, 0]),
                'discounting': {
                    'first_iter': 50000000,
                    'iterations': 10,
                    'interval': 10000
                }
            },
            'side_walk': hp.choice('side_walk', [True, False]),
            'rollin': {
                'estimator': {
                    'my': 'weighted_resp',
                    'opp': 'weighted_resp',
                }
            },
            'convergence': {
                'max_iter': 100000000,
                'timeout_ms': 20000000
            }
        }
    }
    optima = fmin(obj_func, space, algo=tpe.suggest, max_evals=100)
    print(f'optima: {optima}')
