import time
import numpy as np
import json
from hyperopt import hp, fmin, tpe

from common import run_agent


def set_hy_params(hy_params):
    print('-' * 10, end='\n\n')
    print(f'⚙️{json.dumps(hy_params)}\n')
    for pot in ['0pot', 'bigpot', 'midpot']:
        with open(f'config/engine/delta/cfr/cfrs_r1_upoker-p_{pot}.json', 'w') as fh:
            json.dump(hy_params, fh, indent='    ')


def obj_func(args):
    set_hy_params(args)
    res_gen = run_agent(log=False)
    count = 0
    buffer = []
    max_iter = 36000
    try:
        for tup in res_gen:
            if count >= max_iter:
                res_gen.send('term')
            if count >= max_iter - 2500:
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
            'num_threads': 36,
            'algo': 'scalar',
            'regret_matching': {
                'floor': hp.uniformint('floor', -20000.0, 0.0),
                'avg_update_on': hp.choice('avg_update_on', [-1, 0]),
                'discounting': {
                    'first_iter': 50000000,
                    'iterations': 10,
                    'interval': 10000
                }
            },
            'depth_limited': {
                'on': True,
                'reps': hp.uniformint('reps', 10, 200),
                'cache': hp.choice('cache', [True, False]),
            },
            'rollout': {
                'pruning': {
                    'regret_thres': hp.uniform('regret_thres', 0.0, 1.0),
                    'prob': hp.uniform('prob', 0.0, 1.0),
                    'first_iter': hp.uniformint('first_iter', 0, 20000),
                },
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
    optima = fmin(obj_func_mock, space, algo=tpe.suggest, max_evals=2000)
    print(f'optima: {optima}')
