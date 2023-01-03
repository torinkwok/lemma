import time
import numpy as np
import json
from hyperopt import hp, fmin, tpe
import common
import argparse


class LemmaHyperOpt:
    def __init__(self):
        cmd_args_parser = argparse.ArgumentParser(prog='hyopt')
        cmd_args_parser.add_argument('-k', '--mock_response_key', type=str)
        cmd_args_parser.add_argument('-i', '--commit_id', type=str)
        cmd_args_parser.add_argument('-r', '--round', type=int)
        cmd_args_parser.add_argument('-m', '--use_obj_func_mock', action='store_true')
        cmd_args_parser.add_argument('-v', '--verbose', action='store_true')
        cmd_args_parser.add_argument('-nt', '--n_threads', type=int, default=4)
        cmd_args = cmd_args_parser.parse_args()
        self.r = cmd_args.round
        self.commit_id = cmd_args.commit_id
        self.mock_response_key = cmd_args.mock_response_key
        self.use_obj_func_mock = cmd_args.use_obj_func_mock
        self.verbose = cmd_args.verbose
        self.n_threads = cmd_args.n_threads
        del cmd_args_parser, cmd_args
        if self.r == 1:
            self.max_iter, self.avg_window, self.mccfr_var, self.mccfr_var_abbr, self.max_trials = \
                36000, 2500, 'scalar', 'cfrs', 2000
        elif self.r == 2:
            self.max_iter, self.avg_window, self.mccfr_var, self.mccfr_var_abbr, self.max_trials = \
                1400, 50, 'vector_alt', 'cfrv', 500
        elif self.r == 3:
            self.max_iter, self.avg_window, self.mccfr_var, self.mccfr_var_abbr, self.max_trials = \
                1400, 50, 'vector_alt', 'cfrv', 500
        self.is_depth_limited_search = self.r == 1
        self.triggers = ['0pot', 'bigpot', 'midpot']
        self.optima = None
        self.space = {
            'cfr': {
                'num_threads': self.n_threads,
                'algo': 'scalar',
                'regret_matching': {
                    'floor': hp.uniformint('floor', -20000.0, 0.0),
                    'avg_update_on': hp.choice('avg_update_on', [-1, 0]),
                    'discounting': {
                        'first_iter': 50000000,
                        'iterations': 10,
                        'interval': 10000
                    },
                },
                'side_walk': hp.choice('side_walk', [True, False]),
                'rollin': {
                    'estimator': {
                        'my': 'weighted_resp',
                        'opp': 'weighted_resp',
                    },
                },
                'convergence': {
                    'max_iter': 100000000,
                    'timeout_ms': 20000000
                },
            },
        }
        if self.is_depth_limited_search:
            self.space['cfr']['depth_limited'] = {
                'on': True,
                'reps': hp.uniformint('reps', 10, 200),
                'cache': hp.choice('cache', [True, False]),
            }
            self.space['cfr']['rollout'] = {
                'pruning': {
                    'regret_thres': hp.uniform('regret_thres', 0.0, 1.0),
                    'prob': hp.uniform('prob', 0.0, 1.0),
                    'first_iter': hp.uniformint('first_iter', 0, 20000),
                },
            }

    def set_hy_params(self, hy_params):
        print('-' * 10, end='\n\n')
        print(f'⚙️{json.dumps(hy_params)}\n')
        for trigger in self.triggers:
            with open(f'config/engine/delta/cfr/{self.mccfr_var_abbr}_r{self.r}_upoker-p_{trigger}.json', 'w') as fh:
                json.dump(hy_params, fh, indent='    ')

    def obj_func(self, args):
        self.set_hy_params(args)
        res_gen = common.run_agent(log=self.verbose, commit_id=self.commit_id, mock_response_key=self.mock_response_key)
        count = 0
        buffer = []
        try:
            for tup in res_gen:
                if count >= self.max_iter:
                    res_gen.send('term')
                if count >= self.max_iter - self.avg_window:
                    buffer.append(tup[-1])
                count += 1
        except StopIteration:
            avg = np.average(buffer)
            print(f'🎯avg = {avg}\n')
            return avg

    def obj_func_mock(self, args):
        self.set_hy_params(args)
        time.sleep(1)
        return np.random.randint(100)

    def find_optima(self):
        self.optima = fmin(
            self.obj_func_mock if self.use_obj_func_mock else self.obj_func,
            self.space,
            algo=tpe.suggest,
            max_evals=self.max_trials,
            verbose=True,
        )


if __name__ == '__main__':
    hyper_opt = LemmaHyperOpt()
    hyper_opt.find_optima()
    print(f'optima: {hyper_opt.optima}')
