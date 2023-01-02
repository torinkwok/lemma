import subprocess
import pathlib
import os
import re

uint_cap = r'(\d+)'
float_scinot_pat = r'(?:e[+-]\d+)?'
float_decimal_pat = rf'(?:\.\d+{float_scinot_pat})?'
float_cap = rf'(-?\d+{float_decimal_pat})'
explo_line_pat = re.compile(
    rf'^thread_iter_num = {uint_cap}, avg_cfu = {float_cap}, bru_explo = \({float_cap} \+ {float_cap}\)/2 = {float_cap}$'
)

agent_loaded_cfr_pat = re.compile(r'^🧠')

agent_exe_path = os.path.abspath(pathlib.Path() / '..' / 'bin/build/release-min/agent')
os.chdir(pathlib.Path() / '..')


def run_agent(*, log=False):
    cmd = agent_exe_path + \
          " --engine_params=delta/0.json --game=nlh2_200bb.game --connector=1 " \
          "--connector_params=lemma,ckp3t4kkbccHZFmosBZVsGibxz6MnaQ4Heof3uu3nkXtLwn7GVoMDhNrj6qe8ZCU,1000 " \
          "--log_level=trace " \
          "--proxy=http://127.0.0.1:6152"
    agent_proc = subprocess.Popen(
        cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True, encoding='utf-8')
    while True:
        stderr_line = agent_proc.stderr.readline()
        if not stderr_line:
            break
        if log or re.search(agent_loaded_cfr_pat, stderr_line):
            # log unconditionally or only log final CFR config loaded by agent executable
            print(stderr_line, end='')
        m = re.search(explo_line_pat, stderr_line)
        if not m:
            continue
        cmd = yield tuple(map(lambda n: abs(float(n)), m.groups()[-4:]))
        if cmd == 'term':
            agent_proc.terminate()
            break


if __name__ == '__main__':
    for test_line in [
        'thread_iter_num = 3, avg_cfu = 734.506, bru_explo = (-355.634 + 35.6298)/2 = 195.632',
        'thread_iter_num = 253, avg_cfu = 0, bru_explo = (72.8961 + -2.53304)/2 = 35.1815',
        'thread_iter_num = 4735, avg_cfu = -7129.78, bru_explo = (1.12499e-310 + 1.98149e+232)/2 = -5186',
    ]:
        test_match = re.search(explo_line_pat, test_line)
        print(int(test_match.groups()[0]), end=': ')
        print(list(map(lambda s: float(s), test_match.groups()[-3:])))
