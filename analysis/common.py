import subprocess
import pathlib
import os

agent_exe_path = os.path.abspath(pathlib.Path() / '..' / 'bin/build/release-min/agent')
os.chdir(pathlib.Path() / '..')


def run_agent(*, log=False):
    import re
    # remaining_iter = 12499950: avg_cfu = 734.506, bru_explo = (-355.634 + 35.6298)/2 = 195.632
    # remaining_iter = 12499965: avg_cfu = 0, bru_explo = (72.8961 + -2.53304)/2 = 35.1815
    int_cap = r'(\d+)'
    float_cap = r'(-?\d+(?:\.\d+)?)'
    explo_line_pattern = re.compile(
        rf'^remaining_iter = {int_cap}: avg_cfu = {float_cap}, bru_explo = \({float_cap} \+ {float_cap}\)/2 = {float_cap}$'
    )
    agent_loaded_cfr_pattern = re.compile(r'^🧠')
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
        if log or re.search(agent_loaded_cfr_pattern, stderr_line):
            # log unconditionally or only log final CFR config loaded by agent executable
            print(stderr_line, end='')
        m = re.search(explo_line_pattern, stderr_line)
        if not m:
            continue
        cmd = yield tuple(map(lambda n: abs(float(n)), m.groups()[-4:]))
        if cmd == 'term':
            agent_proc.terminate()
            break
