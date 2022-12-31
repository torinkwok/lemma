import subprocess
import pathlib
import os
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import threading
import queue
import time
import numpy as np

buffer_queue = queue.Queue()
avg_bru_explo_data = []
p1_bru_explo_data = []
p2_bru_explo_data = []
moving_average_window_size = 50


def run_agent(*, log=False):
    import re
    # remaining_iter = 12499950: avg_cfu = 734.506, bru_explo = (-355.634 + 35.6298)/2 = 195.632
    # remaining_iter = 12499965: avg_cfu = 0, bru_explo = (72.8961 + -2.53304)/2 = 35.1815
    int_cap = r'(\d+)'
    float_cap = r'(-?\d+(?:\.\d+)?)'
    explo_line_pattern = re.compile(
        rf'^remaining_iter = {int_cap}: avg_cfu = {float_cap}, bru_explo = \({float_cap} \+ {float_cap}\)/2 = {float_cap}$'
    )
    agent_exe_path = os.path.abspath(pathlib.Path() / '..' / 'bin/build/release-dbinfo/agent')
    os.chdir(pathlib.Path() / '..')
    cmd = agent_exe_path + \
          " --engine_params=delta/0.json --game=nlh2_200bb.game --connector=1 " \
          "--connector_params=lemma,ckp3t4kkbccHZFmosBZVsGibxz6MnaQ4Heof3uu3nkXtLwn7GVoMDhNrj6qe8ZCU,1000 " \
          "--log_level=trace " \
          "--proxy=http://127.0.0.1:6152"
    agent_proc = subprocess.Popen(
        cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True, encoding='utf-8')
    while True:
        line = agent_proc.stderr.readline()
        if not line:
            break
        if log:
            print(line, end='')
        m = re.search(explo_line_pattern, line)
        if not m:
            continue
        yield tuple(map(lambda n: abs(float(n)), m.groups()[-3:]))


def run_agent_wrapper():
    for tup in run_agent(log=True):
        print(tup)
        buffer_queue.put(tup)


def run_agent_mock():
    while True:
        for n in np.random.random(100) * 100:
            buffer_queue.put(n)
        time.sleep(3)


def reset_ax():
    ax.cla()
    # plt.title("Curve plotted using the given points")
    plt.xlabel("Iterations")
    plt.ylabel("Exploitability")
    plt.grid(linestyle=':')


def pump_anim_frame(_):
    while not buffer_queue.empty():
        tup = buffer_queue.get()
        p1_bru_explo_data.append(tup[0])
        p2_bru_explo_data.append(tup[1])
        avg_bru_explo_data.append(tup[2])
    reset_ax()
    ax.plot(list(range(len(avg_bru_explo_data))), avg_bru_explo_data)
    ax.plot(list(range(len(p1_bru_explo_data))), p1_bru_explo_data)
    ax.plot(list(range(len(p2_bru_explo_data))), p2_bru_explo_data)
    if len(avg_bru_explo_data) >= moving_average_window_size:
        smoothed = np.convolve(avg_bru_explo_data,
                               np.ones(moving_average_window_size) / moving_average_window_size)
        ax.plot(list(range(len(smoothed))), smoothed)


if __name__ == '__main__':
    matplotlib.use('WebAgg')

    fig, ax = plt.subplots(figsize=(5, 2.7), layout='constrained')
    reset_ax()

    agent_thread = threading.Thread(target=run_agent_wrapper)
    agent_thread.start()

    anim = FuncAnimation(fig, pump_anim_frame, interval=1000, blit=True)
    plt.show()
    agent_thread.join()
