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


def running_agent():
    import re
    # remaining_iter = 12499950: avg_cfu = 734.506, bru_explo = (355.634 + 35.6298)/2 = 195.632
    explo_line_pattern = re.compile(
        r'^remaining_iter = (\d+): avg_cfu = (\d+\.\d+), bru_explo = \((\d+\.\d+) \+ (\d+\.\d+)\)/2 = (\d+\.\d+)$'
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
        print(line, end='')
        m = re.search(explo_line_pattern, line)
        if not m:
            continue
        tup = tuple(map(lambda n: abs(float(n)), m.groups()[-3:]))
        buffer_queue.put(tup)


def running_agent_mock():
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

    agent_thread = threading.Thread(target=running_agent)
    agent_thread.start()

    anim = FuncAnimation(fig, pump_anim_frame, interval=1000, blit=True)
    plt.show()
    agent_thread.join()
