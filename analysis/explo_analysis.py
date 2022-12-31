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
import functools

collected_data = queue.Queue()
# total_count = 0
buffered_data = []
moving_average_window_size = 50


def running_agent():
    import re
    # line = "remaining_iter = 12499950: avg_cfu = 734.506, bru_explo = (355.634 + 35.6298)/2 = 195.632"
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
        collected_data.put(abs(float(m.groups()[-1])))


def running_agent_mock():
    while True:
        for n in np.random.random(100) * 100:
            collected_data.put(n)
        time.sleep(3)


def reset_ax():
    ax.cla()
    # ax.set_xlim(0)
    # ax.set_ylim(0, 2000)
    # plt.title("Curve plotted using the given points")
    plt.xlabel("Iterations")
    plt.ylabel("Exploitability")
    plt.grid(linestyle=':')


def animate(_):
    while not collected_data.empty():
        buffered_data.append(collected_data.get())
    reset_ax()
    ax.plot(list(range(len(buffered_data))), buffered_data)
    if len(buffered_data) >= moving_average_window_size:
        smoothed = np.convolve(buffered_data,
                               np.ones(moving_average_window_size) / moving_average_window_size)
        ax.plot(list(range(len(smoothed))), smoothed)


if __name__ == '__main__':
    matplotlib.use('webagg')

    fig, ax = plt.subplots(figsize=(5, 2.7), layout='constrained')
    reset_ax()

    agent_thread = threading.Thread(target=running_agent)
    agent_thread.start()

    anim = FuncAnimation(fig, animate, interval=1000, blit=True)
    plt.show()
    agent_thread.join()
