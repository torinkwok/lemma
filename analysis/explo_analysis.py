import matplotlib
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import threading
import queue
import time
import numpy as np
from common import run_agent

buffer_queue = queue.Queue()
avg_bru_explo_data = []
p1_bru_explo_data = []
p2_bru_explo_data = []
moving_average_window_size = 50


def run_agent_wrapper():
    for tup in run_agent(log=True):
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
