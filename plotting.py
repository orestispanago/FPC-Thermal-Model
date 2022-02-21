import matplotlib.pyplot as plt
import numpy as np


def plots(time_steps, dtau, t_glass_node, t_air_node, t_abs_node, t_water_node,
          t_insul_node, t_amb, t_in, t_out,
          Q_dot, eff, selected_node, n_nodes):
    T = time_steps * dtau
    fig = plt.figure(figsize=(10, 5))
    fig.suptitle(f"Selected node: {selected_node} of {n_nodes}")
    plt.subplot(2, 2, 1)
    plt.plot(T, t_glass_node, label="t_glass_node")
    plt.plot(T, t_air_node, label="t_air_node")
    plt.plot(T, t_abs_node, label="t_abs_node")
    plt.plot(T, t_water_node, label="t_water_node")
    plt.plot(T, t_insul_node, label="t_insul_node")
    plt.plot(T, t_amb[1:len(time_steps) + 1], label="t_am")
    plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    plt.xlabel("time (s)")
    plt.ylabel("Temperature (K)")
    plt.subplot(2, 2, 2)
    plt.plot(T, t_in, label="t_in")
    plt.plot(T, t_out, label="t_out")
    plt.plot(T, t_water_node, label="t_water_node")
    plt.xlabel("time (s)")
    plt.ylabel("Temperature (K)")
    plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    plt.subplot(2, 2, 3)
    plt.plot(T, Q_dot)
    plt.ylabel("Q_dot")
    plt.xlabel("time (s)")
    plt.subplot(2, 2, 4)
    plt.plot(T, eff)
    plt.ylabel("Efficiency")
    plt.xlabel("time (s)")
    plt.tight_layout()
    plt.show()
