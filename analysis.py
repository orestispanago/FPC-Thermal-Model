import matplotlib.pyplot as plt
import numpy as np
from config import dtau, mdot, n_tubes, area
from model import run, t_amb, t_in, G_r


def namestr(obj, namespace):
    return [name for name in namespace if namespace[name] is obj][0]

def plot_temp_all_nodes(temp_array):
    rows, cols = temp_array.shape
    for n in range(len(temp_array)):
        plt.title(namestr(temp_array, globals()))
        plt.plot(np.arange(cols),temp_array[n], label=n)
    plt.ylabel("Temperature (K)")
    plt.xlabel("Time (s)")
    plt.legend(title="Node #", loc='center left', bbox_to_anchor=(1, 0.5))
    plt.show()

def plot_temp_node(temp_array, node, style=""):
    rows, cols = temp_array.shape
    x = np.arange(cols)*dtau
    y = temp_array[node]
    label = namestr(temp_array, globals())
    plt.plot(x,y,style, label=label)
    plt.ylabel("Temperature (K)")
    plt.xlabel("Time (s)")
    plt.title(f"Selected node: {node}")
    plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))

def plot_1d_array(arr, title=None, ylabel=None):
    x = np.arange(len(arr))*dtau
    y = arr
    plt.plot(x,y)
    plt.ylabel(ylabel)
    plt.xlabel("Time (s)")
    plt.title(title)

def plot_temps_dt_qdot_eff():
    plt.figure(figsize=(10, 5))
    plt.subplot(2, 2, 1)
    plot_temp_node(t_gc, selected_node)
    plot_temp_node(t_absc, selected_node)
    plot_temp_node(t_wc, selected_node)
    plot_temp_node(t_insc, selected_node)
    plot_temp_node(t_am, selected_node, style="--")
    plot_temp_node(t_inp, selected_node, style="--")
    plt.subplot(2, 2, 2)
    plot_1d_array(delta_t, title="$\Delta T$", ylabel="Temperature difference (K)")
    plt.subplot(2, 2, 3)
    plot_1d_array(Q_dot, title="$\dot Q$", ylabel="W")
    plt.subplot(2, 2, 4)
    plot_1d_array(eff, title="Efficiency")
    plt.tight_layout()
    plt.show()

t_gc, t_ac, t_absc, t_wc, t_insc = run()



plot_temp_all_nodes(t_gc)
plot_temp_all_nodes(t_ac)
plot_temp_all_nodes(t_absc)
plot_temp_all_nodes(t_wc)
plot_temp_all_nodes(t_insc)

selected_node = 6
# Create arrays for plots, make sure t_amb and t_in are constant
t_am = np.full_like(t_insc, t_amb)
t_inp = np.full_like(t_insc, t_in)

delta_t = t_wc[-1] - t_in
Q_dot = mdot * n_tubes * 4186 * delta_t
eff = np.divide(Q_dot, G_r * area, out=np.zeros_like(Q_dot), where=G_r!=0)


plot_temps_dt_qdot_eff()


# import pandas as pd

# def save_to_file(time_steps, dtau, t_in, G_r, t_abs_node, t_glass_node, Q_dot, eff):
#     df = pd.DataFrame({"time": time_steps*dtau,
#                        "t_in": t_in,
#                        "irrad": G_r[:-1],
#                        "t_abs": t_abs_node,
#                        "t_glass": t_glass_node,
#                        "Q_dot": Q_dot,
#                        "Eff": eff})
#     df.to_csv("output.csv", index=False)

# TODO: modify save_to_file to accept 2D arrays
# save_to_file(time_steps, dtau, t_in, G_r, t_absc,
#              t_gc, Q_dot, eff)