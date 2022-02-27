import time

import numpy as np

from coefficients import coeff
from output import save_to_file
from plotting import plots
from temp_calculations import (calc_error, calc_t_abs, calc_t_air,
                               calc_t_glass, calc_t_ins, calc_t_w)
from constants import length, width

n_nodes = 8                             # number of nodes along the tube
flowrate = 55                           # total flow rate (Liters per hour)
total_simulation_time = 4000            # total running time (s)
t_in_0 = 298                            # initial input temperature (C)
t_amb_0 = 301                           # initial input temperature (C)
collector_area = length*width           # (m2)
n_tubes = 8                             # number of tubes
d_in = 9.5 / 1000                       # tube inner diameter
lentgh = 1900 / 1000                    # length of tubes
flow_per_tube = flowrate / n_tubes      # mass flow rate per tube (GPM)
Vdot = flow_per_tube / 3.6e+6           # volume flow rate per tube (m**3/s)
mdot = Vdot * 1000                      # mass flow rate per tube (Kg/s)
w_f = 4 * Vdot / (np.pi * d_in**2)      # working fluid velocity
dz = length / (n_nodes - 1)             # spatial step (m)
dtau = dz / w_f                         # maximum time step (s)
n_time_steps = round(total_simulation_time / dtau)  # number of time steps
time_steps = np.arange(n_time_steps)
selected_node = int(n_nodes / 2)
start = time.time()
# if dtau > dz / w_f:
#     print("error in flow rate")
# else:


t_amb = np.zeros(n_time_steps + 1) + t_amb_0  # Ambient temp.
G_r = np.zeros(n_time_steps + 1)  # Heat flux of solar radiation. (W/sqm)
t_out, Q_dot, eff,  = np.zeros((3, n_time_steps))
t_in = np.zeros(n_time_steps) + t_in_0

t_glass, t_air, t_abs, t_water, t_ins = np.zeros((5, n_nodes))+293

t_gc, t_ac, t_absc, t_wc, t_insc = np.zeros((5, n_time_steps))


def step_down(t, trigger_time):
    if t <= trigger_time:
        return 1
    return 0


def check_convergence(n_nodes, n_converge):
    for j in range(n_nodes):
        error = np.zeros(5)
        error[0] = calc_error(j, t_glass, t_glass_old)
        error[1] = calc_error(j, t_air, t_air_old)
        error[2] = calc_error(j, t_abs, t_abs_old)
        error[3] = calc_error(j, t_water, t_water_old)
        error[4] = calc_error(j, t_ins, t_ins_old)
        for i in range(5):
            if error[i] <= 10**-4:
                n_converge += 1
            else:
                break
    return n_converge


for t in range(n_time_steps):
    n_converge = 0
    G_r[t] = 660*step_down(t*dtau, 3600)
    B, C, D, E, F, G, H, K, L, M, O, P, Q, R, S, U, V, W, X = coeff(
        t_glass, t_air, t_abs, t_water, t_ins,
        t_amb[t], dtau, dz, mdot, w_f)
    n_iter = 0
    while n_converge < 5 * n_nodes:
        n_iter = n_iter + 1
        t_glass_old = t_glass.copy()
        t_air_old = t_air.copy()
        t_abs_old = t_abs.copy()
        t_water_old = t_water.copy()
        t_ins_old = t_ins.copy()

        t_glass[0] = calc_t_glass(0, t, dtau, t_glass_old, t_amb, t_abs, t_air, G_r, B, C, D, E, F)
        t_air[0] = calc_t_air(0, t_air_old, t_glass, t_abs, dtau, G, H)
        t_abs[0] = calc_t_abs(0, t, dtau, t_abs_old, t_glass, t_air, t_water, t_ins, G_r, K, L, M, O, P, Q)
        t_ins[0] = calc_t_ins(0, t, dtau, t_ins_old, t_abs, t_amb, V, W, X)
        t_water[0] = t_in[t]
        for j in range(1, n_nodes):
            t_glass[j] = calc_t_glass(j, t, dtau, t_glass_old, t_amb, t_abs, t_air, G_r, B, C, D, E, F)
            t_air[j] = calc_t_air(j, t_air_old, t_glass, t_abs, dtau, G, H)
            t_abs[j] = calc_t_abs(j, t, dtau, t_abs_old, t_glass, t_air, t_water, t_ins, G_r, K, L, M, O, P, Q)
            t_water[j] = calc_t_w(j, dtau, dz, t_water_old, t_water, t_abs, R, S, U)
            t_ins[j] = calc_t_ins(j, t, dtau, t_ins_old, t_abs, t_amb, V, W, X)
        n_converge = check_convergence(n_nodes, n_converge)

    t_gc[t] = t_glass[selected_node]
    t_ac[t] = t_air[selected_node]
    t_absc[t] = t_abs[selected_node]
    t_wc[t] = t_water[selected_node]
    t_insc[t] = t_ins[selected_node]

    t_out[t] = t_water[n_nodes - 1]
    Q_dot[t] = mdot * n_tubes * 4186 * (t_out[t] - t_in[t])
    eff[t] = G_r[t] and Q_dot[t] / (G_r[t] * collector_area) or 0
print('converged')
runtime = time.time() - start
print("Runtime:", runtime)
plots(time_steps, dtau, t_gc, t_ac, t_absc, t_wc,
      t_insc, t_amb, t_in, t_out, Q_dot, eff, selected_node, n_nodes)
save_to_file(time_steps, dtau, t_in, G_r, t_absc,
             t_gc, Q_dot, eff)
