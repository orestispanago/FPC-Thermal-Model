import time

import numpy as np

from coefficients import coeff
from output import save_to_file
from plotting import plots

n_nodes, flowrate, total_simulation_time, t_in_0 = 8, 55, 4000, 25
# n_nodes= number of nodes along the tube.
# flowrate= the total flow rate enter the system in Liters per hour
# total_simulation_time= total running time (s).
# t_in_0= initial input temperature (C)
collector_area = 80 * 39 * 0.0254**2  # (m**2) converted from inches to m
n_tubes = 8  # number of tubes
d_in = 9.5 / 1000  # tube inner diameter
L = 1900 / 1000  # length of tubes
flow_per_tube = flowrate / n_tubes  # mass flow rate per tube (GPM)
Vdot = flow_per_tube / 3.6e+6  # volume flow rate per tube (m**3/s)
mdot = Vdot * 1000  # mass flow rate per tube (Kg/s)
w_f = 4 * Vdot / (np.pi * d_in**2)  # working fluid velocity
dz = L / (n_nodes - 1)  # spatial step (m)
dtau = dz / w_f  # maximum time step (s)
n_time_steps = round(total_simulation_time / dtau)  # number of time steps
time_steps = np.arange(n_time_steps)
selected_node = int(n_nodes / 2)
start = time.time()
# if dtau > dz / w_f:
#     print("error in flow rate")
# else:


t_amb = np.zeros(n_time_steps + 1) + (28+273.15)  # Ambient temp.
G_r = np.zeros(n_time_steps + 1)  # Heat flux of solar radiation. (W/sqm)
t_glass = np.ones(n_nodes) * 293  # initial glass temp.
t_air = np.ones(n_nodes) * 293    # initial air gap temp.
t_abs = np.ones(n_nodes) * 293    # initial absorber temp.
t_water = np.ones(n_nodes) * 293  # initial fluid temp.
t_insul = np.ones(n_nodes) * 293  # initial insulation temp.
t_glass_node = np.zeros(n_time_steps)
t_air_node = np.zeros(n_time_steps)
t_abs_node = np.zeros(n_time_steps)
t_water_node = np.zeros(n_time_steps)
t_insul_node = np.zeros(n_time_steps)
t_out = np.zeros(n_time_steps)
Q_dot = np.zeros(n_time_steps)
eff = np.zeros(n_time_steps)
t_in = np.ones(n_time_steps) * (t_in_0+273.15)


def step_down(t, trigger_time):
    if t <= trigger_time:
        return 1
    return 0


def check_convergence(n_nodes, n_converge):
    for j in range(n_nodes):
        error = np.zeros(5)
        error[0] = abs(t_glass[j] - t_glass_old[j]) / t_glass[j]
        error[1] = abs(t_air[j] - t_air_old[j]) / t_air[j]
        error[2] = abs(t_abs[j] - t_abs_old[j]) / t_abs[j]
        error[3] = abs(t_water[j] - t_water_old[j]) / t_water[j]
        error[4] = abs(t_insul[j] - t_insul_old[j]) / t_insul[j]
        for i in range(5):
            if error[i] <= 10**-4:
                n_converge += 1
            else:
                break
    return n_converge


def calc_t_glass(i, t):
    return ((t_glass_old[i] / dtau) + (B[i] * t_amb[t]) + (C[i] * t_abs[i]) + (D[i] * t_air[i]) + (E[i] * G_r[t])) / F[i]


def calc_t_air(i, t):
    return ((t_air_old[i] / dtau) + (G[i] * (t_glass[i] + t_abs[i]))) / H[i]


def calc_t_abs(i, t):
    return ((t_abs_old[i] / dtau) + (K[i] * G_r[t]) + (L[i] * t_glass[i]) + (M[i] * t_air[i]) + (O[i] * t_water[i]) + (P[i] * t_insul[i])) / Q[i]


def calc_t_insul(i, t):
    return ((t_insul_old[i] / dtau) + (V[i] * t_abs[i]) + (W[i] * t_amb[t])) / X[i]


def calc_t_water(i):
    return ((t_water_old[i] / dtau) + (R[i] * t_abs[i]) + (S[i] * t_water[i - 1] / dz)) / U[i]


for t in range(n_time_steps):
    n_converge = 0
    G_r[t] = 660*step_down(t*dtau, 3600)
    B, C, D, E, F, G, H, K, L, M, O, P, Q, R, S, U, V, W, X = coeff(
        t_glass, t_air, t_abs, t_water, t_insul, t_amb, dtau, dz, n_nodes, mdot, t, w_f)
    n_iter = 0
    while n_converge < 5 * n_nodes:
        n_iter = n_iter + 1
        t_glass_old = t_glass.copy()
        t_air_old = t_air.copy()
        t_abs_old = t_abs.copy()
        t_water_old = t_water.copy()
        t_insul_old = t_insul.copy()

        t_glass[0] = calc_t_glass(0, t)
        t_air[0] = calc_t_air(0, t)
        t_abs[0] = calc_t_abs(0, t)
        t_insul[0] = calc_t_insul(0, t)
        t_water[0] = t_in[t]
        for j in range(1, n_nodes):
            t_glass[j] = calc_t_glass(j, t)
            t_air[j] = calc_t_air(j, t)
            t_abs[j] = calc_t_abs(j, t)
            t_water[j] = calc_t_water(j)
            t_insul[j] = calc_t_insul(j, t)
        n_converge = check_convergence(n_nodes, n_converge)
    t_glass_node[t] = t_glass[selected_node]
    t_air_node[t] = t_air[selected_node]
    t_abs_node[t] = t_abs[selected_node]
    t_water_node[t] = t_water[selected_node]
    t_insul_node[t] = t_insul[selected_node]
    t_out[t] = t_water[n_nodes - 1]
    Q_dot[t] = mdot * n_tubes * 4186 * (t_out[t] - t_in[t])
    eff[t] = G_r[t] and Q_dot[t] / (G_r[t] * collector_area) or 0
print('converged')
runtime = time.time() - start
print("Runtime:", runtime)
plots(time_steps, dtau, t_glass_node, t_air_node, t_abs_node, t_water_node,
      t_insul_node, t_amb, t_in, t_out, Q_dot, eff, selected_node, n_nodes)
save_to_file(time_steps, dtau, t_in, G_r, t_abs_node, t_glass_node, Q_dot, eff)
