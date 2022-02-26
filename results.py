import time

import numpy as np

from coefficients import coeff
from components import Absorber, AirGap, GlassCover, Insulation, Water
from output import save_to_file
from plotting import plots

n_nodes = 8                             # number of nodes along the tube
flowrate = 55                           # total flow rate (Liters per hour)
total_simulation_time = 4000            # total running time (s)
t_in_0 = 25                             # initial input temperature (C)
collector_area = 80 * 39 * 0.0254**2    # (m**2) converted from inches to m
n_tubes = 8                             # number of tubes
d_in = 9.5 / 1000                       # tube inner diameter
L = 1900 / 1000                         # length of tubes
flow_per_tube = flowrate / n_tubes      # mass flow rate per tube (GPM)
Vdot = flow_per_tube / 3.6e+6           # volume flow rate per tube (m**3/s)
mdot = Vdot * 1000                      # mass flow rate per tube (Kg/s)
w_f = 4 * Vdot / (np.pi * d_in**2)      # working fluid velocity
dz = L / (n_nodes - 1)                  # spatial step (m)
dtau = dz / w_f                         # maximum time step (s)
n_time_steps = round(total_simulation_time / dtau)  # number of time steps
time_steps = np.arange(n_time_steps)
selected_node = int(n_nodes / 2)
start = time.time()
# if dtau > dz / w_f:
#     print("error in flow rate")
# else:


t_amb = np.zeros(n_time_steps + 1) + (28+273.15)  # Ambient temp.
G_r = np.zeros(n_time_steps + 1)  # Heat flux of solar radiation. (W/sqm)
t_out = np.zeros(n_time_steps)
Q_dot = np.zeros(n_time_steps)
eff = np.zeros(n_time_steps)
t_in = np.ones(n_time_steps) * (t_in_0+273.15)


def step_down(t, trigger_time):
    if t <= trigger_time:
        return 1
    return 0


glass = GlassCover(n_nodes, 293, n_time_steps)
air = AirGap(n_nodes, 293, n_time_steps)
absorber = Absorber(n_nodes, 293, n_time_steps)
water = Water(n_nodes, 293, n_time_steps)
insulation = Insulation(n_nodes, 293, n_time_steps)


def check_convergence(n_nodes, n_converge, glass, air, absorber):
    for j in range(n_nodes):
        error = np.zeros(5)
        error[0] = glass.error(j)
        error[1] = air.error(j)
        error[2] = absorber.error(j)
        error[3] = water.error(j)
        error[4] = insulation.error(j)
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
        glass.temp, air.temp, absorber.temp, water.temp, insulation.temp,
        t_amb[t], dtau, dz, n_nodes, mdot, w_f)
    n_iter = 0
    while n_converge < 5 * n_nodes:
        n_iter = n_iter + 1
        glass.update_old()
        air.update_old()
        absorber.update_old()
        water.update_old()
        insulation.update_old()
        glass.calc_t(0, t, dtau, t_amb, absorber.temp,
                     air.temp, G_r, B, C, D, E, F)
        air.calc_t(0, glass.temp, absorber.temp, dtau, G, H)
        absorber.calc_t(0, t, dtau, glass.temp, air.temp,
                        water.temp, insulation.temp, G_r, K, L, M, O, P, Q)
        insulation.calc_t(0, t, dtau, absorber.temp, t_amb, V, W, X)
        water.temp[0] = t_in[t]
        for j in range(1, n_nodes):
            glass.calc_t(j, t, dtau, t_amb, absorber.temp,air.temp, G_r, B, C, D, E, F)
            air.calc_t(j, glass.temp, absorber.temp, dtau, G, H)
            absorber.calc_t(j, t, dtau, glass.temp, air.temp,water.temp, insulation.temp, G_r, K, L, M, O, P, Q)
            water.calc_t(j, dtau, dz, absorber.temp, R, S, U)
            insulation.calc_t(j, t, dtau, absorber.temp, t_amb, V, W, X)
        n_converge = check_convergence(n_nodes, n_converge, glass, air, absorber)
    glass.select_node(t, selected_node)
    air.select_node(t, selected_node)
    absorber.select_node(t, selected_node)
    water.select_node(t, selected_node)
    insulation.select_node(t, selected_node)
    t_out[t] = water.temp[n_nodes - 1]
    Q_dot[t] = mdot * n_tubes * 4186 * (t_out[t] - t_in[t])
    eff[t] = G_r[t] and Q_dot[t] / (G_r[t] * collector_area) or 0
print('converged')
runtime = time.time() - start
print("Runtime:", runtime)
plots(time_steps, dtau, glass.selected_node, air.selected_node, absorber.selected_node, water.selected_node,
      insulation.selected_node, t_amb, t_in, t_out, Q_dot, eff, selected_node, n_nodes)
save_to_file(time_steps, dtau, t_in, G_r, absorber.selected_node,
             glass.selected_node, Q_dot, eff)
