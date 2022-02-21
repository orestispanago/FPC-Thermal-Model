import time

import matplotlib.pyplot as plt
import numpy as np

from coefficients import coeff
from present import plots


# def results(n_nodes, flowrate, total_simulation_time, t_in_0):
n_nodes, flowrate, total_simulation_time, t_in_0 = 18,0.2,4000,25
# n_nodes= number of n_nodes along the tube.
# flowrate= the total flow rate enter the system in GPM
# total_simulation_time= total running time (s).
# t_in_0= initial input temperature (C)
# tank volume in litters.
collector_area = 80 * 39 * 0.0254**2  # (m**2) converted from inches to m
n_tubes = 8  # number of tubes
d_in = 9.5 / 1000  # tube inner diameter
L = 1900 / 1000  # length of tubes
flow_per_tube = flowrate / n_tubes  # fluid volume flow rate per tube (GPM)
Vdot = flow_per_tube / 15852  # fluid volume flow_per_tube rate (m**3/s)
mdot = Vdot * 1000  # fluid mass flow_per_tube rate (Kg/s)
w_f = 4 * Vdot / (np.pi * d_in**2)  # working fluid velocity
dz = L / (n_nodes - 1)  # spatial step (m)
dtau = dz / w_f  # maximum time step (s)
n_time_steps = round(total_simulation_time / dtau) # number of time steps
time_steps = np.arange(n_time_steps)  
selected_node = int(n_nodes / 2)
start = time.time()
# if dtau > dz / w_f:
#     print("error in flow rate")
# else:
with open("ptemp.out", "w") as f:
    f.write(
        'number of n_nodes = {}\nflowrate(GPM) = {}\ntotal_simulation_time(s) = {}\ninitial temperature(K) = {}\ntime step(s) = {}\n'
        .format(n_nodes, flowrate, total_simulation_time, t_in_0, dtau))
    f.write(
        'time \t irrad \t t_int \t T_out \t Q_out \t T_g \t T_a     T_ab    T_f     T_i     n_iter \n'
    )
    t_am = np.zeros(n_time_steps + 1) +(28+273.15)  # Ambient temp.
    G_r = np.zeros(n_time_steps + 1)  # Heat flux of solar radiation. (W/sqm)
    t_g = np.ones(n_nodes) * 293  # initial glass temp.
    t_a = np.ones(n_nodes) * 293  # initial air gap temp.
    t_ab = np.ones(n_nodes) * 293  # initial absorber temp.
    t_f = np.ones(n_nodes) * 293  # initial fluid temp.
    t_i = np.ones(n_nodes) * 293  # initial insulation temp.
    t_gc = np.zeros(n_time_steps)
    t_ac = np.zeros(n_time_steps)
    t_abc = np.zeros(n_time_steps)
    t_fc = np.zeros(n_time_steps)
    t_ic = np.zeros(n_time_steps)
    t_out = np.zeros(n_time_steps)
    Q_dot = np.zeros(n_time_steps)
    eff = np.zeros(n_time_steps)
    t_in = np.ones(n_time_steps) * (t_in_0+273.15)
    counter = 1
    for t in range(n_time_steps + 1):
        # Step change in irradiance at 1 hour
        if t * dtau <= 3600:
            G_r[t] = 660
        else:
            G_r[t] = 0
    for t in  range(n_time_steps):
        n_converge = 0
        B, C, D, E, F, G, H, K, L, M, O, P, Q, R, S, U, V, W, X, J = coeff(
            t_g, t_a, t_ab, t_f, t_i, t_am, dtau, dz, n_nodes, mdot, t, w_f)
        n_iter = 0
        while n_converge < 5 * n_nodes:
            n_iter = n_iter + 1
            t_g_old = t_g.copy()
            t_a_old = t_a.copy()
            t_ab_old = t_ab.copy()
            t_f_old = t_f.copy()
            t_i_old = t_i.copy()
            t_g[0] = ((t_g_old[0] / dtau) + (B[0] * t_am[t]) + (C[0] * t_ab[0]) + (D[0] * t_a[0]) + (E[0] * G_r[t])) / F[0]
            t_a[0] = ((t_a_old[0] / dtau) + (G[0] * (t_g[0] + t_ab[0]))) / H[0]
            t_ab[0] = ((t_ab_old[0] / dtau) + (K[0] * G_r[t]) + (L[0] * t_g[0]) + (M[0] * t_a[0]) + (O[0] * t_f[0]) + (P[0] * t_i[0])) / Q[0]
            t_f[0] = t_in[t]
            t_i[0] = ((t_i_old[0] / dtau) + (V[0] * t_ab[0]) + (W[0] * t_am[t])) / X[0]
            for j in range(1, n_nodes):
                t_g[j] = ((t_g_old[j] / dtau) + (B[j] * t_am[t]) + (C[j] * t_ab[j]) + (D[j] * t_a[j]) + (E[j] * G_r[t])) / F[j]
                t_a[j] = ((t_a_old[j] / dtau) + (G[j] *  (t_g[j] + t_ab[j]))) / H[j]
                t_ab[j] = ((t_ab_old[j] / dtau) + (K[j] * G_r[t]) + (L[j] * t_g[j]) + (M[j] * t_a[j]) + (O[j] * t_f[j]) + (P[j] * t_i[j])) / Q[j]
                t_f[j] = ((t_f_old[j] / dtau) + (R[j] * t_ab[j]) + (S[j] * t_f[j - 1] / dz)) / U[j]
                t_i[j] = ((t_i_old[j] / dtau) + (V[j] * t_ab[j]) + (W[j] * t_am[t])) / X[j]
            # check convergence
            ccc = 0
            for j in range(n_nodes):
                if ccc <= 0:
                    error = np.zeros(5)
                    error[0] = abs(t_g[j] - t_g_old[j]) / t_g[j]
                    error[1] = abs(t_a[j] - t_a_old[j]) / t_a[j]
                    error[2] = abs(t_ab[j] - t_ab_old[j]) / t_ab[j]
                    error[3] = abs(t_f[j] - t_f_old[j]) / t_f[j]
                    error[4] = abs(t_i[j] - t_i_old[j]) / t_i[j]
                    for i in range(5):
                        if error[i] <= 10**-4:
                            n_converge = n_converge + 1
                        else:
                            ccc = 1
        t_gc[t] = t_g[selected_node]
        t_ac[t] = t_a[selected_node]
        t_abc[t] = t_ab[selected_node]
        t_fc[t] = t_f[selected_node]
        t_ic[t] = t_i[selected_node]
        t_out[t] = t_f[n_nodes - 1]
        Q_dot[t] = mdot * n_tubes * 4186 * (t_out[t] - t_in[t])
        eff[t] = Q_dot[t] / (G_r[t] * collector_area)
        time_minutes = dtau * t / 60
        print(
            f'time = {time_minutes:.1f} t_in = {t_in[t]:.2f} T_out = {t_out[t]:.2f} Q_dot = {Q_dot[t]:.2f} T_g ={t_gc[t]:.2f} T_a = {t_ac[t]:.2f} T_ab = {t_abc[t]:.2f} T_f = {t_fc[t]:.2f} T_i = {t_ic[t]:.2f} n_iter ={n_iter}'
        )
        if time_minutes-counter >=0:
            f.write(
               f'{time_minutes:.1f}\t\t{G_r[t]:.1f}\t{t_in[t]:.2f}\t{t_out[t]:.2f}\t{Q_dot[t]:.2f}\t{t_gc[t]:.2f}\t{t_ac[t]:.2f}\t{t_abc[t]:.2f}\t{t_fc[t]:.2f}\t{t_ic[t]:.2f}\t\t{n_iter}\n'
            )
            counter=counter+1
    print('converged')
    plots(time_steps, dtau, t_gc, t_ac, t_abc, t_fc, t_ic, t_am, t_in,t_out, Q_dot,eff)

    runtime = time.time() - start
    # count = fprintf(fid, 'run time = #6.1f\n', runtime)
    # status = fclose(fid)
    print("Runtime:", runtime)

# results(18, 0.2, 38, 25)
