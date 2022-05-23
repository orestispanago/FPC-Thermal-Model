import time
import numpy as np
from coefficients import coeff
from temp_calculations import (calc_error, calc_t_abs, calc_t_air,
                               calc_t_glass, calc_t_ins, calc_t_w)

from config import n_time_steps, n_nodes, dtau, dz, mdot, w_f, t_amb_0, t_in_0, t_init, trigger_time


start = time.time()

t_amb = t_amb_0
t_in = np.zeros(n_time_steps) + t_in_0

t_gc, t_ac, t_absc, t_wc, t_insc = [np.zeros([n_nodes, n_time_steps]) for _ in range(5)]


def step_down(t, trigger_time):
    return 1 * (t < trigger_time)

time_steps = np.arange(n_time_steps)
G_r = 660*step_down(time_steps*dtau, trigger_time) # Heat flux of solar radiation. (W/sqm)


def check_convergence(n_converge, t_glass, t_glass_old, t_air, t_air_old, 
                      t_abs, t_abs_old, t_water, t_water_old, t_ins, t_ins_old):
    error_g = calc_error(t_glass, t_glass_old)
    error_air = calc_error(t_air, t_air_old)
    error_abs = calc_error(t_abs, t_abs_old)
    error_w = calc_error(t_water, t_water_old)
    error_ins = calc_error(t_ins, t_ins_old)
    for i in [error_g, error_air, error_abs, error_w, error_ins]:
        if np.any(i > 10**-4):
            return n_converge
        else:
            n_converge+=len(t_glass)
    return n_converge

def run():
    t_glass, t_air, t_abs, t_water, t_ins = np.zeros((5, n_nodes))+t_init
    for t in time_steps:
        n_converge = 0
        B, C, D, E, F, G, H, K, L, M, O, P, Q, R, S, U, V, W, X = coeff(
            t_glass, t_air, t_abs, t_water, t_ins,
            t_amb, dtau, dz, mdot, w_f)
        n_iter = 0
        while n_converge < 5 * n_nodes:
            n_iter = n_iter + 1
            t_glass_old = t_glass.copy()
            t_air_old = t_air.copy()
            t_abs_old = t_abs.copy()
            t_water_old = t_water.copy()
            t_ins_old = t_ins.copy()
    
            t_glass = calc_t_glass(dtau, t_glass_old, t_amb, t_abs, t_air, G_r[t], B, C, D, E, F)
            t_air = calc_t_air(t_air_old, t_glass, t_abs, dtau, G, H)
            t_abs = calc_t_abs(dtau, t_abs_old, t_glass, t_air, t_water, t_ins, G_r[t], K, L, M, O, P, Q)
            t_ins = calc_t_ins(dtau, t_ins_old, t_abs, t_amb, V, W, X)
            t_water[0] = t_in[t]
            t_water[1:] = calc_t_w(dtau, dz, t_water_old[1:], t_water[:-1], t_abs[1:], R[1:], S[1:], U[1:])
            n_converge = check_convergence(n_converge, t_glass, t_glass_old, t_air, t_air_old, t_abs, t_abs_old, t_water, t_water_old, t_ins, t_ins_old)
        t_gc[:, t] = t_glass
        t_ac[:, t] = t_air
        t_absc[:,t] = t_abs
        t_wc[:,t] = t_water
        t_insc[:,t] = t_ins
    print(f"Runtime: {(time.time() - start):.3f}")
    
    return [t_gc, t_ac, t_absc, t_wc, t_insc]


