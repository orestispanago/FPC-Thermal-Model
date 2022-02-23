import numpy as np

from constants import (alpha, cp_abs, cp_air, cp_glass, cp_insul, d_in,
                       delta_a, delta_abs, delta_glass, delta_insul, flow_area,
                       k_insul, p, r_in, r_o, rho_abs, rho_glass, rho_insul,
                       tau_alpha)
from heat_transfer_coefficients import get_h
from materials import Air, Water

air = Air()
water = Water()


def coeff(t_glass, t_air, t_abs, t_water, t_insul, t_amb, dtau, dz, n_nodes, mdot, t, w_f):
    """Coefficients of the transient temperature equations."""
    h_glass_amb, h_r1, h_c1, h_water, h_insul_amb = get_h(t_water, t_air, t_glass, t_abs, t_insul, n_nodes,
                                                          t_amb, delta_a, d_in, t, w_f)
    rho_air = air.rho(t_air)
    rho_water = water.rho(t_water)
    cp_water = water.cp(t_water)
    B = np.zeros((n_nodes))
    C = np.zeros(n_nodes)
    D = np.zeros(n_nodes)
    E = np.zeros(n_nodes)
    F = np.zeros(n_nodes)
    J = np.zeros(n_nodes)
    K = np.zeros(n_nodes)
    L = np.zeros(n_nodes)
    M = np.zeros(n_nodes)
    O = np.zeros(n_nodes)
    P = np.zeros(n_nodes)
    G = np.zeros(n_nodes)
    H = np.zeros(n_nodes)
    Q = np.zeros(n_nodes)
    R = np.zeros(n_nodes)
    S = np.zeros(n_nodes)
    U = np.zeros(n_nodes)
    V = np.zeros(n_nodes)
    W = np.zeros(n_nodes)
    X = np.zeros(n_nodes)
    for j in range(n_nodes):
        B[j] = h_glass_amb[j] / (cp_glass * rho_glass * delta_glass)
        C[j] = h_r1[j] / (cp_glass * rho_glass * delta_glass)
        D[j] = h_c1[j] / (cp_glass * rho_glass * delta_glass)
        E[j] = alpha / (cp_glass * rho_glass * delta_glass)
        F[j] = (1 / dtau) + B[j] + C[j] + D[j]
        J[j] = cp_abs * rho_abs * (p * delta_abs + np.pi * (r_o**2 - r_in**2))
        K[j] = p * (tau_alpha) / J[j]
        L[j] = h_r1[j] * p / J[j]
        M[j] = h_c1[j] * p / J[j]
        O[j] = np.pi * d_in * h_water[j] / J[j]
        P[j] = p * k_insul / (J[j] * delta_insul)
        G[j] = h_c1[j] * p / (cp_air * rho_air[j] * (p * delta_a - np.pi * r_o**2))
        H[j] = (1 / dtau) + (2 * G[j])
        Q[j] = (1 / dtau) + L[j] + M[j] + O[j] + P[j]
        R[j] = np.pi * d_in * h_water[j] / (cp_water[j] * rho_water[j] * flow_area)
        S[j] = mdot / (rho_water[j] * flow_area)
        U[j] = (1 / dtau) + R[j] + (S[j] / dz)
        V[j] = 2 * k_insul / (cp_insul * rho_insul * delta_insul**2)
        W[j] = 2 * h_insul_amb[j] / (cp_insul * rho_insul * delta_insul)
        X[j] = (1 / dtau) + V[j] + W[j]
    return [B, C, D, E, F, G, H, K, L, M, O, P, Q, R, S, U, V, W, X]
