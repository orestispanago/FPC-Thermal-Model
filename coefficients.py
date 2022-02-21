from constants import *
from get_h import get_h
from properties import rho,waterprop
import numpy as np

def coeff(t_g, t_a, t_ab, t_f, t_i, t_am, dtau, dz, nodes, mdot, k, w_f):
    """Coefficients of the transient temperature equations."""
    h_g_am, h_r1, h_c1, h_f, h_i_am = get_h(t_f, t_a, t_g, t_ab, t_i, nodes,
                                            t_am, delta_a, d_in, k, w_f)
    rho_a = rho(t_a)
    rho_f, c_f = waterprop(t_f)
    B = np.zeros(nodes)
    C = np.zeros(nodes)
    D = np.zeros(nodes)
    E = np.zeros(nodes)
    F = np.zeros(nodes)
    J = np.zeros(nodes)
    K = np.zeros(nodes)
    L = np.zeros(nodes)
    M = np.zeros(nodes)
    O = np.zeros(nodes)
    P = np.zeros(nodes)
    G = np.zeros(nodes)
    H = np.zeros(nodes)
    Q = np.zeros(nodes)
    R = np.zeros(nodes)
    S = np.zeros(nodes)
    U = np.zeros(nodes)
    V = np.zeros(nodes)
    W = np.zeros(nodes)
    X = np.zeros(nodes)
    for j in range(nodes):
        B[j] = h_g_am[j] / (c_g * rho_g * delta_g)
        C[j] = h_r1[j] / (c_g * rho_g * delta_g)
        D[j] = h_c1[j] / (c_g * rho_g * delta_g)
        E[j] = alpha / (c_g * rho_g * delta_g)
        F[j] = (1 / dtau) + B[j] + C[j] + D[j]
        J[j] = c_ab * rho_ab * (p * delta_ab + np.pi * (r_o**2 - r_in**2))
        K[j] = p * (tau_alpha) / J[j]
        L[j] = h_r1[j] * p / J[j]
        M[j] = h_c1[j] * p / J[j]
        O[j] = np.pi * d_in * h_f[j] / J[j]
        P[j] = p * K_i / (J[j] * delta_i)
        G[j] = h_c1[j] * p / (c_a * rho_a[j] * (p * delta_a - np.pi * r_o**2))
        H[j] = (1 / dtau) + (2 * G[j])
        Q[j] = (1 / dtau) + L[j] + M[j] + O[j] + P[j]
        R[j] = np.pi * d_in * h_f[j] / (c_f[j] * rho_f[j] * A)
        S[j] = mdot / (rho_f[j] * A)
        U[j] = (1 / dtau) + R[j] + (S[j] / dz)
        V[j] = 2 * K_i / (c_i * rho_i * delta_i**2)
        W[j] = 2 * h_i_am[j] / (c_i * rho_i * delta_i)
        X[j] = (1 / dtau) + V[j] + W[j]
    return [B, C, D, E, F, G, H, K, L, M, O, P, Q, R, S, U, V, W, X, J]
