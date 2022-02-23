import numpy as np

from constants import (alpha, cp_abs, cp_air, cp_glass, cp_insul, d_in,
                       delta_a, delta_abs, delta_glass, delta_insul, flow_area,
                       k_insul, p, r_in, r_o, rho_abs, rho_glass, rho_insul,
                       tau_alpha)
from heat_transfer_coefficients import get_h
from materials import Air, Water

air = Air()
water = Water()


# def coeff_glass(dtau, h_glass_amb, h_r1, h_c1):
#     B = h_glass_amb / (cp_glass * rho_glass * delta_glass)
#     C = h_r1 / (cp_glass * rho_glass * delta_glass)
#     D = h_c1 / (cp_glass * rho_glass * delta_glass)
#     E = alpha / (cp_glass * rho_glass * delta_glass)
#     F = (1 / dtau) + B + C + D
#     return [B, C, D, E, F]


# def coeff_air(h_c1, rho_air, dtau):
#     G = h_c1 * p / (cp_air * rho_air * (p * delta_a - np.pi * r_o**2))
#     H = (1 / dtau) + (2 * G)
#     return [G, H]


# def coeff_abs(h_r1, h_c1, h_water, dtau):
#     J = cp_abs * rho_abs * (p * delta_abs + np.pi * (r_o**2 - r_in**2))
#     K = p * (tau_alpha) / J
#     L = h_r1 * p / J
#     M = h_c1 * p / J
#     O = np.pi * d_in * h_water / J
#     P = p * k_insul / (J * delta_insul)
#     Q = (1 / dtau) + L + M + O + P
#     return [J, K, L, M, O, P, Q]


# def coeff_water(h_water, cp_water, rho_water, mdot, dtau, dz):
#     R = np.pi * d_in * h_water / (cp_water * rho_water * flow_area)
#     S = mdot / (rho_water * flow_area)
#     U = (1 / dtau) + R + (S / dz)
#     return [R, S, U]


# def coeff_ins(h_insul_amb, dtau):
#     V = 2 * k_insul / (cp_insul * rho_insul * delta_insul**2)
#     W = 2 * h_insul_amb / (cp_insul * rho_insul * delta_insul)
#     X = (1 / dtau) + V + W
#     return [V, W, X]


def coeff(t_glass, t_air, t_abs, t_water, t_insul, t_amb, dtau, dz, n_nodes, mdot, t, w_f):
    """Coefficients of the transient temperature equations."""
    h_glass_amb, h_r1, h_c1, h_water, h_insul_amb = get_h(t_water, t_air, t_glass, t_abs, t_insul, n_nodes,
                                                          t_amb, delta_a, d_in, t, w_f)
    rho_air = air.rho(t_air)
    rho_water = water.rho(t_water)
    cp_water = water.cp(t_water)
    """Glass"""
    B = h_glass_amb / (cp_glass * rho_glass * delta_glass)
    C = h_r1 / (cp_glass * rho_glass * delta_glass)
    D = h_c1 / (cp_glass * rho_glass * delta_glass)
    E = alpha / (cp_glass * rho_glass * delta_glass)
    F = (1 / dtau) + B + C + D
    """Air"""
    G = h_c1 * p / (cp_air * rho_air * (p * delta_a - np.pi * r_o**2))
    H = (1 / dtau) + (2 * G)
    """Absorber"""
    J = cp_abs * rho_abs * (p * delta_abs + np.pi * (r_o**2 - r_in**2))
    K = p * (tau_alpha) / J
    L = h_r1 * p / J
    M = h_c1 * p / J
    O = np.pi * d_in * h_water / J
    P = p * k_insul / (J * delta_insul)
    Q = (1 / dtau) + L + M + O + P
    """Water"""
    R = np.pi * d_in * h_water / (cp_water * rho_water * flow_area)
    S = mdot / (rho_water * flow_area)
    U = (1 / dtau) + R + (S / dz)
    """Insulation"""
    V = 2 * k_insul / (cp_insul * rho_insul * delta_insul**2)
    W = 2 * h_insul_amb / (cp_insul * rho_insul * delta_insul)
    X = (1 / dtau) + V + W
    return [B, C, D, E, F, G, H, K, L, M, O, P, Q, R, S, U, V, W, X]
