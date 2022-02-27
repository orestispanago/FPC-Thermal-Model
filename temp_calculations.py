def calc_error(j, temp, temp_old):
    return abs(temp[j] - temp_old[j]) / temp[j]


def calc_t_glass(i, t, dtau, t_glass_old, t_amb, t_abs, t_air, G_r, B, C, D, E, F):
    return ((t_glass_old[i] / dtau) + (B[i] * t_amb[t]) + (C[i] * t_abs[i]) + (D[i] * t_air[i]) + (E * G_r[t])) / F[i]


def calc_t_air(i, t_air_old, t_glass, t_abs, dtau, G, H):
    return ((t_air_old[i] / dtau) + (G[i] * (t_glass[i] + t_abs[i]))) / H[i]


def calc_t_abs(i, t, dtau, t_abs_old, t_glass, t_air, t_water, t_insul, G_r, K, L, M, O, P, Q):
    return ((t_abs_old[i] / dtau) + (K * G_r[t]) + (L[i] * t_glass[i]) + (M[i] * t_air[i]) + (O[i] * t_water[i]) + (P * t_insul[i])) / Q[i]


def calc_t_w(i, dtau, dz, t_w_old, t_w, t_abs, R, S, U):
    return ((t_w_old[i] / dtau) + (R[i] * t_abs[i]) + (S[i] * t_w[i - 1] / dz)) / U[i]


def calc_t_ins(i, t, dtau, t_ins_old, t_abs, t_amb, V, W, X):
    return ((t_ins_old[i] / dtau) + (V * t_abs[i]) + (W[i] * t_amb[t])) / X[i]
