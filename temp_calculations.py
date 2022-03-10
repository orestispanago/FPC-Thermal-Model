def calc_error(temp, temp_old):
    return abs(temp - temp_old) / temp


def calc_t_glass(dtau, t_glass_old, t_amb, t_abs, t_air, G_r, B, C, D, E, F):
    return ((t_glass_old / dtau) + (B * t_amb) + (C * t_abs) + (D * t_air) + (E * G_r)) / F


def calc_t_air(t_air_old, t_glass, t_abs, dtau, G, H):
    return ((t_air_old / dtau) + (G * (t_glass + t_abs))) / H


def calc_t_abs(dtau, t_abs_old, t_glass, t_air, t_water, t_insul, G_r, K, L, M, O, P, Q):
    return ((t_abs_old / dtau) + (K * G_r) + (L * t_glass) + (M * t_air) + (O * t_water) + (P * t_insul)) / Q


def calc_t_w(dtau, dz, t_w_old, t_w, t_abs, R, S, U):
    return ((t_w_old / dtau) + (R * t_abs) + (S * t_w / dz)) / U


def calc_t_ins(dtau, t_ins_old, t_abs, t_amb, V, W, X):
    return ((t_ins_old / dtau) + (V * t_abs) + (W * t_amb)) / X
