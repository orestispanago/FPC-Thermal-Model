import numpy as np

from materials import Air, Water

water = Water()
air = Air()


def get_h(t_water, t_air, t_glass, t_abs, t_insul, n_nodes, t_amb, delta_a, d_in, t, w_f):
    ny_air = air.ny(t_air)
    alpha_air = air.alpha(t_air)
    k_air = air.k(t_air)
    k_water = water.k(t_water)
    Pr_water = water.prandtl(t_water)
    ny_water = water.ny(t_water)
    Re_f = w_f * d_in / ny_water
    segma = 5.6697 * 10**-8
    emi_g = .88
    emi_ab = .1
    emi_i = .05
    g = 9.81
    theta = (np.pi / 4)  # tilt angle
    a = 1.9
    b = .92
    L = a  # collector dimensions
    U_inf = 1.5  # wind velocity
    h_f = np.zeros(n_nodes)
    h_r1 = np.zeros(n_nodes)
    h_c1 = np.zeros(n_nodes)
    h_g_am = np.zeros(n_nodes)
    h_i_am = np.zeros(n_nodes)
    Nu_f = np.zeros(n_nodes)
    Nu_a = np.zeros(n_nodes)
    Ra = np.zeros(n_nodes)
    ny_am = 1.5743 * 10**-5
    k_am = .0262  # ambient thermal conductivity
    Pr_am = 0.71432  # ambiant Prandtl number
    delta = 4 * a * b / np.sqrt(a**2 + b**2)
    Re_am = U_inf * delta / ny_am
    Nu_am = .86 * Re_am**.5 * Pr_am**(1 / 3)
    h_c2 = Nu_am * k_am / delta
    t_sky = .0552 * t_amb**1.5
    Nu_f = 4.4 + (.00398 * (Re_f * Pr_water * (d_in / L))**1.66 /
                  (1 + .0114 * (Re_f * Pr_water * (d_in / L))**1.12))
    h_f = Nu_f * k_water / d_in
    for j in range(n_nodes):
        Ra[j] = abs(t_glass[j] - t_abs[j]) * g * delta_a**3 / \
            (ny_air[j] * alpha_air[j] * t_air[j])
        AA = Ra[j] and 1 - (1708 / (Ra[j] * np.cos(theta))) or 0
        BB = (Ra[j] * np.cos(theta) / 5830)**(1 / 3) - 1
        if AA <= 0:
            if BB <= 0:
                Nu_a[j] = 1
            else:
                Nu_a[j] = 1 + BB
        elif BB <= 0:
            Nu_a[j] = 1 + (1.44 * (1 - (1708 * (np.sin(1.8 * theta))
                           ** 1.6 / (Ra[j] * np.cos(theta)))) * AA)
        else:
            Nu_a[j] = 1 + (1.44 * (1 - (1708 * (np.sin(1.8 * theta))
                           ** 1.6 / (Ra[j] * np.cos(theta)))) * AA) + BB
        h_r1[j] = (segma * (t_abs[j]**2 + t_glass[j]**2) *
                   (t_abs[j] + t_glass[j])) / ((1 / emi_ab) + (1 / emi_g) - 1)
        h_c1[j] = Nu_a[j] * k_air[j] / delta_a
        if t_glass[j] - t_amb[t] == 0:
            h_g_am[j] = h_c2
        else:
            h_g_am[j] = ((segma * emi_g * (t_glass[j]**4 - t_sky[t]**4)
                          ) / (t_glass[j] - t_amb[t])) + h_c2
        if t_insul[j] - t_amb[t] == 0:
            h_i_am[j] = h_c2
        else:
            h_i_am[j] = ((segma * emi_i * (t_insul[j]**4 - t_sky[t]**4)
                          ) / (t_insul[j] - t_amb[t])) + h_c2
    return [h_g_am, h_r1, h_c1, h_f, h_i_am]
