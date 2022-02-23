import numpy as np

from constants import emi_abs, emi_glass, emi_insul, g, length, segma, width
from materials import Air, Ambient, Water

water = Water()
air = Air()
ambient = Ambient(length, width)  # TODO maybe create collector class?


def h_external(j, t, t_surface, t_amb, h_c2, segma, emissivity, t_sky):
    """Calculates heat transfer coefficient of external surface (glass or insulation)"""
    if t_surface[j] - t_amb[t] == 0:
        return h_c2
    return ((segma * emissivity * (t_surface[j]**4 - t_sky[t] ** 4)) / (t_surface[j] - t_amb[t])) + h_c2


def h_radiation_1(j, t_abs, t_glass):
    """Calculates heat transfer coefficient of radiation between absorber and glass"""
    return (segma * (t_abs[j]**2 + t_glass[j]**2) * (t_abs[j] + t_glass[j])) / ((1 / emi_abs) + (1 / emi_glass) - 1)


def h_free_convection_1(j, Nu_a, k_air, delta_a):
    """Calculates heat transfer coefficient of 
    free convection in the inclined gap between absorber and glass"""
    return Nu_a * k_air[j] / delta_a


def rayleigh(j, t_glass, t_abs, delta_a, ny_air, alpha_air, t_air):
    """Calculates Rayleigh number"""
    return abs(t_glass[j] - t_abs[j]) * g * delta_a**3 / (ny_air[j] * alpha_air[j] * t_air[j])


def nu_air(Ra, theta):
    """Calculates Nusselt number for air in gap"""
    AA = Ra and 1 - (1708 / (Ra * np.cos(theta))) or 0
    BB = (Ra * np.cos(theta) / 5830)**(1 / 3) - 1
    if AA <= 0:
        if BB <= 0:
            return 1
        else:
            return 1 + BB
    elif BB <= 0:
        return 1 + (1.44 * (1 - (1708 * (np.sin(1.8 * theta)) ** 1.6 / (Ra * np.cos(theta)))) * AA)
    else:
        return 1 + (1.44 * (1 - (1708 * (np.sin(1.8 * theta)) ** 1.6 / (Ra * np.cos(theta)))) * AA) + BB


def get_h(t_water, t_air, t_glass, t_abs, t_insul, n_nodes, t_amb, delta_a, d_in, t, w_f):
    ny_air = air.ny(t_air)
    alpha_air = air.alpha(t_air)
    k_air = air.k(t_air)
    theta = (np.pi / 4)  # tilt angle
    h_r1 = np.zeros(n_nodes)
    h_c1 = np.zeros(n_nodes)
    h_g_am = np.zeros(n_nodes)
    h_i_am = np.zeros(n_nodes)
    h_water = water.h(t_water, d_in, w_f, length)
    h_c2 = ambient.h()
    t_sky = 0.0552 * t_amb**1.5
    for j in range(n_nodes):
        Ra = rayleigh(j, t_glass, t_abs, delta_a, ny_air, alpha_air, t_air)
        Nu_a = nu_air(Ra, theta)
        h_r1[j] = h_radiation_1(j, t_abs, t_glass)
        h_c1[j] = h_free_convection_1(j, Nu_a, k_air, delta_a)
        h_g_am[j] = h_external(j, t, t_glass, t_amb,
                               h_c2, segma, emi_glass, t_sky)
        h_i_am[j] = h_external(j, t, t_insul, t_amb,
                               h_c2, segma, emi_insul, t_sky)
    return [h_g_am, h_r1, h_c1, h_water, h_i_am]
