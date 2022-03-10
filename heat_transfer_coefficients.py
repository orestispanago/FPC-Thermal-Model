import numpy as np

from constants import emi_abs, emi_glass, emi_insul, g, length, segma
from properties import air, ambient, water


def h_external(t_surface, t_amb, h_c2, segma, emissivity):
    """Calculates heat transfer coefficient of external surface (glass or insulation)"""
    t_sky = 0.0552 * t_amb**1.5
    if abs(t_surface - t_amb) <= 0.1:
        return h_c2
    return segma * emissivity * (t_surface**4 - t_sky ** 4) / (t_surface - t_amb) + h_c2


def h_radiation(t_abs, t_glass):
    """Calculates heat transfer coefficient of radiation between absorber and glass"""
    return (segma * (t_abs**2 + t_glass**2) * (t_abs + t_glass)) / ((1 / emi_abs) + (1 / emi_glass) - 1)


def h_free_convection(Nu_a, k_air, delta_a):
    """Calculates heat transfer coefficient of 
    free convection in the inclined gap between absorber and glass"""
    return Nu_a * k_air / delta_a


def rayleigh(t_glass, t_abs, delta_a, ny_air, alpha_air, t_air):
    """Calculates Rayleigh number"""
    return abs(t_glass - t_abs) * g * delta_a**3 / (ny_air * alpha_air * t_air)


def nu_air(Ra, theta):
    """Calculates Nusselt number for air in gap"""
    AA = Ra and 1 - (1708 / (Ra * np.cos(theta))) or 0
    BB = (Ra * np.cos(theta) / 5830)**(1 / 3) - 1
    if AA <= 0:
        if BB <= 0:
            return 1
        return 1 + BB
    if BB <= 0:
        return 1 + (1.44 * (1 - (1708 * (np.sin(1.8 * theta)) ** 1.6 / (Ra * np.cos(theta)))) * AA)
    return 1 + (1.44 * (1 - (1708 * (np.sin(1.8 * theta)) ** 1.6 / (Ra * np.cos(theta)))) * AA) + BB


h_external_v = np.vectorize(h_external)
nu_air_v = np.vectorize(nu_air)


def get_h(t_water, t_air, t_glass, t_abs, t_insul, t_amb, delta_a, d_in, w_f):
    ny_air = air.ny(t_air)
    alpha_air = air.alpha(t_air)
    k_air = air.k(t_air)
    theta = (np.pi / 4)  # tilt angle
    h_water = water.h(t_water, d_in, w_f, length)
    h_c2 = ambient.h
    Ra = rayleigh(t_glass, t_abs, delta_a, ny_air, alpha_air, t_air)
    h_g_am = h_external_v(t_glass, t_amb, h_c2, segma, emi_glass)
    h_i_am = h_external_v(t_insul, t_amb, h_c2, segma, emi_insul)
    Nu_a = nu_air_v(Ra, theta)
    h_r1 = h_radiation(t_abs, t_glass)
    h_c1 = h_free_convection(Nu_a, k_air, delta_a)
    return [h_g_am, h_r1, h_c1, h_water, h_i_am]
