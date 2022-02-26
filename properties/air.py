import numpy as np
from scipy.interpolate import interp1d

temps = [250, 300, 350, 400, 450]
densities = [1.4235, 1.1771, 1.0085, 0.88213, 0.8770]
nys = np.array([11.44, 15.89, 20.92, 26.41, 32.39]) * 1e-6
alphas = np.array([15.9, 22.5, 29.9, 38.3, 47.2]) * 1e-6
ks = np.array([22.3, 26.3, 30.0, 33.8, 37.3]) * 1e-3
rho_spline = interp1d(temps, densities, kind='cubic')
k_spline = interp1d(temps, ks, kind='cubic')
ny_spline = interp1d(temps, nys, kind='cubic')
alpha_spline = interp1d(temps, alphas, kind='cubic')


def rho(t_a: np.ndarray) -> np.ndarray:
    """Calculates air density from temperature (kg/m3)"""
    return rho_spline(t_a)


def ny(t_a: np.ndarray) -> np.ndarray:
    """Calculates air kinematic viscosity from temperature (m2/s)"""
    return ny_spline(t_a)


def alpha(t_a: np.ndarray) -> np.ndarray:
    "Calculates air thermal diffusivity from temperature (m2/s)"
    return alpha_spline(t_a)


def k(t_a: np.ndarray) -> np.ndarray:
    "Calculates air thermal conductivity from temperature (W/m.k)"
    return k_spline(t_a)
