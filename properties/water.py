import numpy as np
from scipy.interpolate import interp1d

temps = np.array([278.15, 298.15, 323.15, 348.15, 368.15])
cps = np.array([4200, 4183, 4181, 4190, 4210])
densities = np.array([1000, 997.1, 988, 974.9, 961.9])
kfs = np.array([0.5576, 0.5948, 0.6305, 0.653, 0.6634])
# dynamic viscosity (Pa.s)
mus = np.array([0.001519, 0.0008905, 0.0005471, 0.0003779, 0.0002974])
prandtls = np.array([11.44, 6.263, 3.628, 2.425, 1.888])
nys = mus/densities
rho_spline = interp1d(temps, densities, kind='cubic')
cp_spline = interp1d(temps, cps, kind='cubic')
kf_spline = interp1d(temps, kfs, kind='cubic')
nyf_spline = interp1d(temps, nys, kind='cubic')
prandtl_spline = interp1d(temps, prandtls, kind='cubic')


def rho(t_f: np.ndarray) -> np.ndarray:
    """Calculates water density (kg/m3)"""
    return rho_spline(t_f)


def cp(t_f: np.ndarray) -> np.ndarray:
    """Calculates water heat capacity (J/kg.K)"""
    return cp_spline(t_f)


def k(t_f: np.ndarray) -> np.ndarray:
    """Calculates water thermal conductivity (W/m.K)"""
    return kf_spline(t_f)


def _ny(t_f: np.ndarray) -> np.ndarray:
    """Calculates water kinematic viscosity (m2/s)"""
    return nyf_spline(t_f)


def _prandtl(t_f: np.ndarray) -> np.ndarray:
    """Calculates water Prandtl number"""
    return prandtl_spline(t_f)


def _reynolds(t_f: np.ndarray, d_in: float, w_f: float) -> np.ndarray:
    """Calculates water Reynolds number"""
    ny = _ny(t_f)
    return w_f * d_in / ny


def nusselt(t_f, d_in, w_f, length):
    """Calculates water Nusselt number"""
    prandtl = _prandtl(t_f)
    reynolds = _reynolds(t_f, d_in, w_f)
    return 4.4 + (.00398 * (reynolds * prandtl * (d_in / length))**1.66 / (1 + .0114 * (reynolds * prandtl * (d_in / length))**1.12))


def h(t_in, d_in, w_f, length):
    """Calculates air thermal conductivity (W/m.k)"""
    return nusselt(t_in, d_in, w_f, length) * k(t_in) / d_in
