import numpy as np
from scipy import interpolate

def rho(t_air) -> np.float64:
    """Calculates air density from temperature"""
    temperatures = [250, 300, 350, 400, 450]
    densities = [1.4235, 1.1771, 1.0085, 0.88213, 0.8770]
    spline = interpolate.interp1d(temperatures, densities, kind='cubic')
    return spline(t_air)


def kf(t_f):
    """Calculates fluid conductivity, kinematic viscosity and Prandtl number from temperature"""
    # https://onlinelibrary.wiley.com/doi/pdf/10.1002/9781118534892.app2
    # kinematic viscosity here: https://wiki.anton-paar.com/en/water/
    temperatures = np.array([278.15, 298.15, 323.15, 348.15, 368.15])
    densities = np.array([1000, 997.1, 988, 974.9, 961.9])
    mu = np.array([0.001519, 0.0008905, 0.0005471, 0.0003779, 0.0002974]) # dynamic viscosity (Pa.s)
    kf_vec = np.array([0.5576, 0.5948,0.6305, 0.653, 0.6634])
    Prf_vec = np.array([11.44, 6.263, 3.628, 2.425, 1.888])
    nyf_vec = mu/densities
    kf_spline = interpolate.interp1d(temperatures, kf_vec, kind='cubic')
    nyf_spline = interpolate.interp1d(temperatures, nyf_vec, kind='cubic')
    Prf_spline = interpolate.interp1d(temperatures, Prf_vec, kind='cubic')
    k_f = kf_spline(t_f)  # fluid thermal conductivity (W/m2.K)
    ny_f = nyf_spline(t_f)  # fluid kinematic viscosity (m2/s)
    Pr_f = Prf_spline(t_f)  # fluid Prandtl number
    return [k_f, ny_f, Pr_f]


def air_prop(t_a):
    """Calculates air Nussfelt number, thermal diffusivity, and conductivity from temperature"""
    temperatures = [250, 300, 350, 400, 450]
    ny_vec = np.array([11.44, 15.89, 20.92, 26.41, 32.39]) * 1e-6
    alpha_vec = np.array([15.9, 22.5, 29.9, 38.3, 47.2]) * 1e-6
    ka_vec = np.array([22.3, 26.3, 30.0, 33.8, 37.3]) * 1e-3
    k_a_spline = interpolate.interp1d(temperatures, ka_vec, kind='cubic')
    ny_a_spline = interpolate.interp1d(temperatures, ny_vec, kind='cubic')
    alpha_a_spline = interpolate.interp1d(temperatures, alpha_vec, kind='cubic')
    try:
        ny_a = ny_a_spline(t_a)
    except ValueError as e:
        print("ERROR:",e)
        print(t_a)
    alpha_a = alpha_a_spline(t_a)
    k_a = k_a_spline(t_a)
    return [ny_a, alpha_a, k_a]


def waterprop(t_f):
    """Calculates fluid density and heat capacity from temperature"""
    temperatures = np.array([278.15, 298.15, 323.15, 348.15, 368.15])
    cp_f_vec = np.array([4200,4183,4181,4190,4210])
    densities = np.array([1000, 997.1, 988, 974.9, 961.9])
    density_spline = interpolate.interp1d(temperatures, densities, kind='cubic')
    cp_f_spline = interpolate.interp1d(temperatures, cp_f_vec, kind='cubic')
    rho_f = density_spline(t_f)
    c_f = cp_f_spline(t_f)
    return [rho_f, c_f]
