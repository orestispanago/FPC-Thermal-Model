import numpy as np
from scipy import interpolate


class Water:
    def __init__(self):
        self.temps = np.array([278.15, 298.15, 323.15, 348.15, 368.15])
        self.cps = np.array([4200, 4183, 4181, 4190, 4210])
        self.densities = np.array([1000, 997.1, 988, 974.9, 961.9])
        self.kfs = np.array([0.5576, 0.5948, 0.6305, 0.653, 0.6634])
        # dynamic viscosity (Pa.s)
        self.mus = np.array(
            [0.001519, 0.0008905, 0.0005471, 0.0003779, 0.0002974])
        self.prandtls = np.array([11.44, 6.263, 3.628, 2.425, 1.888])
        self.nys = self.mus/self.densities
        self.rho_spline = interpolate.interp1d(
            self.temps, self.densities, kind='cubic')
        self.cp_spline = interpolate.interp1d(
            self.temps, self.cps, kind='cubic')
        self.kf_spline = interpolate.interp1d(
            self.temps, self.kfs, kind='cubic')
        self.nyf_spline = interpolate.interp1d(
            self.temps, self.nys, kind='cubic')
        self.prandtl_spline = interpolate.interp1d(
            self.temps, self.prandtls, kind='cubic')

    def rho(self, t_f: np.ndarray) -> np.ndarray:
        """Calculates water density (kg/m3)"""
        return self.rho_spline(t_f)

    def cp(self, t_f: np.ndarray) -> np.ndarray:
        """Calculates water heat capacity (J/kg.K)"""
        return self.cp_spline(t_f)

    def k(self, t_f: np.ndarray) -> np.ndarray:
        """Calculates water thermal conductivity (W/m.K)"""
        return self.kf_spline(t_f)

    def ny(self, t_f: np.ndarray) -> np.ndarray:
        """Calculates water kinematic viscosity (m2/s)"""
        return self.nyf_spline(t_f)

    def prandtl(self, t_f: np.ndarray) -> np.ndarray:
        """Calculates water Prandtl number"""
        return self.prandtl_spline(t_f)


class Air:
    def __init__(self):
        self.temps = [250, 300, 350, 400, 450]
        self.densities = [1.4235, 1.1771, 1.0085, 0.88213, 0.8770]
        self.nys = np.array([11.44, 15.89, 20.92, 26.41, 32.39]) * 1e-6
        self.alphas = np.array([15.9, 22.5, 29.9, 38.3, 47.2]) * 1e-6
        self.ks = np.array([22.3, 26.3, 30.0, 33.8, 37.3]) * 1e-3
        self.rho_spline = interpolate.interp1d(
            self.temps, self.densities, kind='cubic')
        self.k_spline = interpolate.interp1d(self.temps, self.ks, kind='cubic')
        self.ny_spline = interpolate.interp1d(
            self.temps, self.nys, kind='cubic')
        self.alpha_spline = interpolate.interp1d(
            self.temps, self.alphas, kind='cubic')

    def rho(self, t_a: np.ndarray) -> np.ndarray:
        """Calculates air density from temperature (kg/m3)"""
        return self.rho_spline(t_a)

    def ny(self, t_a: np.ndarray) -> np.ndarray:
        """Calculates air kinematic viscosity from temperature (m2/s)"""
        return self.ny_spline(t_a)

    def alpha(self, t_a: np.ndarray) -> np.ndarray:
        "Calculates air thermal diffusivity from temperature (m2/s)"
        return self.alpha_spline(t_a)

    def k(self, t_a: np.ndarray) -> np.ndarray:
        "Calculates air thermal conductivity from temperature (W/m.k)"
        return self.k_spline(t_a)
