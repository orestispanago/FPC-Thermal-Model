import numpy as np

from constants import length, width

ny = 1.5743 * 10**-5  # kinematic viscosity (m2/s)
U_inf = 1.5  # wind velocity
k = 0.0262  # thermal conductivity (W/m.K)
prandtl = 0.71432  # Prandtl number
# collector characteristic length
delta = 4 * length * width / np.sqrt(length**2 + width**2)
reynolds = U_inf * delta / ny
nusselt = 0.86 * reynolds**.5 * prandtl ** (1 / 3)
h = nusselt * k / delta
