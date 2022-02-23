import numpy as np


class GlassCover:
    def __init__(self, n_nodes, t_glass_0) -> None:
        self.t_glass = np.ones(n_nodes)*t_glass_0
        self.t_glass_old = np.zeros(n_nodes)

    def update_old(self):
        self.t_glass_old = self.t_glass.copy()

    def calc_t(self, i, t, dtau, t_amb, t_abs, t_air, G_r, B, C, D, E, F):
        self.t_glass[i] = ((self.t_glass_old[i] / dtau) + (B[i] * t_amb[t]) +
                           (C[i] * t_abs[i]) + (D[i] * t_air[i]) + (E[i] * G_r[t])) / F[i]

    def check_converge(self, j):
        return abs(self.t_glass[j] - self.t_glass_old[j]) / self.t_glass[j]
