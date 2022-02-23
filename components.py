import numpy as np


class Component:
    def __init__(self, n_nodes, temp_0, n_time_steps):
        self.temp = np.ones(n_nodes)*temp_0
        self.temp_old = np.zeros(n_nodes)
        self.selected_node = np.zeros(n_time_steps)

    def update_old(self):
        self.temp_old = self.temp.copy()

    def error(self, j):
        return abs(self.temp[j] - self.temp_old[j]) / self.temp[j]

    def select_node(self, t, node):
        self.selected_node[t] = self.temp[node]


class GlassCover(Component):
    def __init__(self, n_nodes, temp_0, n_time_steps):
        super().__init__(n_nodes, temp_0, n_time_steps)

    def calc_t(self, i, t, dtau, t_amb, t_abs, t_air, G_r, B, C, D, E, F):
        self.temp[i] = ((self.temp_old[i] / dtau) + (B[i] * t_amb[t]) +
                        (C[i] * t_abs[i]) + (D[i] * t_air[i]) + (E[i] * G_r[t])) / F[i]


class AirGap(Component):
    def __init__(self, n_nodes, temp_0, n_time_steps):
        super().__init__(n_nodes, temp_0, n_time_steps)

    def calc_t(self, i, t_glass, t_abs, dtau, G, H):
        self.temp[i] = ((self.temp_old[i] / dtau) +
                        (G[i] * (t_glass[i] + t_abs[i]))) / H[i]


class Absorber(Component):
    def __init__(self, n_nodes, temp_0, n_time_steps):
        super().__init__(n_nodes, temp_0, n_time_steps)

    def calc_t(self, i, t, dtau, t_glass, t_air, t_water, t_insul, G_r, K, L, M, O, P, Q):
        self.temp[i] = ((self.temp_old[i] / dtau) + (K[i] * G_r[t]) + (L[i] * t_glass[i]) +
                        (M[i] * t_air[i]) + (O[i] * t_water[i]) + (P[i] * t_insul[i])) / Q[i]


class Water(Component):
    def __init__(self, n_nodes, temp_0, n_time_steps):
        super().__init__(n_nodes, temp_0, n_time_steps)

    def calc_t(self, i, dtau, dz, t_abs, R, S, U):
        self.temp[i] = ((self.temp_old[i] / dtau) + (R[i] *
                        t_abs[i]) + (S[i] * self.temp[i - 1] / dz)) / U[i]


class Insulation(Component):
    def __init__(self, n_nodes, temp_0, n_time_steps):
        super().__init__(n_nodes, temp_0, n_time_steps)

    def calc_t(self, i, t, dtau, t_abs, t_amb, V, W, X):
        self.temp[i] = ((self.temp_old[i] / dtau) +
                        (V[i] * t_abs[i]) + (W[i] * t_amb[t])) / X[i]
