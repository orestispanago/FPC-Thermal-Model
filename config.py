from constants import length, width


n_nodes = 8                             # number of nodes along the tube
flowrate = 80                           # total flow rate (Liters per hour)
total_simulation_time = 4000            # total running time (s)
t_in_0 = 293                            # initial input temperature (C)
t_amb_0 = 301                           # initial ambient temperature (C)
t_init = 293
area = length*width                     # (m2)
n_tubes = 8                             # number of tubes
d_in = 9.5 / 1000                       # tube inner diameter
flow_per_tube = flowrate / n_tubes      # mass flow rate per tube (GPM)
Vdot = flow_per_tube / 3.6e+6           # volume flow rate per tube (m**3/s)
mdot = Vdot * 1000                      # mass flow rate per tube (Kg/s)
w_f = 4 * Vdot / (3.14 * d_in**2)      # working fluid velocity
dz = length / (n_nodes - 1)             # spatial step (m)
dtau = dz / w_f                         # maximum time step (s)
n_time_steps = round(total_simulation_time / dtau)  # number of time steps


if dtau > dz / w_f:
    print("error in flow rate")
    exit()
