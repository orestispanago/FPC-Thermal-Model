import pandas as pd


def save_to_file(time_steps, dtau, t_in, G_r, t_abs_node, t_glass_node, Q_dot, eff):
    df = pd.DataFrame({"time": time_steps*dtau,
                       "t_in": t_in,
                       "irrad": G_r[:-1],
                       "t_abs": t_abs_node,
                       "t_glass": t_glass_node,
                       "Q_dot": Q_dot,
                       "Eff": eff})
    df.to_csv("output.csv", index=False)
