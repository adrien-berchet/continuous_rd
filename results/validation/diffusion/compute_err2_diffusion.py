import numpy as np
import pandas as pd
import os

Du = Dv = 1.125
Dw = 13.5

x0 = 64
y0 = 32
sigma = 5

comps = [("u", Du), ("v", Dv), ("w", Dw)]
errors = []

for file in os.listdir("results"):
    if file.endswith(".dat"):
        data = pd.read_csv(
            os.path.join("results", file),
            skiprows=3,
            names=["x", "y", "u", "v", "w", "P_sin_theta"])

        # Compute time from file name
        t = float(os.path.splitext(file)[0])

        # Compute distance of each cell to the gaussian center
        data["r"] = np.sqrt(np.power(data["x"] - x0, 2) + np.power(data["y"] - y0, 2))

        # Compute the theoretical value
        for comp in comps:
            data[comp[0] + "_theo"] = (
                2 * np.pi * np.power(sigma, 2)
            ) / (
                4 * np.pi * comp[1] * (t + np.power(sigma, 2) / (2 * comp[1]))
            ) * np.exp(
                -np.power(data["r"], 2) / (
                    4 * comp[1] * (t + np.power(sigma, 2) / (2 * comp[1]))))

        # Compute L2 error
        tmp = [t]
        for comp in comps:
            tmp.append(np.sqrt(np.power(data[comp[0]] - data[comp[0] + "_theo"], 2).sum() / np.power(data[comp[0]], 2).sum()))
        errors.append(tmp)

filename = "errors.dat"

err_data = pd.DataFrame(errors, columns=["t", "err2_u", "err2_v", "err2_w"])
err_data.sort_values("t").to_csv(filename, index=False)
