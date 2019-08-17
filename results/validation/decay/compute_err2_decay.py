import numpy as np
import pandas as pd
import os

cu = 0.02
cv = 0.025
cw = 0.06

x0 = 64
y0 = 32
sigma = 5

comps = [("u", cu), ("v", cv), ("w", cw)]
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
        data["r"] = np.sqrt(
            np.power(data["x"] - x0, 2) + np.power(data["y"] - y0, 2))

        # Compute the theoretical value
        for comp in comps:
            data[comp[0] + "_theo"] = np.exp(
                -np.power(data["r"], 2) / (2 * np.power(sigma, 2))
            ) * np.exp(-comp[1] * t)

        # Compute L2 error
        tmp = [t]
        for comp in comps:
            tmp.append(np.sqrt(
                np.power(data[comp[0]] - data[comp[0] + "_theo"], 2).sum()
                / np.power(data[comp[0]], 2).sum()))
        errors.append(tmp)

filename = "errors_decay.dat"

err_data = pd.DataFrame(errors, columns=["t", "err2_u", "err2_v", "err2_w"])
err_data.sort_values("t").to_csv(filename, index=False)
