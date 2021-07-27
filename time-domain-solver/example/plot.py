# %%
import matplotlib.pyplot as plt
import numpy as np
import os

plt.style.use("classic")

hbar = 1.0545718e-34
e0 = 1.602e-19
c0 = 299792458
eps0 = 8.854e-12

amp = e0 ** 2 / (16 * np.pi ** 3 * c0 * eps0)

time_data = np.loadtxt(f"spectralangular.dat").T

thetamin = 1 / 5.0
thetamax = -1 / 5.0
wmin = 0
wmax = 500
extent = wmin, wmax, thetamin, thetamax


plt.title(r"$\gamma = 5$", fontsize=20)
plt.xlabel(r"$\omega/\omega_0$", fontsize=20)
plt.ylabel(r"$\theta$ (rad)", fontsize=20)
plt.imshow(
    time_data,
    extent=extent,
    aspect="auto",
    cmap="jet",
)
plt.colorbar()
plt.savefig(
    "spectralangular.pdf",
    box_inches="tight",
)

