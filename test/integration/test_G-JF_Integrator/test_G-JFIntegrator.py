#!/usr/bin/env python3
import matplotlib.pyplot as plt
import numpy as np
import math

traj_hist, traj_bin_edges = np.histogram(\
        [float(line.split()[1]) for line in filter(lambda x: x[0] == 'X', \
               open("data/test_G-JFIntegrator_position.xyz").readlines())],\
        bins=200, density=True)
traj_x = [(traj_bin_edges[i-1] + traj_bin_edges[i]) / 2.0 \
           for i in range(1, len(traj_bin_edges))]

plt.plot(traj_x, traj_hist, "o", label="numerical")

kB   = 0.0019872041
T    = 300.0
beta = 1.0 / (kB * T)

def potential(x):
    return 0.25 * x**4 + math.sin(1.0 + 5 * x)

xmin = -5.0
xmax =  5.0
Nx   = 10000
dx   = (xmax - xmin) / Nx
xs   = [xmin + (i+0.5) * dx for i in range(0, Nx)]
Z    = sum(math.exp(-beta * potential(x))  for x in xs)
P    = [math.exp(-beta * potential(x)) / Z / dx for x in xs]

plt.xlim(-2.5, 2.5)
plt.plot(xs, P, label="theoretical", color="black")
plt.title("G-JFIntegrator")
plt.legend()
plt.savefig("data/test_G-JFIntegrator.png")
plt.show()
