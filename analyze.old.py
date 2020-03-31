import json
import matplotlib.pyplot as plt

import numpy as np
from rebound import SimulationArchive, OrbitPlot

sa = SimulationArchive("water.bin")
with open("water.meta.json") as f:
    meta = json.load(f)
    tmax = meta["tmax"]
    savesteps = meta["savesteps"]
    max_n = meta["max_n"]

times = np.linspace(0., tmax, savesteps)
data = np.zeros((3, savesteps * max_n))
j = 0

for i, t in enumerate(times):
    sim = sa.getSimulation(t=t)
    for particle in range(1, sim.N):
        data[0, j] = i
        data[1, j] = particle
        data[2, j] = sim.particles[particle].a
        j += 1
plt.scatter(data[0, :], data[2, :], c=data[1, :], cmap="tab10", s=.4)
plt.colorbar()
# OrbitPlot(sim, slices=1)
plt.show()
