import matplotlib.pyplot as plt
import numpy as np
from rebound import SimulationArchive, Particle

from extradata import ExtraData
from utils import filename_from_argv

fn = filename_from_argv()
sa = SimulationArchive(str(fn.with_suffix(".bin")))
ed = ExtraData.load(fn.with_suffix(".extra.json"))
print(ed.meta)
times = np.linspace(0., ed.meta.current_time, ed.meta.current_steps)

data = {}
for i, t in enumerate(times):
    try:
        sim = sa.getSimulation(t=t)
    except ValueError:
        print("save is incomplete")
        break
    for pn in range(1, sim.N):
        part: Particle = sim.particles[pn]
        hash = part.hash.value
        if hash not in data:
            data[hash] = ([], [])
        data[hash][0].append(t)
        data[hash][1].append(part.a)

for name, d in data.items():
    times, values = d
    print(list(map(len, [times, values])))
    if False:
        plt.scatter(times, values, label=name, s=.9)
    else:
        plt.plot(times, values, label=name)
# plt.legend()
# OrbitPlot(sim, slices=1)
plt.show()
