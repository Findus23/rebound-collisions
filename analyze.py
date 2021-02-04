import matplotlib.pyplot as plt
from rebound import SimulationArchive, Particle, Simulation

from extradata import ExtraData
from utils import filename_from_argv

fn = filename_from_argv()
sa = SimulationArchive(str(fn.with_suffix(".bin")))
ed = ExtraData.load(fn)
print(ed.meta)

data = {}
sim: Simulation
print(f"{len(sa)} Snapshots found")
for sim in sa:
    t = sim.t
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
