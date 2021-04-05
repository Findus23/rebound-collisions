from sys import argv

from matplotlib import pyplot as plt
from matplotlib.axes import Axes
from matplotlib.figure import Figure
from rebound import SimulationArchive, Simulation

from extradata import ExtraData
from utils import filename_from_argv

fig: Figure = plt.figure()
ax: Axes = plt.gca()

for file in argv[1:]:
    fn = filename_from_argv(file)
    print(fn)
    try:
        ed = ExtraData.load(fn)
        sa = SimulationArchive(str(fn.with_suffix(".bin")))
    except:
        print("skipping")
        continue
    ts = []
    Ns = []
    sim: Simulation
    for sim in sa:
        num_planetesimals = sum(ed.pd(p).type == "planetesimal" for p in sim.particles)
        num_embryos = sum(ed.pd(p).type == "embryo" for p in sim.particles)
        N = num_embryos + num_planetesimals
        Ns.append(N)
        ts.append(sim.t)

    ax.step(ts, Ns, label=fn, where="post")
ax.set_xlabel("time [yr]")
ax.set_ylabel("number of objects")
# ax.set_xscale("log")

plt.legend()
plt.show()
