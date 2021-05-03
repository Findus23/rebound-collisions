from sys import argv

from matplotlib import pyplot as plt
from matplotlib.axes import Axes
from matplotlib.figure import Figure
from rebound import SimulationArchive, Simulation

from extradata import ExtraData
from utils import filename_from_argv, plot_settings, is_ci

plot_settings()

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
    perfect_merging = "pm" in str(fn)
    ax.step(ts, Ns, label=fn, where="post", linestyle="dashed" if perfect_merging else "solid", linewidth=0.7)
ax.set_xlabel("time [yr]")
ax.set_ylabel("number of objects")
# ax.set_xscale("log")

plt.legend()
plt.tight_layout()
if not is_ci():
    plt.savefig("/home/lukas/tmp/particle_numbers.pdf", transparent=True)
plt.show()
