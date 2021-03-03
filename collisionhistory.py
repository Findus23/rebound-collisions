from matplotlib import pyplot as plt
from matplotlib.axes import Axes
from matplotlib.figure import Figure
from rebound import SimulationArchive, Simulation

from extradata import ExtraData, CollisionMeta
from utils import filename_from_argv

fn = filename_from_argv()
ed = ExtraData.load(fn)
sa = SimulationArchive(str(fn.with_suffix(".bin")))

last_sim: Simulation = sa[-1]
print([p.hash.value for p in last_sim.particles])
print(last_sim.t)

fig: Figure = plt.figure()
ax_masses: Axes = fig.add_subplot(2,1, 1)
ax_wmfs: Axes = fig.add_subplot(2,1, 2)

for particle in last_sim.particles:
    if ed.pd(particle).type in ["sun", "gas giant"]:
        continue
    masses = []
    objects = []
    times = []
    hash = particle.hash.value
    objects.append(ed.pdata[hash])
    times.append(ed.meta.current_time)

    while True:
        print(f"looking at {hash}")
        try:
            collision = ed.tree.get_tree()[hash]
        except KeyError:
            print("found end of the tree")
            break
        meta: CollisionMeta = collision["meta"]
        parents = collision["parents"]
        masses.append(ed.pdata[hash].total_mass)
        objects.append(ed.pdata[hash])
        times.append(meta.time)
        # print(collision)
        if ed.pdata[parents[0]].total_mass > ed.pdata[parents[1]].total_mass:
            hash = parents[0]
        else:
            hash = parents[1]
    objects.append(ed.pdata[hash])
    times.append(0)  # TODO: check log-x
    if len(times) < 3:
        continue
    masses = [p.total_mass for p in objects]
    wmfs = [p.water_mass_fraction for p in objects]
    figs = []
    ax_masses.step(times, masses, label=particle.hash.value)
    ax_wmfs.step(times, wmfs, label=particle.hash.value)
    ax_wmfs.set_ylabel("water mass fraction")
    ax_masses.set_ylabel("masses [kg]")
    for ax in [ax_wmfs, ax_masses]:
        ax.set_xlabel("time [yr]")
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.legend()
fig.tight_layout()
plt.show()
