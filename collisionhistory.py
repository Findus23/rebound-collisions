from matplotlib import pyplot as plt
from rebound import SimulationArchive, Simulation

from extradata import ExtraData, CollisionMeta
from utils import filename_from_argv

fn = filename_from_argv()
ed = ExtraData.load(fn.with_suffix(".extra.json"))
sa = SimulationArchive(str(fn.with_suffix(".bin")))

last_sim: Simulation = sa[-1]
print([p.hash.value for p in last_sim.particles])
print(last_sim.t)
for particle in last_sim.particles:
    if ed.pd(particle).type in ["sun", "gas giant"]:
        continue
    masses = []
    objects = []
    times = []
    hash = particle.hash.value

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
    if len(times) < 2:
        continue
    masses = [p.total_mass for p in objects]
    wmfs = [p.water_mass_fraction for p in objects]
    plt.step(times, wmfs, label=particle.hash.value)
plt.xscale("log")
plt.yscale("log")
plt.legend()
plt.show()
