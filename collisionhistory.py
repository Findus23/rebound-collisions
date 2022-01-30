import random
from sys import argv

from matplotlib import pyplot as plt
from matplotlib.axes import Axes
from matplotlib.figure import Figure
from rebound import SimulationArchive, Simulation
from scipy.constants import mega

from extradata import ExtraData, CollisionMeta
from utils import filename_from_argv, earth_mass, earth_water_mass, plot_settings, is_ci, is_potentially_habitable

files = argv[1:]
multifile = len(files) > 1

plot_settings()
fig: Figure = plt.figure()
ax_masses: Axes = fig.add_subplot(2, 1, 1)
ax_wmfs: Axes = fig.add_subplot(2, 1, 2)

ax_wmfs.set_ylabel("water mass fraction")
ax_masses.set_ylabel("masses [kg]")
for ax in [ax_wmfs, ax_masses]:
    ax.set_xlim(1e4, 200 * mega)
    ax.set_xlabel("time [yr]")
    ax.set_xscale("log")
    ax.set_yscale("log")

earth_mass_ax = ax_masses.secondary_yaxis("right",
                                          functions=(lambda x: x / earth_mass, lambda x: x * earth_mass))
earth_mass_ax.set_ylabel('masses [$M_\\oplus$]')
ax_wmfs.axhline(earth_water_mass / earth_mass, linestyle="dotted")

num_formed_planets = 0
num_large_planets = 0
num_habitable_planets = 0
num_water_rich_planets = 0

random.seed(1)
random.shuffle(files)
for file in files:
    fn = filename_from_argv(file)

    ed = ExtraData.load(fn)
    sa = SimulationArchive(str(fn.with_suffix(".bin")))

    last_sim: Simulation = sa[-1]
    print([p.hash.value for p in last_sim.particles])
    print(last_sim.t)

    for particle in last_sim.particles:
        if ed.pd(particle).type in ["sun", "gas giant"]:
            continue
        # if not is_potentially_habitable(particle):
        #     continue
        masses = []
        objects = []
        times = []
        hash = particle.hash.value
        objects.append(ed.pdata[hash])
        times.append(ed.meta.current_time)

        num_formed_planets += 1
        if particle.m > .6 * earth_mass:
            num_large_planets += 1
            if is_potentially_habitable(particle):
                num_habitable_planets += 1
                if ed.pd(particle).water_mass_fraction > 1e-4:
                    num_water_rich_planets += 1

        while True:
            print(f"looking at {hash}")
            try:
                collision = ed.tree.get_tree()[hash]
            except KeyError:
                print("found end of the tree")
                break
            meta: CollisionMeta = collision["meta"]
            parents = collision["parents"]
            print("mass:", ed.pdata[hash].total_mass / earth_mass)
            masses.append(ed.pdata[hash].total_mass)
            objects.append(ed.pdata[hash])
            times.append(meta.time)
            # print(collision)
            if ed.pdata[parents[0]].total_mass > ed.pdata[parents[1]].total_mass:
                hash = parents[0]
            else:
                hash = parents[1]
        objects.append(ed.pdata[hash])
        times.append(0)
        if len(times) < 3:
            continue
        masses = [p.total_mass for p in objects]
        wmfs = [p.water_mass_fraction for p in objects]
        figs = []
        if multifile:
            args = {
                "linewidth": 1,
                "color": "black",
                "alpha": .3
            }
        else:
            args = {}
        ax_masses.step(times, masses, label=particle.hash.value, **args)
        ax_wmfs.step(times, wmfs, label=particle.hash.value, **args)

if not multifile:
    for ax in [ax_wmfs, ax_masses]:
        ax.legend()

print(num_large_planets / num_formed_planets)
habitable_large_planet_fraction = num_habitable_planets / num_large_planets
water_rich_planet_fraction = num_water_rich_planets / num_habitable_planets
print(habitable_large_planet_fraction)
print(water_rich_planet_fraction)

fig.tight_layout()
if not is_ci():
    fig.savefig(f"/home/lukas/tmp/collisionhistory_{fn.name}.pdf", transparent=True)
plt.show()
