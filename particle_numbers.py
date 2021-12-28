import pickle
import random
import re
from os.path import expanduser
from pathlib import Path
from random import shuffle
from sys import argv

from matplotlib import pyplot as plt
from matplotlib.axes import Axes
from matplotlib.figure import Figure
from mpl_toolkits.axes_grid1.inset_locator import mark_inset, inset_axes
from rebound import SimulationArchive, Simulation

from extradata import ExtraData
from utils import filename_from_argv, plot_settings, is_ci

cache_file = Path("particle_numbers_cache.pickle.xz")
if cache_file.exists():
    with cache_file.open("rb") as f:
        cache = pickle.load(f)
else:
    cache = {}

plot_settings()

fig: Figure = plt.figure()
ax: Axes = plt.gca()
# axins2 = zoomed_inset_axes(ax, zoom=1, loc="center")
axins2 = inset_axes(ax, 2,2, loc="center right")  # no zoom
random.seed(1)

colors = {
    "rbf": "C0",
    "nn": "C1",
    "lz": "C2",
    "pm": "C3",
}
files = argv[1:]
shuffle(files)
for file in files:
    fn = filename_from_argv(file)
    print(fn)
    if "bak" in str(fn):
        continue
    if fn.name in cache:
        Ns, ts = cache[fn.name]
    else:
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
        cache[fn.name] = Ns, ts
    mode = re.search(r"final_(\w+)_", str(fn)).group(1)
    for a in [ax, axins2]:
        a.step(ts, Ns, label=mode, where="post", color=colors[mode], linewidth=0.7)
    # break

with cache_file.open("wb") as f:
    pickle.dump(cache, f)

ax.set_xlabel("time [yr]")
ax.set_ylabel("number of objects")
# ax.set_xscale("log")

axins2.set_xlim(1e8, 2e8)
axins2.set_ylim(0, 15)
axins2.tick_params(labelleft=False, labelbottom=False)
mark_inset(ax, axins2, loc1=2, loc2=4, fc="none", ec="0.5")

handles, labels = ax.get_legend_handles_labels()
by_label = dict(zip(labels, handles))
ax.legend(by_label.values(), by_label.keys())
# plt.tight_layout()
if not is_ci():
    plt.savefig(expanduser(f"~/tmp/particle_numbers.pdf"))
plt.show()
