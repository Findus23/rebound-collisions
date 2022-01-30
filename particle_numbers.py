import pickle
import random
from os.path import expanduser
from pathlib import Path
from sys import argv

from matplotlib import pyplot as plt
from matplotlib.axes import Axes
from matplotlib.figure import Figure
from rebound import SimulationArchive, Simulation

from extradata import ExtraData
from utils import filename_from_argv, plot_settings, is_ci, scenario_colors, mode_from_fn

cache_file = Path("particle_numbers_cache.pickle.xz")
if cache_file.exists():
    with cache_file.open("rb") as f:
        cache = pickle.load(f)
else:
    cache = {}

plot_settings()

fig: Figure = plt.figure()
ax: Axes = plt.gca()

files = argv[1:]
random.seed(1)
random.shuffle(files)
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
    mode = mode_from_fn(fn)
    ax.step(ts, Ns, label=mode, where="post", color=scenario_colors[mode], linewidth=0.7,alpha=.5)

with cache_file.open("wb") as f:
    pickle.dump(cache, f)

ax.set_xlabel("time [yr]")
ax.set_ylabel("number of objects")
ax.set_xscale("log")

handles, labels = ax.get_legend_handles_labels()
by_label = dict(zip(labels, handles))
new_labels = {
    "RBF": by_label["rbf"],
    "NN": by_label["nn"],
    "LZ": by_label["lz"],
    "PM": by_label["pm"],
}
ax.legend(new_labels.values(), new_labels.keys())
plt.tight_layout()
if not is_ci():
    plt.savefig(expanduser(f"~/tmp/particle_numbers.pdf"))
plt.show()
