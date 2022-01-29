import random
from os.path import expanduser
from sys import argv

import numpy as np
from matplotlib import pyplot as plt
from scipy.constants import mega

from extradata import ExtraData, CollisionMeta
from utils import filename_from_argv, create_figure, plot_settings, mode_from_fn, scenario_colors

plot_settings()

dotsize = 1.5

angle_label = "Impact angle [deg]"
v_label = "v/v_esc"
time_label = "Time [Myr]"

fig1, ax1 = create_figure()

ax1.set_xlabel(angle_label)
ax1.set_ylabel(v_label)

fig2, ax2 = create_figure()

ax2.set_xlabel(time_label)
ax2.set_ylabel(angle_label)
ax2.set_xscale("log")

fig3, ax3 = create_figure()

ax3.set_xlabel(time_label)
ax3.set_ylabel(v_label)
ax3.set_xscale("log")

fig4, ax4 = create_figure()
ax4.set_xscale("log")
ax4.set_yscale("log")

ax4.set_xlabel(time_label)
ax4.set_ylabel("water loss")

fig5, ax5 = create_figure()
ax5.set_xscale("log")
ax5.set_yscale("log")

ax5.set_xlabel(time_label)
ax5.set_ylabel("mantle loss")

fig6, ax6 = create_figure()
ax6.set_xscale("log")
ax6.set_yscale("log")

ax6.set_xlabel(time_label)
ax6.set_ylabel("core loss")

fig7, ax7 = create_figure()
ax7.set_xlabel(angle_label)

sum = 0
all_angles = []
files = argv[1:]
random.seed(1)
random.shuffle(files)
for file in files:
    fn = filename_from_argv(file)
    print(fn)
    try:
        ed = ExtraData.load(fn)
    except:
        print("skipping")
        continue
    vs = []
    angles = []
    water_loss = []
    mantle_loss = []
    core_loss = []
    times = []

    for collision in ed.tree.get_tree().values():
        sum += 1
        meta: CollisionMeta = collision["meta"]
        vs.append(meta.input.velocity_esc)
        angles.append(meta.input.alpha)
        loss = 1 - meta.water_retention
        if loss == 0:
            loss = 1e-5  # TODO: proper fix for log-log
        water_loss.append(loss)
        mantle_loss.append(1 - meta.mantle_retention)
        core_loss.append(1 - meta.core_retention)
        times.append(meta.time / mega)

    mode = mode_from_fn(fn)

    kwargs = {
        "s": dotsize,
        "color": scenario_colors[mode],
        "alpha": .5
    }
    ax1.scatter(angles, vs, **kwargs)
    ax2.scatter(times, angles, **kwargs)
    ax3.scatter(times, vs, **kwargs)
    if mode != "pm":
        ax4.scatter(times, water_loss, **kwargs)
        ax5.scatter(times, mantle_loss, **kwargs)
        ax6.scatter(times, core_loss, **kwargs)
    all_angles.extend(angles)
ax4.autoscale(enable=True, axis='y')
all_angles = np.asarray(all_angles)

hist, bins = np.histogram(all_angles, bins=50)
width = bins[1] - bins[0]
center = (bins[:-1] + bins[1:]) / 2
ax7.bar(center, hist, align="center", width=width)
xs = np.linspace(0, 90, 1000)
ys = np.sin(np.radians(xs) * 2) * np.sum(np.radians(width) * hist)
ax7.plot(xs, ys, c="C1")

print()
print(all_angles.mean())

for i, fig in enumerate([fig1, fig2, fig3, fig4, fig5, fig6, fig7]):
    fig.tight_layout()
    fig.savefig(expanduser(f"~/tmp/collision{i}.pdf"))

print(sum)
plt.show()
