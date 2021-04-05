from pathlib import Path
from sys import argv

from matplotlib import pyplot as plt
from scipy.constants import mega

from extradata import ExtraData, CollisionMeta
from utils import filename_from_argv, create_figure

dotsize = 1.5

angle_label = "Impact angle (deg)"
v_label = "v/v_esc"
time_label = "Time (Myr)"

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

for file in argv[1:]:
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
    times = []

    for collision in ed.tree.get_tree().values():
        meta: CollisionMeta = collision["meta"]
        vs.append(meta.input.velocity_esc)
        angles.append(meta.input.alpha)
        loss = 1 - meta.water_retention
        if loss == 0:
            loss = 1e-5  # TODO: proper fix for log-log
        water_loss.append(loss)
        times.append(meta.time / mega)

    ax1.scatter(angles, vs,s=dotsize)
    ax2.scatter(times, angles,s=dotsize)
    ax3.scatter(times, vs,s=dotsize)
    ax4.scatter(times, water_loss,s=dotsize)

ax4.autoscale(enable=True, axis='y')

for i, fig in enumerate([fig1, fig2, fig3, fig4]):
    fig.savefig(Path("plots") / fn.with_suffix(f".collision{i}.pdf").name)

plt.show()
