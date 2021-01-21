from matplotlib import pyplot as plt
from scipy.constants import mega

from extradata import ExtraData, CollisionMeta
from utils import filename_from_argv, create_figure

fn = filename_from_argv()
ed = ExtraData.load(fn.with_suffix(".extra.json"))

vs = []
angles = []
times = []

for collision in ed.tree.get_tree().values():
    meta: CollisionMeta = collision["meta"]
    vs.append(meta.input.velocity_esc)
    angles.append(meta.input.alpha)
    times.append(meta.time / mega)

angle_label = "Impact angle (deg)"
v_label = "v/v_esc"
time_label = "Time (Myr)"

fig1, ax1 = create_figure()

ax1.scatter(angles, vs)
ax1.set_xlabel(angle_label)
ax1.set_ylabel(v_label)

fig2, ax2 = create_figure()

ax2.scatter(times, angles)
ax2.set_xlabel(time_label)
ax2.set_ylabel(angle_label)
ax2.set_xscale("log")

fig3, ax3 = create_figure()

ax3.scatter(times, vs)
ax3.set_xlabel(time_label)
ax3.set_ylabel(v_label)
ax3.set_xscale("log")


plt.show()
