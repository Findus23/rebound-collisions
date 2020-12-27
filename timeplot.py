import matplotlib
import matplotlib.animation as animation
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.cm import ScalarMappable
from matplotlib.collections import PathCollection
from matplotlib.colors import Normalize, Colormap
from matplotlib.text import Text
from rebound import SimulationArchive, Particle

from extradata import ExtraData, ParticleData
from utils import filename_from_argv

matplotlib.use('Qt4Agg')
plt.style.use('dark_background')
fps = 5
duration = 10  # s
total_frames = fps * duration
logtime = False
cmap: Colormap = matplotlib.cm.get_cmap('Blues')


def update_line(num: int, sa: SimulationArchive, ed: ExtraData, dots: PathCollection, title: Text):
    # line = num * int((len(data) / total_frames))
    # print(line, num)
    if logtime:
        timestep = total_frames

    else:
        timestep = ed.meta.current_time / total_frames
        # timestep = 10e6 / total_frames
        time = num * timestep
    print(time)
    sim = sa.getSimulation(t=time)
    title.set_text(f"({len(sim.particles)}) {time / 1e6:.2f}M Years")

    p: Particle
    a = [p.a for p in sim.particles[1:]]
    e = [p.e for p in sim.particles[1:]]

    water_fractions = []
    for p in sim.particles[3:]:
        # try:
        pd: ParticleData = ed.pdata[p.hash.value]
        wf = pd.water_mass_fraction
        # except KeyError:  # gas planet
        #     print(p.hash.value)
        #     wf = 0
        water_fractions.append(wf)
    # a, e, i, M, M_rat = data[line]
    # title.set_text(f"({len(a)}) {ages[line]:.2f}K Years")
    bla = np.array([a, e])
    dots.set_offsets(bla.T)
    water_fractions = np.array(water_fractions)
    color_val = (np.log10(water_fractions) + 5) / 5
    print(color_val)
    colors = cmap(color_val)
    print(colors)
    dots.set_color(colors)
    plt.savefig("tmp/" + str(num) + ".pdf",transparent=True)
    return dots, title


fig1 = plt.figure()

l: PathCollection = plt.scatter([1], [1])

fn = filename_from_argv()
sa = SimulationArchive(str(fn.with_suffix(".bin")))
ed = ExtraData.load(fn.with_suffix(".extra.json"))

plt.xlim(0, 10)
plt.xlabel("a")
title: Text = plt.title("0")
plt.ylim(-0.1, 1)  # e
plt.ylabel("e")

# plt.ylim(0, 90)  # i
# plt.ylabel("i")

# plt.yscale('log')
# plt.ylim(1e-7,1e-3)

# plt.ylim(-0.05, 0.2)  # i
# plt.ylabel("water_fraction")

fig1.colorbar(ScalarMappable(norm=Normalize(vmin=-5, vmax=0), cmap=cmap), label="log(water fraction)")

plt.tight_layout()

line_ani = animation.FuncAnimation(fig1, update_line, total_frames, fargs=(sa, ed, l, title),
                                   interval=1000 / fps, repeat=False)
line_ani.save(str(fn.with_suffix(".mp4")), dpi=300)
# plt.show()
