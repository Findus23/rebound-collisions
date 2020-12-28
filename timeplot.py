import argparse
from collections import namedtuple

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


class MyProgramArgs(argparse.Namespace):
    save_video: bool
    fps: int
    duration: int
    y_axis: str


matplotlib.use("Qt5Agg")
plt.style.use("dark_background")
cmap: Colormap = matplotlib.cm.get_cmap('Blues')
# fps = 5
# duration = 10  # s
logtime = False  # TODO: support logarithmic time


def update_line(num: int, args: MyProgramArgs, sa: SimulationArchive, ed: ExtraData, dots: PathCollection, title: Text):
    # line = num * int((len(data) / total_frames))
    total_frames = args.fps * args.duration
    # print(line, num)
    if logtime:
        timestep = total_frames
    else:
        timestep = ed.meta.current_time / total_frames
        # timestep = 10e6 / total_frames
        time = num * timestep
    print(time)
    sim = sa.getSimulation(t=time)
    if time < 1e3:
        timestr = f"{time:.0f}"
    elif time < 1e6:
        timestr = f"{time / 1e3:.2f}K"
    else:
        timestr = f"{time / 1e6:.2f}M"
    title.set_text(f"({len(sim.particles)}) {timestr} Years")

    p: Particle

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
    a = [p.a for p in sim.particles[1:]]

    if args.y_axis == "e":
        bla = np.array([a, [p.e for p in sim.particles[1:]]])
    elif args.y_axis == "i":
        bla = np.array([a, [p.inc for p in sim.particles[1:]]])
    elif args.y_axis == "Omega":
        bla = np.array([a, [p.Omega for p in sim.particles[1:]]])
    else:
        raise ValueError("invalid y-axis")
    dots.set_offsets(bla.T)
    water_fractions = np.array(water_fractions)
    color_val = (np.log10(water_fractions) + 5) / 5
    print(color_val)
    colors = cmap(color_val)
    print(colors)
    dots.set_color(colors)
    # plt.savefig("tmp/" + str(num) + ".pdf",transparent=True)
    return dots, title


def main(args: MyProgramArgs):
    total_frames = args.fps * args.duration

    logtime = False
    cmap: Colormap = matplotlib.cm.get_cmap("Blues")

    fig1 = plt.figure()

    l: PathCollection = plt.scatter([1], [1])

    fn = filename_from_argv()
    sa = SimulationArchive(str(fn.with_suffix(".bin")))
    ed = ExtraData.load(fn.with_suffix(".extra.json"))

    plt.xlim(0, 10)
    plt.xlabel("a")
    title: Text = plt.title("0")

    if args.y_axis == "e":
        plt.ylim(-0.1, 1)  # e
        plt.ylabel("e")
    elif args.y_axis == "i":
        plt.ylim(0, 90)  # i
        plt.ylabel("i")
    elif args.y_axis == "Omega":
        plt.ylim(0, 360)  # i
        plt.ylabel("Omega")

    # plt.yscale("log")
    # plt.ylim(1e-7,1e-3)

    # plt.ylim(-0.05, 0.2)  # i
    # plt.ylabel("water_fraction")

    fig1.colorbar(ScalarMappable(norm=Normalize(vmin=-5, vmax=0), cmap=cmap), label="log(water fraction)")

    plt.tight_layout()
    line_ani = animation.FuncAnimation(fig1, update_line, total_frames, fargs=(args, sa, ed, l, title),
                                       interval=1000 / args.fps, repeat=False)
    if args.save_video:
        line_ani.save(str(fn.with_suffix(".mp4")), dpi=300)
    else:
        plt.show()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="create a video showing ",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("file")
    parser.add_argument("-v", "--save-video", action="store_true",
                        help="save video as .mp4 file")
    parser.add_argument("--fps", default=5, type=int,
                        help="frames per second")
    parser.add_argument("--duration", default=10, type=int,
                        help="duration in seconds")
    parser.add_argument("--y-axis", default="e", type=str, choices=["e", "i", "Omega"],
                        help="what to show on the y-axis")
    ArgNamespace = namedtuple('ArgNamespace', ['some_arg', 'another_arg'])
    args = parser.parse_args()
    print(vars(args))
    print(args.save_video)
    # noinspection PyTypeChecker
    main(args)
