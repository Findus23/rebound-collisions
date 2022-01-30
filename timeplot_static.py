from os.path import expanduser

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.axes import Axes
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize
from matplotlib.figure import Figure
from rebound import SimulationArchive
from scipy.constants import mega

from extradata import ExtraData
from utils import filename_from_argv, get_water_cmap

mean_mass = 5.208403167890638e+24  # constant between plots
size_mult = 100


def plot_file(file, time, ax: Axes, mode):
    fn = filename_from_argv(file)
    sa = SimulationArchive(str(fn.with_suffix(".bin")))
    ed = ExtraData.load(fn)
    sim = sa.getSimulation(t=time)

    water_fractions = [ed.pd(p).water_mass_fraction for p in sim.particles[1:]]
    a = [p.a for p in sim.particles[1:]]
    e = [p.e for p in sim.particles[1:]]
    m = np.array([p.m for p in sim.particles[1:]])
    # m[:2] /= 1e2 # reduce size of gas giants
    sizes = (np.array(m) / mean_mass) ** (2 / 3) * size_mult
    with np.errstate(divide='ignore'):  # allow 0 water (becomes -inf)
        color_val = (np.log10(water_fractions) + 5) / 5
    cmap = get_water_cmap()
    colors = cmap(color_val)
    ax.set_xlim(0, 5.5)
    ax.set_ylim(-0.05, .8)
    ax.text(0.05, 0.93, f"{time / mega:.0f} Myr", size=10,
            horizontalalignment='left', verticalalignment='top', transform=ax.transAxes)
    ax.text(0.95, 0.93, mode, size=10,
            horizontalalignment='right', verticalalignment='top', transform=ax.transAxes)
    ax.scatter(a, e, s=sizes, zorder=10, c=colors)


def main():
    fig: Figure
    modes = ["RBF", "NN", "PM", "LZ"]
    fig, axes = plt.subplots(5, 4, sharex=True, sharey=True, figsize=(10, 13))
    for col_nr, file in enumerate(
            ["data/final_rbf_1.bin", "data/final_nn_1.bin", "data/final_pm_1.bin", "data/final_lz_1.bin"]):
        for row_nr, time in enumerate([1 * mega, 5 * mega, 20 * mega, 50 * mega, 100 * mega]):
            plot_file(file, time, axes[row_nr][col_nr], mode=modes[col_nr])
    fig_colorbar = plt.figure()
    fig_colorbar.colorbar(ScalarMappable(norm=Normalize(vmin=-5, vmax=0), cmap=get_water_cmap()), aspect=20,
                          label="log(water mass fraction)", orientation="horizontal")
    fig.supxlabel("semi-major axis [AU]")
    fig.supylabel("eccentricity")
    fig.tight_layout()
    fig_colorbar.gca().remove()
    fig.savefig(expanduser(f"~/tmp/timeplot.pdf"), dpi=300)
    fig_colorbar.savefig(expanduser(f"~/tmp/timeplot_colorbar.pdf"), bbox_inches='tight')
    plt.show()


if __name__ == '__main__':
    main()
