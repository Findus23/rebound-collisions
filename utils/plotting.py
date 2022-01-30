import re
from pathlib import Path
from typing import Tuple

import matplotlib
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.axes import Axes
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Colormap, LinearSegmentedColormap, Normalize
from matplotlib.figure import Figure
from mpl_toolkits.axes_grid1 import make_axes_locatable

christoph_cmap: Colormap = LinearSegmentedColormap.from_list("mycmap", ["#FF4500", "blue"])

scenario_colors = {
    "rbf": "C3",
    "nn": "C1",
    "lz": "C2",
    "pm": "C0",
}


def create_figure() -> Tuple[Figure, Axes]:
    """
    helper function for matplotlib OOP interface with proper typing
    """
    fig: Figure = plt.figure()
    ax: Axes = fig.gca()
    return fig, ax


def plot_settings() -> None:
    ...
    # plt.style.use("dark_background")


def get_water_cmap() -> Colormap:
    # return christoph_cmap
    min_val, max_val = 0.2, 1.0

    orig_cmap: Colormap = matplotlib.cm.get_cmap('plasma_r')
    colors = orig_cmap(np.linspace(min_val, max_val, 20))
    cmap: Colormap = LinearSegmentedColormap.from_list("mycmap", colors)
    return cmap


def add_water_colormap(fig: Figure, ax: Axes, cmap: Colormap) -> None:
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("bottom", size="10%", pad=0.5)

    fig.colorbar(ScalarMappable(norm=Normalize(vmin=-5, vmax=0), cmap=cmap),
                 label="log(water mass fraction)", cax=cax, orientation="horizontal")


def add_au_e_label(ax: Axes) -> None:
    ax.set_xlabel("semi-major axis [AU]")
    ax.set_ylabel("eccentricity")


def mode_from_fn(fn: Path) -> str:
    return re.search(r"final_(\w+)_", str(fn)).group(1)
