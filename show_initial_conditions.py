from pathlib import Path

import matplotlib
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Colormap, Normalize, LinearSegmentedColormap
from matplotlib.figure import Figure
from mpl_toolkits.axes_grid1 import make_axes_locatable
from rebound import Simulation

from extradata import ExtraData
from water_sim import add_particles_from_conditions_file


def cm_to_inch(cm: float) -> float:
    return cm / 2.54


sim = Simulation()
ed = ExtraData()
add_particles_from_conditions_file(sim, ed, "initcon/conditions_final.input")

size_factor = 2
min_val, max_val = 0.2, 1.0
orig_cmap: Colormap = matplotlib.cm.get_cmap('plasma_r')
colors = orig_cmap(np.linspace(min_val, max_val, 20))
# colors = ["black", "lightblue"]
cmap: Colormap = LinearSegmentedColormap.from_list("mycmap", colors)
a_list = []
e_list = []
m_list = []
wf_list = []
for p in sim.particles[3:]:
    orbit = p.calculate_orbit(primary=sim.particles[0])
    a_list.append(orbit.a)
    e_list.append(orbit.e)
    m_list.append(p.m)
    wf_list.append(ed.pd(p).water_mass_fraction)
m_list = np.asarray(m_list)
with np.errstate(divide='ignore'):  # allow 0 water (becomes -inf)
    color_val = (np.log10(wf_list) + 5) / 5
colors = cmap(color_val)

plt.figure(figsize=(cm_to_inch(16.5), cm_to_inch(6)))
fig: Figure = plt.gcf()
print(fig.get_size_inches())

plt.scatter(a_list, e_list, s=size_factor * m_list / m_list.mean(), c=colors)
plt.xlabel("a [AU]")
plt.ylabel("e")
divider = make_axes_locatable(plt.gca())
cax = divider.append_axes("bottom", size="10%", pad=0.5)

plt.colorbar(ScalarMappable(norm=Normalize(vmin=-5, vmax=0), cmap=cmap),
             label="log(water mass fraction)", cax=cax, orientation="horizontal")

plt.tight_layout()
plt.savefig(Path("~/tmp/out.pdf").expanduser())
plt.show()
