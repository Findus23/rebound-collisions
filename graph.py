from os import environ

import numpy as np
from graphviz import Digraph
from matplotlib.colors import to_hex

from extradata import ExtraData, CollisionMeta
from utils import filename_from_argv, get_water_cmap

fn = filename_from_argv()
ed = ExtraData.load(fn)

dot = Digraph(comment='Collisions')

# dot.engine = "neato"

if dot.engine == "neato":
    dot.attr("graph", overlap="false")
interacting_objects = set()
for merged, originals in ed.tree.get_tree().items():
    first_parent = True
    for parent in originals["parents"]:
        meta: CollisionMeta = originals["meta"]
        water_ret = meta.water_retention
        mantle_ret = meta.mantle_retention
        core_ret = meta.core_retention
        if first_parent:

            label = f"{water_ret:.2f}/{mantle_ret:.2f}/{core_ret:.2f}"
            first_parent = False
            dot.node(label=label, name=f"{merged}-collision", shape="diamond")
            dot.edge(f"{merged}-collision", str(merged))
        else:
            label = None
        dot.edge(str(parent), f"{merged}-collision")
        interacting_objects.add(parent)
        interacting_objects.add(int(merged))

cmap = get_water_cmap()

for name in ed.pdata.keys():
    object = ed.pdata[name]
    if object.type == "sun":
        displayname = f"{name} (Sun)"
    elif object.type == "gas giant":
        displayname = f"{name} (gas giant)"
    else:
        displayname = str(name)
    try:
        mass = object.total_mass
    except KeyError:
        mass = 0
    with np.errstate(divide='ignore'):  # allow 0 water (becomes -inf)
        color_val = (np.log10(object.water_mass_fraction) + 5) / 5
    color = cmap(color_val)
    if color_val > 0.55:
        textcolor = "white"
    else:
        textcolor = "black"
    print(to_hex(color), color_val)
    dot.node(name=str(name), label=f"m={mass:.1e}\nwmf={object.water_mass_fraction:.1e}",
             shape='box' if object.type == "planetesimal" else "ellipse", style="filled", fillcolor=to_hex(color),
             fontcolor=textcolor)
    # if object.escaped:
    #     dot.edge(str(name), str("escaped"))
    # if object.collided_with_sun:
    #     dot.edge(str(name), str("collided with sun"))

if environ.get("CI"):
    dot.save(fn.with_suffix(".gv"))
else:
    dot.render(fn.with_suffix(".gv"), view=True, format="svg")
