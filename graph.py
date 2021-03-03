from os import environ

from graphviz import Digraph

from extradata import ExtraData, CollisionMeta
from utils import filename_from_argv

fn = filename_from_argv()
ed = ExtraData.load(fn)

dot = Digraph(comment='Collisions')

dot.engine = "neato"

if dot.engine == "neato":
    dot.attr("graph", overlap="false")
interacting_objects = set()
for merged, originals in ed.tree.get_tree().items():
    first_parent = True
    for parent in originals["parents"]:
        meta: CollisionMeta = originals["meta"]
        water_ret = meta.water_retention
        mass_ret = meta.mass_retention
        if first_parent:

            label = f"{water_ret:.2f}/{mass_ret:.2f}"
            first_parent = False
        else:
            label = None
        dot.edge(str(parent), str(merged), xlabel=label)
        interacting_objects.add(parent)
        interacting_objects.add(int(merged))

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
    dot.node(name=str(name), label=f"{displayname} ({object.water_mass_fraction:.2e}; {mass:.2e})",
             shape='box' if object.type == "planetesimal" else "ellipse")
    if object.escaped:
        dot.edge(str(name), str("escaped"))
    if object.collided_with_sun:
        dot.edge(str(name), str("collided with sun"))

if environ.get("CI"):
    dot.save(fn.with_suffix(".gv"))
else:
    dot.render(fn.with_suffix(".gv"), view=True, format="svg")
