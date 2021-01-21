from graphviz import Digraph

from extradata import ExtraData, CollisionMeta
from utils import filename_from_argv

fn = filename_from_argv()
ed = ExtraData.load(fn.with_suffix(".extra.json"))

dot = Digraph(comment='Collisions')

interacting_objects = set()
for merged, originals in ed.tree._tree.items():
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

for name in interacting_objects:
    object = ed.pdata[name]
    if object.type == "sun":
        name = "Sun"
    if object.type == "gas giant":
        name = "gas giant"
    dot.node(name=str(name), label=f"{name} ({object.water_mass_fraction:.2e})",
             shape='box' if object.type == "planetesimal" else "ellipse")

# dot.engine = 'neato'
dot.render(fn.with_suffix(".gv"), view=True, format="svg")
