from graphviz import Digraph

from extradata import ExtraData
from utils import filename_from_argv

fn = filename_from_argv()
ed = ExtraData.load(fn.with_suffix(".extra.json"))

dot = Digraph(comment='Collisions')
for merged, originals in ed.tree._tree.items():
    first_parent = True
    for parent in originals["parents"]:
        meta = originals["meta"]
        water_ret = meta["water_retention"]
        mass_ret = meta["mass_retention"]
        if first_parent:

            label = f"{water_ret:.2f}/{mass_ret:.2f}"
            first_parent = False
        else:
            label = None
        dot.edge(str(parent), str(merged), xlabel=label)

# dot.engine = 'neato'
dot.render(fn.with_suffix(".gv"), view=False, format="svg")
