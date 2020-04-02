from graphviz import Digraph

from extradata import ExtraData
from utils import filename_from_argv

ed = ExtraData.load(filename_from_argv().with_suffix(".extra.json"))

dot = Digraph(comment='Collisions')
for merged, originals in ed.tree._tree.items():
    for parent in originals["parents"]:
        dot.edge(str(parent), str(merged))

dot.render('graph.gv', view=True)
