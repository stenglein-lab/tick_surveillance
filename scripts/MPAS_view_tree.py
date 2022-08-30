#!/usr/bin/env python

import toytree
import toyplot
import toyplot.pdf
import numpy as np
import os
import sys

# setup command one line arguments:

fasttree: sys.argv[1]

# create tree pdf with toytree python tool. Documentation found here: https://toytree.readthedocs.io/en/latest/


tree = toytree.tree(sys.argv[1])
canvas, axes, makr = tree.draw(
        tip_labels_style={"font-size":"11px", "-toyplot-anchor-shift": "20px",},
        tip_labels_align = False,
        height=1000,
        width = 1000,
        node_hover=True, 
        node_sizes=5, 
        node_colors='teal', 
        tree_style='n',
        scalebar=True);
toyplot.pdf.render(canvas, f'{sys.argv[1]}.pdf')
