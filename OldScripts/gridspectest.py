#!/usr/bin/env python
# from https://stackoverflow.com/questions/24738578/control-wspace-for-matplotlib-subplots

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np

fig = plt.figure()
# create a 1-row 3-column container as the left container
gs_left = gridspec.GridSpec(1, 3)

# create a 1-row 1-column grid as the right container
gs_right = gridspec.GridSpec(1, 1)

# add plots to the nested structure
ax1 = fig.add_subplot(gs_left[0,0])
ax2 = fig.add_subplot(gs_left[0,1])
ax3 = fig.add_subplot(gs_left[0,2])

# create a 
ax4 = fig.add_subplot(gs_right[0,0])

# now the plots are on top of each other, we'll have to adjust their edges so that they won't overlap
gs_left.update(right=0.65)
gs_right.update(left=0.7)

# also, we want to get rid of the horizontal spacing in the left gridspec
gs_left.update(wspace=0)

plt.show()
