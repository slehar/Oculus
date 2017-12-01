# -*- coding: utf-8 -*-
"""
TestReturnImg.py

Created on Thu Nov 30 10:44:58 2017

@author: slehar
"""

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from   PIL import Image
import numpy as np

plt.close('all')
fig = plt.figure(figsize=(6,5))
fig.canvas.set_window_title('TestReturnImg')
ax = fig.add_axes([.1,.1,.8,.8])
ax.axes.set_xticks([])
ax.axes.set_yticks([])

circ = mpatches.Circle((.5,.5), .25, color='r')
ax.add_patch(circ)

fig.canvas.draw()

def returnImage():
    return Image.frombytes('RGB', fig.canvas.get_width_height(),
                           fig.canvas.tostring_rgb())
                        
img = returnImage()
        
plt.show()


