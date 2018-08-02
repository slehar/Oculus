# -*- coding: utf-8 -*-
"""
Created on Wed Feb 14 09:39:08 2018

@author: slehar
"""
print '*** In callblob ***'

import matplotlib.pyplot as plt
import numpy as np
#import matplotlib.image as mpimg

import blob

xSize, ySize = 512, 512
hafX, hafY = int(xSize/2), int(ySize/2)


fig = plt.figure(figsize=(6, 6))
fig.canvas.set_window_title('CallBlob')
ax = fig.add_axes([.1, .1, .8, .8])

imgPSF = np.zeros([2*hafY, 2*hafX])



blob.openBlobWindow(imgPSF)

imgplot = plt.imshow(imgPSF)



