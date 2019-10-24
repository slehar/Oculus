# TestFourier.py
'''
    Test the forward and inverse transform process stream with graphics
'''

import matplotlib
matplotlib.use('TkAgg')

import matplotlib.pyplot as plt
from   PIL import Image, ImageDraw
import numpy as np
import utils

radius = 5.
snr = 100.

# Open figure window
winXSize = 15  # approx inches
winYSize = 8  # in width & height
winAspect = winXSize / winYSize
fig = plt.figure(figsize=(winXSize, winYSize))  # Open the Fig
fig.canvas.set_window_title('TestFourier')
# fig.canvas.mpl_connect('key_press_event', press)

plt.ion()
plt.show()

# Axes for Images
axBefore = utils.get_axes(fig, [.05, .6, .25 / winAspect, .35], 'Before')
axAfter  = utils.get_axes(fig, [.35, .6, .25 / winAspect, .35], 'After')
axFour   = utils.get_axes(fig, [.65, .6, .25 / winAspect, .35], 'Fourier')
axPSF    = utils.get_axes(fig, [.05, .2, .25 / winAspect, .35], 'PSF')
axDiag   = utils.get_axes(fig, [.65, .2, .25 / winAspect, .35], 'Diagnostic')

# Input image for Before
imgPil = Image.open('blr.png').convert('LA')  # LA = 'luminance + alpha
imgNp = np.array(imgPil.convert('L')) / 255.  # convert to L, luminance only
ySize, xSize = imgNp.shape
hafY, hafX = int(ySize / 2), int(xSize / 2)
ySize, xSize = imgNp.shape

# Display Before image
plt.sca(axBefore)
imgPlot = plt.imshow(imgNp, cmap='gray')

# First Pass Hanning
han = np.outer(np.hanning(ySize), np.hanning(xSize))
imgHan = imgNp * han  # apply hanning window

# Display hanning image
plt.sca(axBefore)
imgplot = plt.imshow(imgHan, cmap='gray')
x, y = hafX, hafY

# Generate PSF image
K = np.zeros(imgNp.shape)  # array for inverse filter
yy, xx = np.mgrid[-hafY:hafY, -hafX:hafX]
distImg = np.sqrt(xx ** 2 + yy ** 2)
imgPSF = (distImg < radius)  # This is a Disc PSF

# Display PSF image
plt.sca(axPSF)
psfPlot = plt.imshow(imgPSF, cmap='gray')

# Fourier Transform
fourImg = np.fft.fft2(imgHan)  # set dc term to 1 to control contrast
fourImg[0, 0] = 1.0 + 0j
fourShft = np.fft.fftshift(fourImg)
fourLog = np.log(np.abs(fourShft))
fourLog = fourLog / complex(fourLog.max())

# Display Fourier image
plt.sca(axFour)
fourPlot = plt.imshow(fourLog.real, cmap='gray')
# plt.pause(.001)

# take transform of psf
fourPSF = np.fft.fft2(imgPSF)
fourShftPSF = np.fft.fftshift(fourPSF)
psfLog = np.log(np.maximum(np.abs(fourShftPSF), 1.))
psfLog = psfLog / complex(psfLog.max())

# Display PSF transform image
fourPlot.set_data(psfLog.real)

# Create the Linear MAP filter, K(u,v)
isnr = 1. / snr
conjfourPSF = np.conj(fourPSF)
K = (conjfourPSF + isnr) / ((conjfourPSF * fourPSF) + isnr)
KLog = np.log(np.maximum(np.abs(K), 1.))
KLog = KLog / complex(KLog.max())  # normalizing 0 to 1


# Display KLog image
plt.sca(axDiag)
diagPlot = plt.imshow(KLog.real, cmap='gray')

# do the inverse filtering
fourResult = fourShft * K  # convolution in the fourier domain

# Inverse Fourier Transform
fourIshft = np.fft.ifftshift(fourResult)
fourIshft[0, 0] = 0.5 + 0.0j  # set d.c. term for display
fourInv = np.fft.ifft2(fourIshft)

#   make sure fourReal scales 0.to 1.0 for display
fourInv = np.fft.ifftshift(fourInv)

plt.sca(axAfter)
invPlot = plt.imshow(fourInv.real, cmap='gray')

plt.show(block=True)

