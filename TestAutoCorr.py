""" TestAutoCorr.py """

import matplotlib
matplotlib.use('TkAgg')  # Bypass use('macos') default incompatible with PyCharm
import matplotlib.pyplot as plt
from   PIL import Image, ImageDraw
from scipy import signal
import numpy as np
np.seterr(divide='ignore', invalid='ignore') # Ignore divide-by-zero errors


# Open figure window
winXSize = 15.  # approx inches
winYSize = 8.  # in width & height
winAspect = winXSize / winYSize
fig = plt.figure(figsize=(winXSize, winYSize))  # Open the Fig
fig.canvas.set_window_title('Oculus')

plt.ion()

# Get 'axes' (rectangular regions) for each image in the figure
def get_axes(fig, rect, title=None):
    ax = fig.add_axes(rect)
    ax.axes.set_xticks([])
    ax.axes.set_yticks([])
    if title:
        ax.set_title(title)
    return ax


# Axes for Images
axPSF    = get_axes(fig, [.05, .6, .36 / winAspect, .36], 'PSF')
axFour   = get_axes(fig, [.65, .6, .36 / winAspect, .36], 'Fourier')
axBefore = get_axes(fig, [.05, .2, .36 / winAspect, .36], 'PSF')
axAfter  = get_axes(fig, [.35, .6, .36 / winAspect, .36], 'After')
axDiag   = get_axes(fig, [.65, .2, .36 / winAspect, .36], 'Diagnostic')

# Generate PSF image
xSize, ySize = 256, 256
hafX, hafY   = xSize/2, ySize/2
imgPSF = np.zeros([2 * hafY, 2 * hafX])
pilPSF = Image.fromarray(imgPSF, 'L')

# Draw circular ellipse
draw = ImageDraw.Draw(pilPSF)
draw.ellipse(((-16 + hafX, -16 + hafY),
              ( 16 + hafX,  16 + hafY)),
              fill=255, outline=None)

# Normalize to range [0-1]
imgPSF = np.asarray(pilPSF) / 255.

# Display PSF image
plt.sca(axPSF)
psfPlot = plt.imshow(imgPSF, cmap='gray')

# Autocorrelation of the PSF
floatPSF = np.array(imgPSF, dtype=float)
autocorr = signal.correlate2d(floatPSF, floatPSF) # works - but takes 30 seconds!!!

# Copy to diagnostic image
imgDiag = np.array(autocorr, dtype=int)

# Display After image
plt.sca(axAfter)
diagPlot = plt.imshow(imgDiag, cmap='gray')

plt.show(block=True)


