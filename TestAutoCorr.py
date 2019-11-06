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

def xcorr(x):
    """FFT based autocorrelation function, which is faster than numpy.correlate"""
    # x is supposed to be an array of sequences, of shape (totalelements, length)
    # https://stackoverflow.com/questions/4503325/autocorrelation-of-a-multidimensional-array-in-numpy
    fftx = np.fft.fft2(x)
    ret = np.fft.ifft(fftx * np.conjugate(fftx), axis=1)
    ret = np.fft.fftshift(ret, axes=1)
    return ret

from itertools import product
from numpy import empty, roll

def autocorrelate(x):
    """
    Compute the multidimensional autocorrelation of an nd array.
    input: an nd array of floats
    output: an nd array of autocorrelations
    """

    # used for transposes
    t = roll(range(x.ndim), 1)

    # pairs of indexes
    # the first is for the autocorrelation array
    # the second is the shift
    ii = [list(enumerate(range(1, s - 1))) for s in x.shape]

    # initialize the resulting autocorrelation array
    acor = empty(shape=[len(s0) for s0 in ii])

    # iterate over all combinations of directional shifts
    for i in product(*ii):
        # extract the indexes for
        # the autocorrelation array
        # and original array respectively
        i1, i2 = np.asarray(i).T

        x1 = x.copy()
        x2 = x.copy()

        for i0 in i2:
            # clip the unshifted array at the end
            x1 = x1[:-i0]
            # and the shifted array at the beginning
            x2 = x2[i0:]

            # prepare to do the same for
            # the next axis
            x1 = x1.transpose(t)
            x2 = x2.transpose(t)

        # normalize shifted and unshifted arrays
        x1 -= x1.mean()
        x1 /= x1.std()
        x2 -= x2.mean()
        x2 /= x2.std()

        # compute the autocorrelation directly
        # from the definition
        acor[tuple(i1)] = (x1 * x2).mean()

    return acor

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

choice = 'signal.correlate2d'
# choice = 'xcorr'
# choice = 'autocorrelate'


if choice == 'signal.correlate2d':
    # Using signal.correlate2d
    autocorr = signal.correlate2d(floatPSF, floatPSF) # works - but takes 30 seconds!!!

elif choice == 'xcorr':
    # Using xcorr
    autocorr = xcorr(floatPSF)

elif choice == 'autocorrelate':
    # Using autocorrelate
    autocorr = autocorrelate(floatPSF) / 255. # Takes 1:15 min:sec and is wrong!


# Copy to diagnostic image
imgDiag = np.array(autocorr, dtype=int)

# Display After image
plt.sca(axAfter)
diagPlot = plt.imshow(imgDiag, cmap='gray')

plt.show(block=True)


