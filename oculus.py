# -*- coding: utf-8 -*-
"""
oculus.py

Created on Wed Sep 28 16:36:45 2016

@author: slehar
"""

import matplotlib.pyplot as plt
from   matplotlib.widgets import RadioButtons
from   PIL import Image, ImageDraw 
import Tkinter, tkFileDialog
import numpy as np
import sys
import blob
import utils

# Global Variables
radius = 5.
# rad2=0
snr = 100.
isnr=1./snr
slidersLocked = False
psfMode = 'Disc'
#psfMode = 'Line'
lineOri = 0.
lineLength = 100.
lineWidth = 5


# Get image using finder dialog
root = Tkinter.Tk()
root.withdraw() # Hide the root window
imgFile = tkFileDialog.askopenfilename(title='Select Image',
                                       initialdir = ".",
                                       initialfile = 'blr.png')
# Open figure window
winXSize = 15    # approx inches
winYSize = 8     # in width & height
winAspect = winXSize/winYSize
fig = plt.figure(figsize=(winXSize, winYSize))
fig.canvas.set_window_title('Oculus')

# Keypress 'q' to quit callback function
def press(event):
    global ptList, data
    sys.stdout.flush()
    print 'Keypress detected'
    if event.key == 'q':
        print "key is 'q'"
        plt.close()
# Connect keypress event to callback function
fig.canvas.mpl_connect('key_press_event', press)


# 
# Disc / Line / blob Checkbox
rax = plt.axes([0.41, 0.1, 0.06/winAspect, 0.1])
radio = RadioButtons(rax, ['Disc', 'Line', 'Blob'])


# Axes for Images
axBefore = utils.get_axes(fig, [.05, .6, .35/winAspect, .35], 'Before')
axFour   = utils.get_axes(fig, [.65, .6, .35/winAspect, .35], 'Fourier')
axPSF    = utils.get_axes(fig, [.05, .2, .35/winAspect, .35], 'PSF')
axAfter  = utils.get_axes(fig, [.35, .6, .35/winAspect, .35], 'After')
axDiag   = utils.get_axes(fig, [.65, .2, .35/winAspect, .35], 'Diagnostic')

# Turn on interactive mode
plt.ion()

# Read image and display
imgPil = Image.open(imgFile).convert('LA') #LA = 'luminance + alpha
imgNp = np.array(imgPil.convert('L'))/255. # convert to L, luminance only
ySize, xSize = imgNp.shape
print 'original image size', xSize, ySize

imgNp = utils.resize_even_dim(imgNp)
ySize, xSize = imgNp.shape

print 'resized image', 'xSize=', xSize, 'ySize=', ySize, 'pixels'



# Slider 1
hafRadiusMax = min(ySize,xSize) # Setting upper limit to radius focus blur
slider1 = utils.get_slider(fig, [0.41, 0.5, 0.234, 0.04],   'radius', 0.5, hafRadiusMax/10, valinit=5.)
rad1 = np.log(slider1.val)      #make log control
slider1.valtext.set_text(rad1)

# Slider 2
slider2 = utils.get_slider(fig, [0.41, 0.4509, 0.234, 0.04], 'angle', -np.pi, np.pi, valinit=0.)
lineOri = slider2.val

# Slider 3
slider3 = utils.get_slider(fig, [0.41, 0.40, 0.234, 0.04], 'SNR', 1., 1000., valinit=100. )
snr = slider3.val

# Slider 4
slider4 = utils.get_slider(fig, [0.41, 0.35, 0.234, 0.04],'linewidth', 1, 50, valinit=5 )
lineWidth = slider4.val

# Slider 5
slider5 = utils.get_slider(fig, [0.41, 0.3, 0.234, 0.04], 'skew', 1, 50, valinit=0)
skew = slider5.val

#########################[ First Pass Initialization ]########################################

hafY, hafX = int(ySize/2), int(xSize/2)
#print 'input image', 'xSize=', xSize, 'ySize=', ySize, 'pixels'
plt.sca(axBefore)
beforePlot = plt.imshow(imgNp, cmap='gray',vmin=0.,vmax=1.)

# First Pass Hanning
han = np.outer(np.hanning(ySize),np.hanning(xSize))
imgHan = imgNp*han #apply hanning window
totalpixels = xSize*ySize


#print 'total pixels = ', totalpixels
#radius = min(ySize,xSize)/4
plt.sca(axDiag)
imgplot = plt.imshow(imgNp, cmap='gray')
x,y=hafX,hafY
#sigma=xSize/100.

K=np.zeros(imgNp.shape)#array for inverse filter

#psf image

yy, xx = np.mgrid[-hafY:hafY, -hafX:hafX]
distImg = np.sqrt(xx**2 + yy**2)
imgPSF = (distImg < radius)  # This is a Disc PSF

psfPlot = plt.imshow(imgPSF, cmap='gray')
plt.sca(axFour)

# Fourier Transform
fourImg  = np.fft.fft2(imgHan) #set dc term to 1 to control contrast
fourImg[0,0] = 1.0+0j
print 'DC term is', fourImg[1:2,1:2]
fourShft = np.fft.fftshift(fourImg)
fourLog  = np.log(np.abs(fourShft))
fourLog = fourLog/complex(fourLog.max())

plt.sca(axFour)
fourPlot = plt.imshow(fourLog.real, cmap='gray')
# plt.pause(.001)

plt.sca(axDiag)
diagPlot = plt.imshow(imgHan, cmap='gray',vmin=0.,vmax=1.)

imgPSF = imgPSF.astype(float)
plt.sca(axPSF) # Set imgPSF the "current axes"
psfPlot = plt.imshow(imgPSF, cmap='gray')


#########################[ def update(): ]########################################
# This is loop where all the action happens
def update():
    global filtImg, imgPSF, K, KLog, psfLog, resultReal, isnr, snr, fourImg, totalpixels, imgNp, fourReal

    # PSF in axPSF
    if psfMode == 'Disc':
        imgPSF = (distImg < radius)  # This is a Disc PSF
    elif psfMode == 'Line':
        lineLength = radius * 3.
        imgPSF = np.zeros([2 * hafY, 2 * hafX])
        pilPSF = Image.fromarray(imgPSF, 'L')
        draw = ImageDraw.Draw(pilPSF)
        draw.line(((lineLength * -np.cos(lineOri) + hafX, lineLength * -np.sin(lineOri) + hafY,
                    lineLength * np.cos(lineOri) + hafX, lineLength * np.sin(lineOri) + hafY)),
                  fill=255, width=lineWidth)
        imgPSF = np.asarray(pilPSF) / 255.
    elif psfMode == 'Blob':
        imgPSF = blob.returnBlobImage()

    realimgPSF = imgPSF.astype(float)
    psfPlot.set_data(realimgPSF)

    # take transform of psf and displaying it
    fourPSF = np.fft.fft2(imgPSF)
    fourShftPSF = np.fft.fftshift(fourPSF)
    psfLog = np.log(np.maximum(np.abs(fourShftPSF), 1.))
    psfLog = psfLog / complex(psfLog.max())
    fourPlot.set_data(psfLog.real)

    # Create the Linear MAP filter, K(u,v)
    isnr = 1. / snr
    # -------
    conjfourPSF = np.conj(fourPSF)
    K = (conjfourPSF + isnr) / ((conjfourPSF * fourPSF) + isnr)
    #    print K[hafY:hafY+1,hafX:hafX+1]    #is this the d.c. term?
    KLog = np.log(np.maximum(np.abs(K), 1.))
    KLog = KLog / complex(KLog.max())  # normalizing 0 to 1
    diagPlot.set_data(KLog.real)  # Plotting diag data
    norm = np.sum(K)  # for normalizing K

    # do the inverse filtering
    fourResult = fourShft * K  # convolution in the fourier domain

    # Inverse Fourier Transform
    fourIshft = np.fft.ifftshift(fourResult)
    fourIshft[0, 0] = 0.5 + 0.0j  # set d.c. term for display

    fourInv = np.fft.ifft2(fourIshft)

    #   make sure fourReal scales 0.to 1.0 for display
    fourInv = np.fft.ifftshift(fourInv)
    fourReal = np.real(fourInv)

    #    lmax=fourReal.max()
    #    lmin=fourReal.min()

    plt.sca(axAfter)
    invPlot = plt.imshow(fourReal, cmap='gray')

    plt.pause(.001)


####################[ end def update(): ] ###############@#@



# radio button callback function to switch PSF mode
def modefunc(label):
    global psfMode
    
    if   label == 'Disc':
        psfMode = 'Disc'
        imgPSF = (distImg < radius)# This is a Disc PSF
        slider1.label = 'radius'
    if   label == 'Line':
        psfMode = 'Line'
        lineLength = radius * 3
        imgPSF = np.zeros([2*hafY, 2*hafX])
        pilPSF = Image.fromarray(imgPSF, 'L')
        slider1.label = 'length'
        print slider1.label
        draw = ImageDraw.Draw(pilPSF)
        draw.line(((lineLength * -np.cos(lineOri)+ hafX, lineLength * -np.sin(lineOri)+ hafY,
                    lineLength *  np.cos(lineOri)+ hafX, lineLength *  np.sin(lineOri)+ hafY)), 
                   fill=255, width=lineWidth)
        imgPSF = np.asarray(pilPSF)/255.
        
    if label == 'Blob':
        psfMode = 'Blob'
        blobFig, blobAx = blob.openBlobWindow()
        fig.canvas.draw()
        imgPSF = blob.returnBlobImage()
        
        
    imgPSF = imgPSF.astype(float)
    plt.show()
    plt.pause(.001)
    update()
    plt.draw()    
radio.on_clicked(modefunc)

#%%
def update1(val):
    global radius
    radius = slider1.val
    update()

def update2(val):
    global lineOri
    lineOri = slider2.val
    update()

def update3(val):
    global snr
    snr = slider3.val
    update()

def update4(val):
    global lineWidth
    lineWidth = int(slider4.val)
    update()

#    fig.canvas.draw()
slider1.on_changed(update1)
slider2.on_changed(update2)
slider3.on_changed(update3)
slider4.on_changed(update4)

update()

# Show [block = leave it open, don't close]
plt.show(block=True)



