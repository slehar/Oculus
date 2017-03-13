# -*- coding: utf-8 -*-
"""
FourInteractive.py

Created on Wed Sep 28 16:36:45 2016

@author: slehar
"""

import matplotlib.pyplot as plt
#import matplotlib.patches as patches #for drawing lines, rectangles
from   matplotlib.widgets import Slider
from   matplotlib.widgets import RadioButtons
from   PIL import Image, ImageDraw 
import Tkinter, tkFileDialog
import numpy as np
#import numpy.ma as ma
import sys

# Global Variables
radius = 5.
rad2=0
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
imgFile = tkFileDialog.askopenfilename(initialfile = 'blr.png')
#imgFile = 'blr.png'
# Open figure window
winXSize = 18
winYSize = 14
winAspect = winXSize/winYSize
plt.close('all')
fig = plt.figure(figsize=(winXSize, winYSize))
fig.canvas.set_window_title('First Oculus Attempt')


# Keypress 'q' to quit callback function
def press(event):
    global ptList, data
    sys.stdout.flush()
    if event.key == 'q':
        plt.close()
# Connect keypress event to callback function
fig.canvas.mpl_connect('key_press_event', press)


# 
# Disc / Line Checkbox
rax = plt.axes([0.41, 0.1, 0.06/winAspect, 0.1])
radio = RadioButtons(rax, ['Disc', 'Line'])


# Axes for Before Image
axBefore = fig.add_axes([.05, .6, .35/winAspect, .35])
axBefore.axes.set_xticks([])
axBefore.axes.set_yticks([])
axBefore.set_title('before')

# Read image and display
imgPil = Image.open(imgFile).convert('LA') #LA = 'luminance + alpha
imgNp = np.array(imgPil.convert('L'))/255. # convert to L, luminance only
ySize, xSize = imgNp.shape

print 'original image size', xSize, ySize
#here we might like to resize the array if it is not even, by cutting out
# a line on any odd-dimension.   
if ySize % 2 == 0:
    print 'y dimension is even'
else:
    print 'y dimension is odd, so drop a row'
    ySize = ySize-1
    imgNp = imgNp[0:ySize]
    
if xSize % 2 == 0:
    print 'x dimension is even'
else:
    print 'x dimension is odd, so drop a column'
    xSize = xSize-1
    imgNp=imgNp[:, :xSize]
    

print 'input image', 'xSize=', xSize, 'ySize=', ySize, 'pixels'

hafY, hafX = int(ySize/2), int(xSize/2)
#print 'input image', 'xSize=', xSize, 'ySize=', ySize, 'pixels'
plt.sca(axBefore)
beforePlot = plt.imshow(imgNp, cmap='gray',vmin=0.,vmax=1.)

# Axes for Diagnostic window
axDiag = fig.add_axes([.65, .2, .35/winAspect, .35])
axDiag.axes.set_xticks([])
axDiag.axes.set_yticks([])
axDiag.set_title('Diagnostic')

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

# Axes for Fourier Image
axFour = fig.add_axes([.65, .6, .35/winAspect, .35])
axFour.axes.set_xticks([])
axFour.axes.set_yticks([])
axFour.set_title('Fourier')

# Axes for PSF Image
axPSF = fig.add_axes([.05, .2, .35/winAspect, .35])
axPSF.axes.set_xticks([])
axPSF.axes.set_yticks([])
axPSF.set_title('PSF')

yy, xx = np.mgrid[-hafY:hafY, -hafX:hafX]

# Fourier Transform
fourImg  = np.fft.fft2(imgHan) #set dc term to 1 to control contrast
fourImg[0,0] = 1.0+0j
print 'DC term is', fourImg[1:2,1:2]
fourShft = np.fft.fftshift(fourImg)
fourLog  = np.log(np.abs(fourShft))
fourLog = fourLog/complex(fourLog.max())

plt.sca(axFour)
fourPlot = plt.imshow(fourLog.real, cmap='gray')
plt.pause(.001)

#### Fourier Filtering ####
distImg = np.sqrt(xx**2 + yy**2)

# Generate PSF Image
#psfMode = 'Line'
if psfMode == 'Disc':
    imgPSF = (distImg < radius)# This is a Disc PSF
elif psfMode == 'line':
    lineLength = radius*3.
    imgPSF = np.zeros([2*hafY, 2*hafX])
    pilPSF = Image.fromarray(imgPSF, 'L')
    draw = ImageDraw.Draw(pilPSF)
    draw.Line(((lineLength * -np.cos(lineOri)+ hafX, lineLength * -np.sin(lineOri)+ hafY,
                lineLength *  np.cos(lineOri)+ hafX, lineLength *  np.sin(lineOri)+ hafY)), 
               fill=255, width=lineWidth)
    imgPSF = np.asarray(pilPSF)/255.
    
imgPSF = imgPSF.astype(float) 
plt.sca(axPSF) # Set imgPSF the "current axes"
psfPlot = plt.imshow(imgPSF, cmap='gray')
    
filtImg=fourShft * imgPSF
filtLog = np.log(np.maximum(np.abs(filtImg),1.))


# Axes for Inverse After Image
axAfter = fig.add_axes([.35, .6, .35/winAspect, .35])
axAfter.axes.set_xticks([])
axAfter.axes.set_yticks([])
axAfter.set_title('After')

# Inverse Fourier Transform
fourIshft = np.fft.ifftshift(filtImg)
fourInv   = np.fft.ifft2(fourIshft)
fourReal  = np.real(fourInv)
plt.sca(axAfter)
invPlot = plt.imshow(fourReal, cmap='gray', vmin=0, vmax=1)

# Filter radius sliders
axSlider1 = fig.add_axes([0.41, 0.5, 0.234, 0.04])
axSlider1.set_xticks([])
axSlider1.set_yticks([])

axSlider2 = fig.add_axes([0.41, 0.4509, 0.234, 0.04])
axSlider2.set_xticks([])
axSlider2.set_yticks([])

hafRadiusMax = min(ySize,xSize) # Setting upper limit to radius focus blur
slider1 = Slider(axSlider1, 'radius', 0.5, hafRadiusMax/10, valinit=5.)
rad1 = np.log(slider1.val)      #make log control
slider1.valtext.set_text(rad1)

slider2 = Slider(axSlider2, 'angle', -np.pi, np.pi, valinit=0.)
lineOri = slider2.val

# Filter angular sliders
axSlider3 = fig.add_axes([0.41, 0.40, 0.234, 0.04])
axSlider3.set_xticks([])
axSlider3.set_yticks([])
slider3 = Slider(axSlider3, 'snr',  1., 1000., valinit=100.)
snr = slider3.val

axSlider4 = fig.add_axes([0.41, 0.35, 0.234, 0.04])
axSlider4.set_xticks([])
axSlider4.set_yticks([])
slider4 = Slider(axSlider4, 'linewidth', 1, 50, valinit=5, valfmt='%d')
lineWidth = slider4.val

plt.sca(axDiag)
diagPlot = plt.imshow(imgHan, cmap='gray',vmin=0.,vmax=1.)

#%%
# This is loop where all the action happens
def update():
    global filtImg, imgPSF, K, KLog, psfLog, resultReal, isnr, snr, fourImg, totalpixels, imgNp, fourReal
    
    # PSF in axPSF
    if psfMode == 'Disc':
        imgPSF = (distImg < radius)# This is a Disc PSF
    elif psfMode == 'Line':
        lineLength = radius*3.
        imgPSF = np.zeros([2*hafY, 2*hafX])
        pilPSF = Image.fromarray(imgPSF, 'L')
        draw = ImageDraw.Draw(pilPSF)
        draw.line(((lineLength * -np.cos(lineOri)+ hafX, lineLength * -np.sin(lineOri)+ hafY,
                    lineLength *  np.cos(lineOri)+ hafX, lineLength *  np.sin(lineOri)+ hafY)), 
                   fill=255, width=lineWidth)
        imgPSF = np.asarray(pilPSF)/255.
    
    realimgPSF = imgPSF.astype(float) 
    psfPlot.set_data(realimgPSF)
          
    #take transform of psf and displaying it 
    fourPSF = np.fft.fft2(imgPSF)
    fourShftPSF = np.fft.fftshift(fourPSF)
    psfLog = np.log(np.maximum(np.abs(fourShftPSF),1.))
    psfLog = psfLog/complex(psfLog.max())
    fourPlot.set_data(psfLog.real)
          
    # Create the Linear MAP filter, K(u,v) 
    isnr=1./snr
    #-------
    conjfourPSF = np.conj(fourPSF)
    K=(conjfourPSF + isnr)/((conjfourPSF*fourPSF)+isnr)
#    print K[hafY:hafY+1,hafX:hafX+1]    #is this the d.c. term?
    KLog = np.log(np.maximum(np.abs(K),1.))
    KLog = KLog/complex(KLog.max()) #normalizing 0 to 1
    diagPlot.set_data(KLog.real) #Plotting diag data
    norm = np.sum(K)    #for normalizing K

    #do the inverse filtering
    fourResult=fourShft*K    #convolution in the fourier domain
    
    # Inverse Fourier Transform
    fourIshft = np.fft.ifftshift(fourResult)
    fourIshft[0,0] = 1.0+0.0j #set d.c. term for display

    fourInv   = np.fft.ifft2(fourIshft)
   
#   make sure fourReal scales 0.to 1.0 for display
    fourInv=np.fft.ifftshift(fourInv)
    fourReal  = np.real(fourInv)

#    lmax=fourReal.max()
#    lmin=fourReal.min()

    plt.sca(axAfter)
    invPlot = plt.imshow(fourReal, cmap='gray')

    plt.pause(.001)

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

plt.show()

# Pop fig window to top]]
figmgr=plt.get_current_fig_manager()
figmgr.canvas.manager.window.raise_()
geom=figmgr.window.geometry()
(xLoc,yLoc,dxWidth,dyHeight)=geom.getRect()
figmgr.window.setGeometry(10,10,dxWidth,dyHeight)
 

