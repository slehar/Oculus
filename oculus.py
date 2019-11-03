# -*- coding: utf-8 -*-
"""
oculus.py

Created on Wed Sep 28 16:36:45 2016

@author: slehar
"""

import matplotlib
matplotlib.use('TkAgg')  # Bypass use('macos') default incompatible with PyCharm
import matplotlib.pyplot as plt
from   matplotlib.widgets import RadioButtons
from   PIL import Image, ImageDraw
import Tkinter, tkFileDialog
import numpy as np
np.seterr(divide='ignore', invalid='ignore') # Ignore divide-by-zero errors
import utils
import blob

class pltFig():
    ''' Plot Figure class '''

    radius = 5.
    snr = 100.
    isnr=1./snr
    slidersLocked = False
    psfMode = 'Disc'
    lineOri = 0.
    lineLength = 100.
    lineWidth = 5.
    lineSkew = 0.

    def __init__(self):

        # Open figure window
        winXSize = 15.  # approx inches
        winYSize = 8.   # in width & height
        self.winAspect = winXSize / winYSize
        self.fig = plt.figure(figsize=(winXSize, winYSize))  # Open the Fig
        self.fig.canvas.set_window_title('Oculus')

        # Axes for Images
        self.axBefore = utils.get_axes(self.fig, [.05, .6, .36 / self.winAspect, .36], 'Before')
        self.axFour   = utils.get_axes(self.fig, [.65, .6, .36 / self.winAspect, .36], 'Fourier')
        self.axPSF    = utils.get_axes(self.fig, [.05, .2, .36 / self.winAspect, .36], 'PSF')
        self.axAfter  = utils.get_axes(self.fig, [.35, .6, .36 / self.winAspect, .36], 'After')
        self.axDiag   = utils.get_axes(self.fig, [.65, .2, .36 / self.winAspect, .36], 'Diagnostic')

        # Disc / Line / Parallelogram Checkbox
        self.rax = plt.axes([0.41, 0.05, 0.2 / self.winAspect, 0.2])
        self.radio = RadioButtons(self.rax, ['Disc', 'Line', 'Parallelogram'])

        # Read Before image
        imgFile = self.get_image()
        imgPil = Image.open(imgFile).convert('LA')  # LA = 'luminance + alpha
        self.imgNp = np.array(imgPil.convert('L')) / 255.  # convert to L, luminance only
        self.ySize, self.xSize = self.imgNp.shape
        self.hafY, self.hafX = int(self.ySize / 2), int(self.xSize / 2)

        # Convert to numpy array
        self.imgNp = utils.resize_even_dim(self.imgNp)
        self.ySize, self.xSize = self.imgNp.shape

        # Display before image
        plt.sca(self.axBefore)
        self.beforePlot = plt.imshow(self.imgNp, cmap='gray')
        plt.pause(.15)

        # First Pass Hanning
        han = np.outer(np.hanning(self.ySize), np.hanning(self.xSize))
        self.imgHan = self.imgNp * han  # apply hanning window

        # Display hanning image
        self.beforePlot.set_data(self.imgHan)

        # Slider 1
        hafRadiusMax = min(self.ySize, self.xSize)  # Setting upper limit to radius focus blur
        self.slider1 = utils.get_slider(self.fig, [0.35, 0.5, 0.234, 0.04], 'radius', 0.5, hafRadiusMax / 10, valinit=5.)
        rad1 = np.log(self.slider1.val)  # make log control
        self.slider1.valtext.set_text(rad1)
        self.slider1.on_changed(self.update1)

        # Set PSF image based on mode
        self.imgPSF = np.zeros([2 * self.hafY, 2 * self.hafX])
        self.pilPSF = Image.fromarray(self.imgPSF, 'L')
        # self.draw   = ImageDraw.Draw(self.pilPSF)

        # Generate PSF image
        self.modefunc('Disc')

        # Display PSF image
        plt.sca(self.axPSF)
        self.psfPlot = plt.imshow(self.imgPSF, cmap='gray')

        # Do Fourier filtering
        self.fourier_filter()

        # Display Fourier image
        plt.sca(self.axFour)
        self.fourPlot = plt.imshow(self.fourLog.real, cmap='gray')

        # Display Fourier inverse image
        plt.sca(self.axAfter)
        self.invPlot = plt.imshow(self.fourInv.real, cmap='gray')

        # Display KLog image
        plt.sca(self.axDiag)
        self.diagPlot = plt.imshow(self.KLog.real, cmap='gray')

        # Display After image
        plt.sca(self.axAfter)

        # Slider 2
        self.slider2 = utils.get_slider(self.fig, [0.35, 0.4509, 0.234, 0.04], 'angle', -np.pi, np.pi, valinit=0.)
        self.lineOri = self.slider2.val
        self.slider2.on_changed(self.update2)

        # Slider 3
        self.slider3 = utils.get_slider(self.fig, [0.35, 0.40, 0.234, 0.04], 'SNR', 1., 1000., valinit=100.)
        self.snr = self.slider3.val
        self.slider3.on_changed(self.update3)

        # Slider 4
        self.slider4 = utils.get_slider(self.fig, [0.35, 0.35, 0.234, 0.04], 'linewidth', 1, 50, valinit=5)
        self.lineWidth = self.slider4.val
        self.slider4.on_changed(self.update4)

        # Slider 5
        self.slider5 = utils.get_slider(self.fig, [0.35, 0.3, 0.234, 0.04], 'skew',  -np.pi, np.pi, valinit=0)
        self.lineSkew = self.slider5.val
        self.slider5.on_changed(self.update5)

        self.beforePlot.set_data(self.imgHan)

        self.diagPlot.set_data(self.imgHan)

        self.radio.on_clicked(self.modefunc)

        self.fig.canvas.mpl_connect('key_press_event', self.press)


    def update(self):
        ''' Update display when sliders change '''


        self.fourier_filter()


        # Update filtered Fourier image
        self.fourPlot.set_data(self.fourLog.real)

        # Update PSF image
        self.psfPlot.set_data(self.imgPSF)

        # Update inverse after image
        self.invPlot.set_data(self.fourInv.real)

        # Display diagnostic KLog image
        self.diagPlot.set_data(self.KLog.real)








    def modefunc(self, label):


        if label == 'Disc':
            self.psfMode = 'Disc'

            # Generate disc PSF
            self.K = np.zeros(self.imgNp.shape)  # array for inverse filter
            yy, xx = np.mgrid[-self.hafY:self.hafY, -self.hafX:self.hafX]
            self.distImg = np.sqrt(xx ** 2 + yy ** 2)

            self.imgPSF = (self.distImg < self.radius)  # This is a Disc PSF
            self.slider1.label = 'radius'

        if label == 'Line':
            self.psfMode = 'Line'
            self.lineLength = self.radius * 3.
            imgPSF = np.zeros([2 * self.hafY, 2 * self.hafX])
            pilPSF = Image.fromarray(imgPSF, 'L')
            self.slider1.label = 'length'
            draw = ImageDraw.Draw(pilPSF)
            draw.line(((
                self.lineLength * -np.cos(self.lineOri) + self.hafX,
                self.lineLength * -np.sin(self.lineOri) + self.hafY,
                self.lineLength * np.cos(self.lineOri) + self.hafX,
                self.lineLength * np.sin(self.lineOri) + self.hafY)),
                fill=255, width=int(self.lineWidth))
            self.imgPSF = np.asarray(pilPSF) / 255.

        if label == 'Parallelogram':
            self.psfMode = 'Parallelogram'
            self.slider1.label = 'length'

            self.lineLength = self.radius
            imgPSF = np.zeros([2 * self.hafY, 2 * self.hafX])
            pilPSF = Image.fromarray(imgPSF, 'L')
            draw = ImageDraw.Draw(pilPSF)
            draw.polygon((
                       (self.lineLength * -np.cos(self.lineOri) + self.lineWidth * np.sin(self.angleMod(self.lineSkew - self.lineOri)) + self.hafX,
                        self.lineLength * -np.sin(self.lineOri) + self.lineWidth * np.cos(self.angleMod(self.lineSkew - self.lineOri)) + self.hafY),
                       (self.lineLength * -np.cos(self.lineOri) - self.lineWidth * np.sin(self.angleMod(self.lineSkew - self.lineOri)) + self.hafX,
                        self.lineLength * -np.sin(self.lineOri) - self.lineWidth * np.cos(self.angleMod(self.lineSkew - self.lineOri)) + self.hafY),
                       (self.lineLength * np.cos(self.lineOri)  - self.lineWidth * np.sin(self.angleMod(self.lineSkew - self.lineOri)) + self.hafX,
                        self.lineLength * np.sin(self.lineOri)  - self.lineWidth * np.cos(self.angleMod(self.lineSkew - self.lineOri)) + self.hafY),
                       (self.lineLength * np.cos(self.lineOri)  + self.lineWidth * np.sin(self.angleMod(self.lineSkew - self.lineOri)) + self.hafX,
                        self.lineLength * np.sin(self.lineOri)  + self.lineWidth * np.cos(self.angleMod(self.lineSkew - self.lineOri)) + self.hafY)),
                        fill=255)
            self.imgPSF = np.asarray(pilPSF) / 255.

        print '\nIn Modefunc Psf mode = %s' % self.psfMode

        # self.update()
        plt.draw()

    def fourier_filter(self):

        # Fourier Transform
        fourImg = np.fft.fft2(self.imgHan)  # set dc term to 1 to control contrast
        fourImg[0, 0] = 1.0 + 0j
        self.fourShft = np.fft.fftshift(fourImg)
        self.fourLog = np.log(np.abs(self.fourShft))
        self.fourLog = self.fourLog / complex(self.fourLog.max())

        # take transform of PSF
        self.fourPSF = np.fft.fft2(self.imgPSF)
        fourShftPSF = np.fft.fftshift(self.fourPSF)
        self.psfLog = np.log(np.maximum(np.abs(fourShftPSF), 1.))
        self.psfLog = self.psfLog / complex(self.psfLog.max())

        # Create the Linear MAP filter, K(u,v)
        isnr = 1. / self.snr
        conjfourPSF = np.conj(self.fourPSF)
        self.K = (conjfourPSF + isnr) / ((conjfourPSF * self.fourPSF) + isnr)
        self.KLog = np.log(np.maximum(np.abs(self.K), 1.))
        self.KLog = self.KLog / complex(self.KLog.max())  # normalizing 0 to 1

        # do the inverse filtering
        self.fourResult = self.fourShft * self.K  # convolution in the fourier domain

        # Inverse Fourier Transform
        self.fourIshft = np.fft.ifftshift(self.fourResult)
        self.fourIshft[0, 0] = 0.5 + 0.0j  # set d.c. term for display
        self.fourInv = np.fft.ifft2(self.fourIshft)

        #   Swap quadrants
        self.fourInv = np.fft.ifftshift(self.fourInv)




    def press(self, event):

        if event.key == 'q':
            plt.close()

    def get_image(self):
        ''' get_image '''


        # Get image using finder dialog
        root = Tkinter.Tk()
        root.withdraw()  # Hide the root window
        imgFile = tkFileDialog.askopenfilename(title='Select Image',
                                               initialdir=".",
                                               initialfile='blr.png')

        return imgFile

    def angleMod(self, angle):

        return angle % (2 * np.pi)

    def update1(self, val):
        figPlot.radius     = self.slider1.val
        figPlot.lineLength = self.slider1.val
        self.modefunc(self.psfMode)
        self.update()

    def update2(self, val):
        figPlot.lineOri = self.slider2.val
        self.modefunc(self.psfMode)
        self.update()

    def update3(self, val):
        figPlot.snr = self.slider3.val
        self.modefunc(self.psfMode)
        self.update()

    def update4(self, val):
        figPlot.lineWidth = self.slider4.val
        self.modefunc(self.psfMode)
        self.update()


    def update5(self, val):
        figPlot.lineSkew = self.slider5.val
        self.modefunc(self.psfMode)
        self.update()


figPlot = pltFig()

plt.show(block=True)


