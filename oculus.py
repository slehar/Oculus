# -*- coding: utf-8 -*-
"""
oculus.py

Created on Wed Sep 28 16:36:45 2016

@author: slehar
"""

import matplotlib
matplotlib.use('TkAgg')

import matplotlib.pyplot as plt
from   matplotlib.widgets import RadioButtons
# import Tkinter, tkFileDialog
from   PIL import Image, ImageDraw
import numpy as np
import utils
import blob
# import utils

class pltFig():
    ''' Plot Figure class '''

    radius = 5.
    snr = 100.
    isnr=1./snr
    slidersLocked = False
    psfMode = 'Disc'
    lineOri = 0.
    lineLength = 100.
    lineWidth = 5

    def __init__(self):

        # Open figure window
        winXSize = 15  # approx inches
        winYSize = 8   # in width & height
        self.winAspect = winXSize / winYSize
        self.fig = plt.figure(figsize=(winXSize, winYSize))  # Open the Fig
        self.fig.canvas.set_window_title('Oculus')
        # self.fig.canvas.mpl_connect('key_press_event', self.press)

        # Axes for Images
        self.axBefore = utils.get_axes(self.fig, [.05, .6, .35 / self.winAspect, .35], 'Before')
        self.axFour   = utils.get_axes(self.fig, [.65, .6, .35 / self.winAspect, .35], 'Fourier')
        self.axPSF    = utils.get_axes(self.fig, [.05, .2, .35 / self.winAspect, .35], 'PSF')
        self.axAfter  = utils.get_axes(self.fig, [.35, .6, .35 / self.winAspect, .35], 'After')
        self.axDiag   = utils.get_axes(self.fig, [.65, .2, .35 / self.winAspect, .35], 'Diagnostic')

        # Disc / Line / blob Checkbox
        self.rax = plt.axes([0.41, 0.1, 0.06 / self.winAspect, 0.1])
        self.radio = RadioButtons(self.rax, ['Disc', 'Line', 'Blob'])


        # Read image and display
        imgFile = self.get_image()
        imgPil = Image.open(imgFile).convert('LA')  # LA = 'luminance + alpha
        self.imgNp = np.array(imgPil.convert('L')) / 255.  # convert to L, luminance only
        self.ySize, self.xSize = self.imgNp.shape
        self.hafY, self.hafX = int(self.ySize / 2), int(self.xSize / 2)

        imgNp = utils.resize_even_dim(self.imgNp)
        self.ySize, self.xSize = imgNp.shape

        # First Pass Hanning
        han = np.outer(np.hanning(self.ySize), np.hanning(self.xSize))
        imgHan = self.imgNp * han  # apply hanning window
        totalpixels = self.xSize * self.ySize

        plt.sca(self.axDiag)
        imgplot = plt.imshow(self.imgNp, cmap='gray')
        x, y = self.hafX, self.hafY

        K = np.zeros(self.imgNp.shape)  # array for inverse filter
        yy, xx = np.mgrid[-self.hafY:self.hafY, -self.hafX:self.hafX]
        self.distImg = np.sqrt(xx ** 2 + yy ** 2)
        imgPSF = (self.distImg < self.radius)  # This is a Disc PSF

        self.psfPlot = plt.imshow(imgPSF, cmap='gray')
        plt.sca(self.axFour)

        # Fourier Transform
        fourImg = np.fft.fft2(imgHan)  # set dc term to 1 to control contrast
        fourImg[0, 0] = 1.0 + 0j
        self.fourShft = np.fft.fftshift(fourImg)
        fourLog = np.log(np.abs(self.fourShft))
        fourLog = fourLog / complex(fourLog.max())

        plt.sca(self.axFour)
        self.fourPlot = plt.imshow(fourLog.real, cmap='gray')
        # plt.pause(.001)

        # Slider 1
        hafRadiusMax = min(self.ySize, self.xSize)  # Setting upper limit to radius focus blur
        self.slider1 = utils.get_slider(self.fig, [0.41, 0.5, 0.234, 0.04], 'radius', 0.5, hafRadiusMax / 10, valinit=5.)
        rad1 = np.log(self.slider1.val)  # make log control
        self.slider1.valtext.set_text(rad1)

        # Slider 2
        self.slider2 = utils.get_slider(self.fig, [0.41, 0.4509, 0.234, 0.04], 'angle', -np.pi, np.pi, valinit=0.)
        self.lineOri = self.slider2.val

        # Slider 3
        self.slider3 = utils.get_slider(self.fig, [0.41, 0.40, 0.234, 0.04], 'SNR', 1., 1000., valinit=100.)
        self.snr = self.slider3.val

        # Slider 4
        self.slider4 = utils.get_slider(self.fig, [0.41, 0.35, 0.234, 0.04], 'linewidth', 1, 50, valinit=5)
        self.lineWidth = self.slider4.val

        # Slider 5
        self.slider5 = utils.get_slider(self.fig, [0.41, 0.3, 0.234, 0.04], 'skew', 1, 50, valinit=0)
        skew = self.slider5.val

        plt.sca(self.axBefore)
        beforePlot = plt.imshow(imgNp, cmap='gray', vmin=0., vmax=1.)

        plt.sca(self.axDiag)
        self.diagPlot = plt.imshow(imgHan, cmap='gray', vmin=0., vmax=1.)

        imgPSF = imgPSF.astype(float)
        plt.sca(self.axPSF)  # Set imgPSF the "current axes"
        self.psfPlot = plt.imshow(imgPSF, cmap='gray')

        self.radio.on_clicked(self.modefunc)

        self.fig.canvas.mpl_connect('key_press_event', self.press)

    def modefunc(self, label):
        if label == 'Disc':
            psfMode = 'Disc'
            imgPSF = (self.distImg < figPlot.radius)  # This is a Disc PSF
            self.slider1.label = 'radius'
        if label == 'Line':
            figPlot.psfMode = 'Line'
            figPlot.lineLength = figPlot.radius * 3.
            imgPSF = np.zeros([2 * self.hafY, 2 * self.hafX])
            pilPSF = Image.fromarray(imgPSF, 'L')
            self.slider1.label = 'length'
            print self.slider1.label
            draw = ImageDraw.Draw(pilPSF)
            draw.line(((
                figPlot.lineLength * -np.cos(figPlot.lineOri) + self.hafX,
                figPlot.lineLength * -np.sin(figPlot.lineOri) + self.hafY,
                figPlot.lineLength * np.cos(figPlot.lineOri) + self.hafX,
                figPlot.lineLength * np.sin(figPlot.lineOri) + self.hafY)),
                fill=255, width=figPlot.lineWidth)
            imgPSF = np.asarray(pilPSF) / 255.

        if label == 'Blob':
            figPlot.psfMode = 'Blob'
            blobFig, blobAx = blob.openBlobWindow()
            figPlot.fig.canvas.draw()
            imgPSF = blob.returnBlobImage()

        imgPSF = imgPSF.astype(float)
        plt.show()
        plt.pause(.001)
        self.update()
        plt.draw()


    # Keypress 'q' to quit callback function
    def press(self, event):

        print 'Keypress detected'
        if event.key == 'q':
            print "key is 'q'"
            plt.close()


        # plt.sca(self.axAfter)
        # invPlot = plt.imshow(fourReal, cmap='gray')


    def get_image(self):
        ''' get_image '''

        import Tkinter, tkFileDialog

        # Get image using finder dialog
        root = Tkinter.Tk()
        root.withdraw()  # Hide the root window
        imgFile = tkFileDialog.askopenfilename(title='Select Image',
                                               initialdir=".",
                                               initialfile='blr.png')

        return imgFile





    #########################[ def update(): ]########################################
    # This is loop where all the action happens
    def update(self):

        # PSF in axPSF
        if figPlot.psfMode == 'Disc':
            self.imgPSF = (self.distImg < figPlot.radius)  # This is a Disc PSF
        elif figPlot.psfMode == 'Line':
            figPlot.lineLength = figPlot.radius * 3.
            imgPSF = np.zeros([2 * self.hafY, 2 * self.hafX])
            pilPSF = Image.fromarray(imgPSF, 'L')
            draw = ImageDraw.Draw(pilPSF)
            draw.line(((figPlot.lineLength * -np.cos(figPlot.lineOri) + self.hafX,
                        figPlot.lineLength * -np.sin(figPlot.lineOri) + self.hafY,
                        figPlot.lineLength *  np.cos(figPlot.lineOri) + self.hafX,
                        figPlot.lineLength *  np.sin(figPlot.lineOri) + self.hafY)),
                        fill=255, width=figPlot.lineWidth)
            imgPSF = np.asarray(pilPSF) / 255.
        elif figPlot.psfMode == 'Blob':
            imgPSF = blob.returnBlobImage()

        realimgPSF = self.imgPSF.astype(float)
        self.psfPlot.set_data(realimgPSF)

        # take transform of psf and displaying it
        fourPSF = np.fft.fft2(self.imgPSF)
        self.fourShftPSF = np.fft.fftshift(fourPSF)
        psfLog = np.log(np.maximum(np.abs(self.fourShftPSF), 1.))
        psfLog = psfLog / complex(psfLog.max())
        self.fourPlot.set_data(psfLog.real)

        # Create the Linear MAP filter, K(u,v)
        isnr = 1. / figPlot.snr
        # -------
        conjfourPSF = np.conj(fourPSF)
        K = (conjfourPSF + isnr) / ((conjfourPSF * fourPSF) + isnr)
        #    print K[hafY:hafY+1,hafX:hafX+1]    #is this the d.c. term?
        KLog = np.log(np.maximum(np.abs(K), 1.))
        KLog = KLog / complex(KLog.max())  # normalizing 0 to 1
        self.diagPlot.set_data(KLog.real)  # Plotting diag data
        norm = np.sum(K)  # for normalizing K

        # do the inverse filtering
        fourResult = self.fourShft * K  # convolution in the fourier domain

        # Inverse Fourier Transform
        fourIshft = np.fft.ifftshift(fourResult)
        fourIshft[0, 0] = 0.5 + 0.0j  # set d.c. term for display

        fourInv = np.fft.ifft2(fourIshft)

        #   make sure fourReal scales 0.to 1.0 for display
        fourInv = np.fft.ifftshift(fourInv)
        fourReal = np.real(fourInv)

        #    lmax=fourReal.max()
        #    lmin=fourReal.min()

        plt.sca(self.axAfter)
        invPlot = plt.imshow(fourReal, cmap='gray')


     # radio button callback function to switch PSF mode


    # %%
    def update1(self, val):
        figPlot.radius = self.slider1.val
        self.update()


    def update2(self, val):
        figPlot.lineOri = self.slider2.val
        self.update()


    def update3(self, val):
        figPlot.snr = self.slider3.val
        self.update()


    def update4(self, val):
        figPlot.lineWidth = int(self.slider4.val)
        self.update()


        #    fig.canvas.draw()
        self.slider1.on_changed(self.update1)
        self.slider2.on_changed(self.update2)
        self.slider3.on_changed(self.update3)
        self.slider4.on_changed(self.update4)


figPlot = pltFig()
figPlot.update()

# Show [block = leave it open, don't close]
plt.show(block=True)
# plt.show()


