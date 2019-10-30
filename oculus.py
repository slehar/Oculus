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

        # Disc / Line / Paragram Checkbox
        self.rax = plt.axes([0.41, 0.1, 0.15 / self.winAspect, 0.15])
        self.radio = RadioButtons(self.rax, ['Disc', 'Line', 'Paragram'])

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
        plt.pause(1.)

        # First Pass Hanning
        han = np.outer(np.hanning(self.ySize), np.hanning(self.xSize))
        self.imgHan = self.imgNp * han  # apply hanning window

        # Display hanning image
        self.beforePlot.set_data(self.imgHan)

        # Generate PSF image
        K = np.zeros(self.imgNp.shape)  # array for inverse filter
        yy, xx = np.mgrid[-self.hafY:self.hafY, -self.hafX:self.hafX]
        self.distImg = np.sqrt(xx ** 2 + yy ** 2)
        imgPSF = (self.distImg < self.radius)  # This is a Disc PSF

        # Display PSF image
        plt.sca(self.axPSF)
        self.psfPlot = plt.imshow(imgPSF, cmap='gray')

        # Fourier Transform
        fourImg = np.fft.fft2(self.imgHan)  # set dc term to 1 to control contrast
        fourImg[0, 0] = 1.0 + 0j
        self.fourShft = np.fft.fftshift(fourImg)
        self.fourLog = np.log(np.abs(self.fourShft))
        self.fourLog = self.fourLog / complex(self.fourLog.max())

        # Display Fourier image
        plt.sca(self.axFour)
        self.fourPlot = plt.imshow(self.fourLog.real, cmap='gray')

        # take transform of PSF
        fourPSF = np.fft.fft2(imgPSF)
        fourShftPSF = np.fft.fftshift(fourPSF)
        psfLog = np.log(np.maximum(np.abs(fourShftPSF), 1.))
        psfLog = psfLog / complex(psfLog.max())

        # Display PSF transform image
        # self.fourPlot.set_data(psfLog.real)

        # Create the Linear MAP filter, K(u,v)
        isnr = 1. / self.snr
        conjfourPSF = np.conj(fourPSF)
        K = (conjfourPSF + isnr) / ((conjfourPSF * fourPSF) + isnr)
        KLog = np.log(np.maximum(np.abs(K), 1.))
        KLog = KLog / complex(KLog.max())  # normalizing 0 to 1

        # Display KLog image
        plt.sca(self.axDiag)
        self.diagPlot = plt.imshow(KLog.real, cmap='gray')

        # do the inverse filtering
        fourResult = self.fourShft * K  # convolution in the fourier domain

        # Inverse Fourier Transform
        fourIshft = np.fft.ifftshift(fourResult)
        fourIshft[0, 0] = 0.5 + 0.0j  # set d.c. term for display
        self.fourInv = np.fft.ifft2(fourIshft)

        #   Swap quadrants
        self.fourInv = np.fft.ifftshift(self.fourInv)

        # Display After image
        plt.sca(self.axAfter)
        self.invPlot = plt.imshow(self.fourInv.real, cmap='gray')

        # Slider 1
        hafRadiusMax = min(self.ySize, self.xSize)  # Setting upper limit to radius focus blur
        self.slider1 = utils.get_slider(self.fig, [0.35, 0.5, 0.234, 0.04], 'radius', 0.5, hafRadiusMax / 10, valinit=5.)
        rad1 = np.log(self.slider1.val)  # make log control
        self.slider1.valtext.set_text(rad1)
        self.slider1.on_changed(self.update1)

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

        imgPSF = imgPSF.astype(float)
        self.psfPlot.set_data(imgPSF)

        self.radio.on_clicked(self.modefunc)

        self.fig.canvas.mpl_connect('key_press_event', self.press)


    def update(self):
        ''' Update display when sliders change '''

        linelen = self.lineLength
        linewidth = self.lineWidth
        lineori = self.lineOri
        skew    = self.lineSkew
        radius  = self.radius
        hafx, hafy = self.hafX, self.hafY

        # PSF in axPSF
        if self.psfMode == 'Disc':
            self.imgPSF = (self.distImg < radius)  # This is a Disc PSF
        elif self.psfMode == 'Line':
            self.lineLength = radius * 3.
            self.imgPSF = np.zeros([2 * hafy, 2 * hafx])
            self.pilPSF = Image.fromarray(self.imgPSF, 'L')
            draw = ImageDraw.Draw(self.pilPSF)
            draw.line(((linelen * -np.cos(lineori) + hafx,
                        linelen * -np.sin(lineori) + hafy,
                        linelen *  np.cos(lineori) + hafx,
                        linelen *  np.sin(lineori) + hafy)),
                        fill=255, width=self.lineWidth)
            self.imgPSF = np.asarray(self.pilPSF) / 255.
        elif self.psfMode == 'Paragram':
            linelen = self.radius
            self.imgPSF = np.zeros([2 * hafy, 2 * hafx])
            self.pilPSF = Image.fromarray(self.imgPSF, 'L')
            draw = ImageDraw.Draw(self.pilPSF)

            # Draw line across the  center length linelen orientation lineori
            # draw.line(((linelen * -np.cos(lineori) + hafx,
            #             linelen * -np.sin(lineori) + hafy,
            #             linelen *  np.cos(lineori) + hafx,
            #             linelen *  np.sin(lineori) + hafy)),
            #             fill=255, width=1)


            # skew = self.angleMod(lineori - skew)

            # print 'linori = %5.2f, skew = %5.2f, sin(%5.2f) = %5.2f' % (lineori, skew, np.sin(lineori - skew), skew)

            # Draw four lines from the two ends of the oriented line to the four points of the parallelogram
            # draw.line(((linelen * -np.cos(lineori) + hafx,
            #             linelen * -np.sin(lineori) + hafy,
            #             linelen * -np.cos(lineori) + linewidth * np.sin(self.angleMod(skew - lineori)) + hafx,
            #             linelen * -np.sin(lineori) + linewidth * np.cos(self.angleMod(skew - lineori)) + hafy)),
            #             fill=255, width=1)
            #
            # draw.line(((linelen * -np.cos(lineori) + hafx,
            #             linelen * -np.sin(lineori) + hafy,
            #             linelen * -np.cos(lineori) - linewidth * np.sin(self.angleMod(skew - lineori)) + hafx,
            #             linelen * -np.sin(lineori) - linewidth * np.cos(self.angleMod(skew - lineori)) + hafy)),
            #             fill=255, width=1)
            #
            # draw.line(((linelen * np.cos(lineori) + hafx,
            #             linelen * np.sin(lineori) + hafy,
            #             linelen * np.cos(lineori) + linewidth * np.sin(self.angleMod(skew - lineori)) + hafx,
            #             linelen * np.sin(lineori) + linewidth * np.cos(self.angleMod(skew - lineori)) + hafy)),
            #           fill=255, width=1)
            #
            # draw.line(((linelen * np.cos(lineori) + hafx,
            #             linelen * np.sin(lineori) + hafy,
            #             linelen * np.cos(lineori) - linewidth * np.sin(self.angleMod(skew - lineori)) + hafx,
            #             linelen * np.sin(lineori) - linewidth * np.cos(self.angleMod(skew - lineori)) + hafy)),
            #           fill=255, width=1)
            #

            draw.polygon((
                       (linelen * -np.cos(lineori) + linewidth * np.sin(self.angleMod(skew - lineori)) + hafx,
                        linelen * -np.sin(lineori) + linewidth * np.cos(self.angleMod(skew - lineori)) + hafy),
                       (linelen * -np.cos(lineori) - linewidth * np.sin(self.angleMod(skew - lineori)) + hafx,
                        linelen * -np.sin(lineori) - linewidth * np.cos(self.angleMod(skew - lineori)) + hafy),
                       (linelen * np.cos(lineori)  - linewidth * np.sin(self.angleMod(skew - lineori)) + hafx,
                        linelen * np.sin(lineori)  - linewidth * np.cos(self.angleMod(skew - lineori)) + hafy),
                       (linelen * np.cos(lineori)  + linewidth * np.sin(self.angleMod(skew - lineori)) + hafx,
                        linelen * np.sin(lineori)  + linewidth * np.cos(self.angleMod(skew - lineori)) + hafy)),
                        fill=255)

            self.imgPSF = np.asarray(self.pilPSF) / 255.

            # self.imgPSF = blob.returnBlobImage()

        realimgPSF = self.imgPSF.astype(float)
        self.psfPlot.set_data(realimgPSF)

        # take transform of psf and display it
        fourPSF = np.fft.fft2(self.imgPSF)
        self.fourShftPSF = np.fft.fftshift(fourPSF)
        psfLog = np.log(np.maximum(np.abs(self.fourShftPSF), 1.))
        psfLog = psfLog / complex(psfLog.max())

        # Create the Linear MAP filter, K(u,v)
        isnr = 1. / self.snr
        conjfourPSF = np.conj(fourPSF)
        K = (conjfourPSF + isnr) / ((conjfourPSF * fourPSF) + isnr)
        #    print K[hafY:hafY+1,hafX:hafX+1]    #is this the d.c. term?
        KLog = np.log(np.maximum(np.abs(K), 1.))
        KLog = KLog / complex(KLog.max())  # normalizing 0 to 1
        self.diagPlot.set_data(KLog.real)  # Plotting diag data
        norm = np.sum(K)  # for normalizing K

        # Display filtered Fourier image
        self.fourPlot.set_data(self.fourLog.real)

        # do the inverse filtering
        fourResult = self.fourShft * K  # convolution in the fourier domain

        # Inverse Fourier Transform
        fourIshft = np.fft.ifftshift(fourResult)
        fourIshft[0, 0] = 0.5 + 0.0j  # set d.c. term for display
        self.fourInv = np.fft.ifft2(fourIshft)

        # Swap quadrants
        self.fourInv = np.fft.ifftshift(self.fourInv)
        self.invPlot.set_data(self.fourInv.real)




    def modefunc(self, label):
        if label == 'Disc':
            self.psfMode = 'Disc'
            self.imgPSF = (self.distImg < self.radius)  # This is a Disc PSF
            self.slider1.label = 'radius'
        if label == 'Line':
            self.psfMode = 'Line'
            self.lineLength = self.radius * 3.
            self.imgPSF = np.zeros([2 * self.hafY, 2 * self.hafX])
            pilPSF = Image.fromarray(self.imgPSF, 'L')
            self.slider1.label = 'length'
            print self.slider1.label
            draw = ImageDraw.Draw(pilPSF)
            draw.line(((
                self.lineLength * -np.cos(self.lineOri) + self.hafX,
                self.lineLength * -np.sin(self.lineOri) + self.hafY,
                self.lineLength * np.cos(self.lineOri) + self.hafX,
                self.lineLength * np.sin(self.lineOri) + self.hafY)),
                fill=255, width=self.lineWidth)
            self.imgPSF = np.asarray(pilPSF) / 255.

        if label == 'Paragram':
            self.psfMode = 'Paragram'
            # blobFig, blobAx = blob.openBlobWindow()
            # self.fig.canvas.draw()
            # self.imgPSF = blob.returnBlobImage()

            self.imgPSF = np.zeros([2 * self.hafY, 2 * self.hafX])
            self.imgPSF = self.imgPSF.astype(float)


        self.update()
        plt.draw()

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
        figPlot.radius = self.slider1.val
        self.update()

    def update2(self, val):
        figPlot.lineOri = self.slider2.val
        self.update()

    def update3(self, val):
        figPlot.snr = self.slider3.val
        self.update()

    def update4(self, val):
        figPlot.lineWidth = self.slider4.val
        self.update()


    def update5(self, val):
        figPlot.lineSkew = self.slider5.val
        self.update()


figPlot = pltFig()

plt.show()


