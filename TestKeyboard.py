# -*- coding: utf-8 -*-
"""
TestKeyboard.py

Created on Fri Sep 22 14:37:47 2017

@author: slehar
"""


import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt

fig = plt.figure()
ax = fig.add_axes([.1, .1, .8, .8])

# Keypress 'q' to quit
def keypress(event):
    print 'event key = %s'%event.key
    if event.key == 'q':
        plt.close()
fig.canvas.mpl_connect('key_press_event', keypress)

#
plt.show()
