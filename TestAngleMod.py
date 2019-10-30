""" TestAngleMod.py """
''' Test the AngleMod() function '''


import numpy as np

def AngleMod(angle):

    sign = np.sign(angle)
    return angle % (2 * np.pi) * sign

for a1 in np.linspace(-np.pi, np.pi, 16):
    print '\na1 = %5.2f pi' % (a1/np.pi)
    for a2 in np.linspace(-np.pi, np.pi, 16):
        print 'a1 %5.2f a2 %5.2f sum %5.2f AngleMod %5.2f pi' % (a1/np.pi, a2/np.pi, (a1+a2)/np.pi, AngleMod(a1+a2)/np.pi)

