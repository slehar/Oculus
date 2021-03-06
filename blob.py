"""
blob.py
"""


import numpy as np
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
from matplotlib.path import Path
from   PIL import Image 


# global variables
ptList = []
patchList = []
inAPoint = False
selectedPt = None
scale = 1.
ptRad = 0.02
buttonState = False
altKeyPressed = False
shiftKeyPressed = False
curve4Pt = None
curveType = 'line'
nPts = 0
nAnchors = 0
fig = None
ax = None



# Keypress 'q' to quit
def keypress(event):
    global ptList, data, mute, altKeyPressed, shiftKeyPressed
    global curveType
    global fig, ax
    print 'event key = %s'%event.key    
    if event.key == 'alt':
        altKeyPressed = True
        print 'alt key pressed'
    if event.key == 'shift':
        shiftKeyPressed = True
        print 'shift key pressed'
    if event.key == '2':
        curveType = 'line'
        curveTypeText.set_text(curveType)
        print '"2" key pressed'
    if event.key == '3':
        curveType = 'curve3'
        curveTypeText.set_text(curveType)
        print '"3" key pressed'
        curveType = 'curve3'
    if event.key == '4':
        curveType = 'curve4'
        curveTypeText.set_text(curveType)
        print '"4" key pressed'
        curveType = 'curve4'
    if event.key == 'q':
        plt.close() 
    if event.key == 'return':
        print "RETURN key detected!"
        
    fig.canvas.draw()
    fig.canvas.mpl_connect('key_press_event', keypress)

# Key release
def keyrelease(event):
    global ptList, data, mute, altKeyPressed, shiftKeyPressed
    global fig, ax
        
    if event.key == 'alt':
        altKeyPressed = False
        print 'shift key UN-pressed'
    if event.key == 'shift':
        shiftKeyPressed = False
        print 'alt key UN-pressed'
    if event.key == '2':
        print '"2" key UN-pressed'
    if event.key == '3':
        print '"3" key UN-pressed'
    if event.key == '4':
        print '"4" key UN-pressed'
        
    fig.canvas.mpl_connect('key_release_event', keyrelease)

transMat = np.array([[scale,   0,     0],
                     [0,       scale, 0],
                     [0,       0,     1]])

invMat = np.linalg.inv(transMat)

def generatePolyPath(): 
    global patchList
    print '\nin generatePolyPath'
    
    path_data = []
    for pt in ptList:
        path_data.append(pt['pathcode'])
    for pth in path_data:
        print '  [%2d [%5.2f, %5.2f]]' % (pth[0], pth[1][0], pth[1][1])
    
    for pa in patchList:    # Remove previous polygon before adding new path
        pa.remove()
    patchList = []
    codes, verts = zip(*path_data)
    path = Path(verts, codes)
    patch = mpatches.PathPatch(path, facecolor='r', alpha=.5)
    patchList.append(patch)
    ax.add_patch(patch)
    fig.canvas.draw()

    
def appendPtMOVETO(xyLoc):
    global nPts, nAnchors
    xdata, ydata = xyLoc
    transPos     = [xdata, ydata, 1,]
    absPos       = np.matmul(transPos, invMat)
    circ         = mpatches.Circle(transPos, ptRad)
    ax.add_patch(circ)
    print 'in appendPtMOVETO [%5.2f, %5.2f]' % (xdata, ydata)
    xdata, ydata = xyLoc
    pathCode     = [Path.MOVETO, [xdata, ydata]]
    isAnchor     = False
    nPts += 1
    ptList.append({'xPos'     : transPos[0],
                   'yPos'     : transPos[1],
                   'absPos'   : absPos,
                   'transPos' : transPos,
                   'selected' : True,
                   'circle'   : circ,
                   'pathcode' : pathCode,
                   'isanchor' : isAnchor,
                   })
    
def appendPtLINETO(xyLoc):
    global nPts, nAnchors
    xdata, ydata = xyLoc
    transPos     = [xdata, ydata, 1,]
    absPos       = np.matmul(transPos, invMat)
    circ         = mpatches.Circle(transPos, ptRad)
    ax.add_patch(circ)
    print 'in appendPtLINETO [%5.2f, %5.2f]' % (xdata, ydata)
    xdata, ydata = xyLoc    
    pathCode     = [Path.LINETO, [xdata, ydata]]
    isAnchor     = False
    nPts += 1
    ptList.append({'xPos'     : transPos[0],
                   'yPos'     : transPos[1],
                   'absPos'   : absPos,
                   'transPos' : transPos,
                   'selected' : True,
                   'circle'   : circ,
                   'pathcode' : pathCode,
                   'isanchor' : isAnchor,
                   })

def appendPtLINECLOSE(xyLoc):
    global nPts, nAnchors
    xdata, ydata = xyLoc
    transPos     = [xdata, ydata, 1,]
    absPos       = np.matmul(transPos, invMat)
    circ         = mpatches.Circle(transPos, ptRad)
    ax.add_patch(circ)
    print 'in appendPtLINECLOSE [%5.2f, %5.2f]' % (xdata, ydata)
    xdata, ydata = xyLoc    
    pathCode     = [Path.LINETO, [xdata, ydata]]
    isAnchor     = False
    nPts += 1
    ptList.append({'xPos'     : transPos[0],
                   'yPos'     : transPos[1],
                   'absPos'   : absPos,
                   'transPos' : transPos,
                   'selected' : True,
                   'circle'   : circ,
                   'pathcode' : pathCode,
                   'isanchor' : isAnchor,
                   })
    pathCode = [Path.CLOSEPOLY, [xdata, ydata]]
    isAnchor = False
    ptList.append({'xPos'     : transPos[0],
                   'yPos'     : transPos[1],
                   'absPos'   : absPos,
                   'transPos' : transPos,
                   'selected' : False,
                   'circle'   : circ,
                   'pathcode' : pathCode,
                   'isanchor' : isAnchor,
                   })
                   
def appendPtLINEAPPEND(xyLoc):
    global nPts, nAnchors
    xdata, ydata = xyLoc
    transPos     = [xdata, ydata, 1,]
    absPos       = np.matmul(transPos, invMat)
    circ         = mpatches.Circle(transPos, ptRad)
    ax.add_patch(circ)
    print 'in appendPtLINEAPPEND [%5.2f, %5.2f]' % (xdata, ydata)
    xdata, ydata = xyLoc    
    pathCode     = [Path.LINETO, [xdata, ydata]]
    isAnchor     = False
    nPts += 1
    ptList.insert(-1, {'xPos'     : transPos[0],
                   'yPos'     : transPos[1],
                   'absPos'   : absPos,
                   'transPos' : transPos,
                   'selected' : True,
                   'circle'   : circ,
                   'pathcode' : pathCode,
                   'isanchor' : isAnchor,
                   })

#### CURVE3 ####

def setAnchors3(fromLoc, toLoc):
    (dx, dy) = (toLoc[0] - fromLoc[0],
                toLoc[1] - fromLoc[1])
    h = (fromLoc[0] + dx/2., fromLoc[1] + dy/2.)
    return h

def appendPtCURVE3TO(fromLoc, toLoc):
    global nPts, nAnchors
    (xdataFr, ydataFr)    = fromLoc
    (xdataTo, ydataTo)  = toLoc
    (handleX, handleY) = setAnchors3(fromLoc, toLoc)
                                          
    transPos1 = [handleX, handleY, 1,]
    absPos1   = np.matmul(transPos1, invMat)
    circ      = mpatches.Circle(transPos1, ptRad/2, fc='g')
    ax.add_patch(circ)
    print 'in appendPtCURVE3TO [%5.2f, %5.2f] [%5.2f, %5.2f]' % \
                              (xdataFr, ydataFr, xdataTo, ydataTo)
    pathCode  = [Path.CURVE3, [handleX, handleY]]
    isAnchor  = True
    nAnchors += 1
    ptList.append({'xPos'     : transPos1[0],
                   'yPos'     : transPos1[1],
                   'absPos'   : absPos1,
                   'transPos' : transPos1,
                   'selected' : True,
                   'circle'   : circ,
                   'pathcode' : pathCode,
                   'isanchor' : isAnchor,
                   })
                                      
    transPos2 = [xdataTo, ydataTo, 1,]
    absPos2   = np.matmul(transPos2, invMat)
    circ      = mpatches.Circle(transPos2, ptRad, fc='b')
    ax.add_patch(circ)
    print '  [Path.CURVE3 [%5.2f, %5.2f]]' % (xdataTo, ydataTo)
    pathCode  = [Path.CURVE3, [xdataTo, ydataTo]]
    isAnchor  = False
    nPts += 1
    ptList.append( {'xPos'     : transPos2[0],
                    'yPos'     : transPos2[1],
                    'absPos'   : absPos2,
                    'transPos' : transPos2,
                    'selected' : True,
                    'circle'   : circ,
                    'pathcode' : pathCode,
                    'isanchor' : isAnchor,
                   })

def appendPtCURVE3CLOSE(fromLoc, toLoc):
    global nPts, nAnchors
    (xdataFr, ydataFr) = fromLoc
    (xdataTo, ydataTo) = toLoc
    (xdataClose, ydataClose) = (ptList[0]['xPos'], ptList[0]['yPos'])
    closeLoc = (xdataClose, ydataClose)
    (handle1X, handle1Y) = setAnchors3(fromLoc, toLoc)
    (handle2X, handle2Y) = setAnchors3(toLoc,   closeLoc)

    transPos1 = [handle1X, handle1Y, 1, ]
    absPos1 = np.matmul(transPos1, invMat)
    circ = mpatches.Circle(transPos1, ptRad / 2, fc='g')
    ax.add_patch(circ)
    print 'in appendPtCURVE3CLOSE [%5.2f, %5.2f] [%5.2f, %5.2f]' % (
        xdataFr, ydataFr, xdataTo, ydataTo)
    print '  [Path.CURVE3 [%5.2f, %5.2f]]' % (handle1X, handle1Y)
    pathCode = [Path.CURVE3, [handle1X, handle1Y]]
    isAnchor = True
    nAnchors += 1
    ptList.append({'xPos': transPos1[0],
                   'yPos': transPos1[1],
                   'absPos': absPos1,
                   'transPos': transPos1,
                   'selected': True,
                   'circle': circ,
                   'pathcode': pathCode,
                   'isanchor': isAnchor,
                   })

    transPos2 = [xdataTo, ydataTo, 1, ]
    absPos2 = np.matmul(transPos2, invMat)
    circ = mpatches.Circle(transPos2, ptRad, fc='b')
    ax.add_patch(circ)
    print '  [Path.CURVE3 [%5.2f, %5.2f]]' % (xdataTo, ydataTo)
    pathCode = [Path.CURVE3, [xdataTo, ydataTo]]
    isAnchor = False
    nPts += 1
    ptList.append({'xPos': transPos2[0],
                   'yPos': transPos2[1],
                   'absPos': absPos2,
                   'transPos': transPos2,
                   'selected': True,
                   'circle': circ,
                   'pathcode': pathCode,
                   'isanchor': isAnchor,
                   })

    transPos3 = [handle2X, handle2Y, 1, ]
    absPos3 = np.matmul(transPos3, invMat)
    circ = mpatches.Circle(transPos3, ptRad / 2, fc='g')
    ax.add_patch(circ)
    print '  [Path.CURVE3 [%5.2f, %5.2f]]' % (handle2X, handle2Y)
    pathCode = [Path.CURVE3, [handle2X, handle2Y]]
    isAnchor = True
    nAnchors += 1
    ptList.append({'xPos': transPos3[0],
                   'yPos': transPos3[1],
                   'absPos': absPos3,
                   'transPos': transPos3,
                   'selected': True,
                   'circle': circ,
                   'pathcode': pathCode,
                   'isanchor': isAnchor,
                   })

    transPos4 = [xdataClose, ydataClose, 1, ]
    absPos4 = np.matmul(transPos4, invMat)
    circ = mpatches.Circle(transPos4, ptRad, fc='b')
    ax.add_patch(circ)
    print '  [Path.CURVE3 [%5.2f, %5.2f]]' % (xdataClose, ydataClose)
    pathCode = [Path.CURVE3, [xdataClose, ydataClose]]
    isAnchor = False
    nPts += 1
    ptList.append({'xPos': transPos4[0],
                   'yPos': transPos4[1],
                   'absPos': absPos4,
                   'transPos': transPos4,
                   'selected': True,
                   'circle': circ,
                   'pathcode': pathCode,
                   'isanchor': isAnchor,
                   })

    print '  [Path.CLOSEPOLY [%5.2f, %5.2f]]' % (xdataClose, ydataClose)
    pathCode = [Path.CLOSEPOLY, [xdataClose, ydataClose]]
    isAnchor = False
    #    nPts += 1
    ptList.append({'xPos': transPos4[0],
                   'yPos': transPos4[1],
                   'absPos': absPos4,
                   'transPos': transPos4,
                   'selected': True,
                   'circle': circ,
                   'pathcode': pathCode,
                   'isanchor': isAnchor,
                   })

def appendPtCURVE3APPEND(fromLoc, toLoc):
    global nPts, nAnchors

#    (xdata, ydata) = fromLoc
#    (xdata2, ydata2) = toLoc
#    (xdata0, ydata0) = (ptList[0]['xPos'], ptList[0]['yPos'])
#    (dx, dy) = (xdata2 - xdata, ydata2 - ydata)
#    (xdata1, ydata1) = ((xdata2 + xdata)  / 2, (ydata2 + ydata) / 2)
#    (xdata3, ydata3) = ((xdata2 + xdata0) / 2, (ydata2 + ydata) / 2)

    (xdataFr, ydataFr) = fromLoc
    (xdataTo, ydataTo) = toLoc
    (xdataClose, ydataClose) = (ptList[0]['xPos'], ptList[0]['yPos'])
    closeLoc = (xdataClose, ydataClose)
    (handle1X, handle1Y) = setAnchors3(fromLoc, toLoc)
    (handle2X, handle2Y) = setAnchors3(toLoc,   closeLoc)


    ptList[-3]['circle'].remove()
    del ptList[-3]
    nAnchors -= 1
    
    transPos1 = [handle1X, handle1Y, 1, ]
    absPos1 = np.matmul(transPos1, invMat)
    circ = mpatches.Circle(transPos1, ptRad / 2, fc='g')
    ax.add_patch(circ)
    print 'in appendPtCURVE3APPEND [%5.2f, %5.2f] [%5.2f, %5.2f]' % (
        xdataFr, ydataFr, xdataTo, ydataTo)
    print '  [Path.CURVE3 [%5.2f, %5.2f]]' % (handle1X, handle1Y)
    pathCode = [Path.CURVE3, [handle1X, handle1Y]]
    isAnchor = True
    nAnchors += 1
    ptList.insert(-2, {'xPos': transPos1[0],
                       'yPos': transPos1[1],
                       'absPos': absPos1,
                       'transPos': transPos1,
                       'selected': True,
                       'circle': circ,
                       'pathcode': pathCode,
                       'isanchor': isAnchor,
                       })

    transPos2 = [xdataTo, ydataTo, 1, ]
    absPos2 = np.matmul(transPos2, invMat)
    circ = mpatches.Circle(transPos2, ptRad, fc='b')
    ax.add_patch(circ)
    print '  [Path.CURVE3 [%5.2f, %5.2f]]' % (xdataTo, ydataTo)
    pathCode = [Path.CURVE3, [xdataTo, ydataTo]]
    isAnchor = False
    nPts += 1
    ptList.insert(-2, {'xPos': transPos2[0],
                       'yPos': transPos2[1],
                       'absPos': absPos2,
                       'transPos': transPos2,
                       'selected': True,
                       'circle': circ,
                       'pathcode': pathCode,
                       'isanchor': isAnchor,
                       })

    transPos3 = [handle2X, handle2Y, 1, ]
    absPos3 = np.matmul(transPos3, invMat)
    circ = mpatches.Circle(transPos3, ptRad / 2, fc='g')
    ax.add_patch(circ)
    print '  [Path.CURVE3 [%5.2f, %5.2f]]' % (handle2X, handle2Y)
    pathCode = [Path.CURVE3, [handle2X, handle2Y]]
    isAnchor = True
    nAnchors += 1
    ptList.insert(-2, {'xPos': transPos3[0],
                       'yPos': transPos3[1],
                       'absPos': absPos3,
                       'transPos': transPos3,
                       'selected': True,
                       'circle': circ,
                       'pathcode': pathCode,
                       'isanchor': isAnchor,
                       })


#### CURVE4 ####
                                                                           
def setAnchors4(fromLoc, toLoc):
    (dx, dy) = (toLoc[0] - fromLoc[0],
                toLoc[1] - fromLoc[1])
    h1 = (fromLoc[0] + dx*.25, fromLoc[1] + dy*.25)
    h2 = (fromLoc[0] + dx*.75, fromLoc[1] + dy*.75)
    return h1, h2

def appendPtCURVE4TO(fromLoc, toLoc):
    global nPts, nAnchors
    (xdataFr, ydataFr)  = fromLoc
    (xdataTo, ydataTo)  = toLoc
    print '\nin appendPtCURVE4TO'
    
    (handle1X, handle1Y), (handle2X, handle2Y)  = setAnchors4(fromLoc, toLoc)
    print '  fromLoc = [%5.2f, %5.2f]' % (xdataFr, ydataFr)
    print '  handle1 = [%5.2f, %5.2f]' % (handle1X, handle1Y)
    print '  handle2 = [%5.2f, %5.2f]' % (handle2X, handle2Y)
    print '  toLoc   = [%5.2f, %5.2f]' % (xdataTo, ydataTo)
                                             
    transPos1 = [handle1X, handle1Y, 1, ]
    absPos1   = np.matmul(transPos1, invMat)
    circ      = mpatches.Circle(transPos1, ptRad/2, fc='g')
    ax.add_patch(circ)
    print 'in appendPtCURVE4TO [%5.2f, %5.2f] [%5.2f, %5.2f]' % \
                              (xdataFr, ydataFr, handle2X, handle2Y)
    print '  [Path.CURVE4 [%5.2f, %5.2f]]' % (handle1X, handle1Y)
    pathCode  = [Path.CURVE4, [handle1X, handle1Y]]
    isAnchor  = True
    nAnchors += 1
    ptList.append({'xPos'    : handle1X,
                   'yPos'    : handle1Y,
                   'absPos'  : absPos1,
                   'transPos': transPos1,
                   'selected': True,
                   'circle'  : circ,
                   'pathcode': pathCode,
                   'isanchor': isAnchor,
                   })
                                      
    transPos2 = [handle2X, handle2Y, 1, ]
    absPos2 = np.matmul(transPos2, invMat)
    circ = mpatches.Circle(transPos2, ptRad / 2, fc='r')
    ax.add_patch(circ)
    print '  [Path.CURVE4 [%5.2f, %5.2f]]' % (handle2X, handle2Y,)
    pathCode = [Path.CURVE4, [handle2X, handle2Y,]]
    isAnchor = True
    nAnchors += 1
    ptList.append({'xPos': handle2X,
                   'yPos': handle2Y,
                   'absPos': absPos2,
                   'transPos': transPos2,
                   'selected': True,
                   'circle': circ,
                   'pathcode': pathCode,
                   'isanchor': isAnchor,
                   })

    transTo = [xdataTo, ydataTo, 1, ]
    absPosTo = np.matmul(transTo, invMat)
    circ = mpatches.Circle(transTo, ptRad, fc='b')
    ax.add_patch(circ)
    print '  [Path.CURVE4 [%5.2f, %5.2f]]' % (xdataTo, ydataTo)
    pathCode = [Path.CURVE4, [xdataTo, ydataTo]]
    isAnchor = False
    nPts += 1
    ptList.append({'xPos': xdataTo,
                   'yPos': ydataTo,
                   'absPos': absPosTo,
                   'transPos': transTo,
                   'selected': True,
                   'circle': circ,
                   'pathcode': pathCode,
                   'isanchor': isAnchor,
                   })

def appendPtCURVE4CLOSE(fromLoc, toLoc):
    global nPts, nAnchors
    (xdataFr, ydataFr) = fromLoc
    (xdataTo, ydataTo) = toLoc
    (xdataClose, ydataClose) = (ptList[0]['xPos'], ptList[0]['yPos'])
    closeLoc = (xdataClose, ydataClose)
    print '\nin appendPtCURVE4CLOSE'

    (handle1X, handle1Y), (handle2X, handle2Y) = setAnchors4(fromLoc, toLoc)
    print '  fromLoc = [%5.2f, %5.2f]' % (xdataFr, ydataFr)
    print '  handle1 = [%5.2f, %5.2f]' % (handle1X, handle1Y)
    print '  handle2 = [%5.2f, %5.2f]' % (handle2X, handle2Y)
    print '  toLoc   = [%5.2f, %5.2f]' % (xdataTo, ydataTo)

    (handle3X, handle3Y), (handle4X, handle4Y) = setAnchors4(toLoc, closeLoc)
    print '  handle3  = [%5.2f, %5.2f]' % (handle3X, handle3Y)
    print '  handle4  = [%5.2f, %5.2f]' % (handle4X, handle4Y)
    print '  closeLoc = [%5.2f, %5.2f]' % (closeLoc[0], closeLoc[1])

    transPos1 = [handle1X, handle1Y, 1, ]
    absPos1 = np.matmul(transPos1, invMat)
    circ = mpatches.Circle(transPos1, ptRad / 2, fc='g')
    ax.add_patch(circ)
    print '  [Path.CURVE4 [%5.2f, %5.2f]]' % (handle1X, handle1Y)
    pathCode = [Path.CURVE4, [handle1X, handle1Y]]
    isAnchor = True
    nAnchors += 1
    ptList.append({'xPos': handle1X,
                   'yPos': handle1Y,
                   'absPos': absPos1,
                   'transPos': transPos1,
                   'selected': True,
                   'circle': circ,
                   'pathcode': pathCode,
                   'isanchor': isAnchor,
                   })

    transPos2 = [handle2X, handle2Y, 1, ]
    absPos2 = np.matmul(transPos2, invMat)
    circ = mpatches.Circle(transPos2, ptRad / 2, fc='r')
    ax.add_patch(circ)
    print '  [Path.CURVE4 [%5.2f, %5.2f]]' % (handle2X, handle2Y,)
    pathCode = [Path.CURVE4, [handle2X, handle2Y,]]
    isAnchor = True
    nAnchors += 1
    ptList.append({'xPos': handle2X,
                   'yPos': handle2Y,
                   'absPos': absPos2,
                   'transPos': transPos2,
                   'selected': True,
                   'circle': circ,
                   'pathcode': pathCode,
                   'isanchor': isAnchor,
                   })

    transTo = [xdataTo, ydataTo, 1, ]
    absPosTo = np.matmul(transTo, invMat)
    circ = mpatches.Circle(transTo, ptRad, fc='b')
    ax.add_patch(circ)
    print '  [Path.CURVE4 [%5.2f, %5.2f]]' % (xdataTo, ydataTo)
    pathCode = [Path.CURVE4, [xdataTo, ydataTo]]
    isAnchor = False
    nPts += 1
    ptList.append({'xPos': xdataTo,
                   'yPos': ydataTo,
                   'absPos': absPosTo,
                   'transPos': transTo,
                   'selected': True,
                   'circle': circ,
                   'pathcode': pathCode,
                   'isanchor': isAnchor,
                   })
                   
    transPos3 = [handle3X, handle3Y, 1, ]
    absPos3 = np.matmul(transPos3, invMat)
    circ = mpatches.Circle(transPos3, ptRad / 2, fc='g')
    ax.add_patch(circ)
    print '  [Path.CURVE4 [%5.2f, %5.2f]]' % (handle3X, handle3Y)
    pathCode = [Path.CURVE4, [handle3X, handle3Y]]
    isAnchor = True
    nAnchors += 1
    ptList.append({'xPos': handle3X,
                   'yPos': handle3Y,
                   'absPos': absPos3,
                   'transPos': transPos3,
                   'selected': True,
                   'circle': circ,
                   'pathcode': pathCode,
                   'isanchor': isAnchor,
                   })

    transPos4 = [handle4X, handle4Y, 1, ]
    absPos4 = np.matmul(transPos4, invMat)
    circ = mpatches.Circle(transPos4, ptRad / 2, fc='r')
    ax.add_patch(circ)
    print '  [Path.CURVE4 [%5.2f, %5.2f]]' % (handle4X, handle4Y,)
    pathCode = [Path.CURVE4, [handle4X, handle4Y,]]
    isAnchor = True
    nAnchors += 1
    ptList.append({'xPos': handle4X,
                   'yPos': handle4Y,
                   'absPos': absPos4,
                   'transPos': transPos4,
                   'selected': True,
                   'circle': circ,
                   'pathcode': pathCode,
                   'isanchor': isAnchor,
                   })

    transClose = [xdataClose, ydataClose, 1, ]
    absPosClose = np.matmul(transClose, invMat)
    circ = mpatches.Circle(transClose, ptRad, fc='b')
    ax.add_patch(circ)
    print '  [Path.CURVE4 [%5.2f, %5.2f]]' % (xdataClose, ydataClose)
    pathCode = [Path.CURVE4, [xdataClose, ydataClose]]
    isAnchor = False
    nPts += 1
    ptList.append({'xPos': xdataClose,
                   'yPos': ydataClose,
                   'absPos': absPosClose,
                   'transPos': transClose,
                   'selected': True,
                   'circle': circ,
                   'pathcode': pathCode,
                   'isanchor': isAnchor,
                   })
                  
    (xdata0, ydata0) = (ptList[0]['xPos'], ptList[0]['yPos'])
    transPos0 = [xdata0, ydata0, 1]
    absPos0 = np.matmul(transPos0, invMat)
    print '  [Path.CLOSEPOLY [%5.2f, %5.2f]]' % (xdata0, ydata0)
    pathCode = [Path.CLOSEPOLY, [xdata0, ydata0]]
    isAnchor = False
    #    nPts += 1
    ptList.append({'xPos': xdata0,
                   'yPos': ydata0,
                   'absPos': absPos0,
                   'transPos': transPos0,
                   'selected': True,
                   'circle': circ,
                   'pathcode': pathCode,
                   'isanchor': isAnchor,
                   })

def appendPtCURVE4APPEND(fromLoc, toLoc):
    global nPts, nAnchors

    (xdata, ydata)    = fromLoc
    (xdata2, ydata2)  = toLoc
    (xdata0, ydata0)  = (ptList[0]['xPos'], ptList[0]['yPos'])
    (dx, dy) = (xdata2 - xdata, ydata2 - ydata)
    (xdata1, ydata1) = ((xdata2 + xdata)/2 + .2, (ydata2 + ydata)/2 +.2)
    (xdata3, ydata3) = ((xdata2 + xdata0)/2, (ydata2 + ydata)/2 +.2)
    
    ptList[-3]['circle'].remove()
    del ptList[-3]
    nAnchors -= 1
                                          
    transPos1 = [xdata1, ydata1, 1,]
    absPos1   = np.matmul(transPos1, invMat)
    circ      = mpatches.Circle(transPos1, ptRad/2, fc='g')
    ax.add_patch(circ)
    print 'in appendPtCURVE4APPEND [%5.2f, %5.2f] [%5.2f, %5.2f]' % (
                    xdata, ydata, xdata2, ydata2)
    print '  [Path.CURVE4 [%5.2f, %5.2f]]' % (xdata1, ydata1)
    pathCode  = [Path.CURVE4, [xdata1, ydata1]]
    isAnchor  = True
    nAnchors += 1
    ptList.insert(-2, {'xPos'     : transPos1[0],
                       'yPos'     : transPos1[1],
                       'absPos'   : absPos1,
                       'transPos' : transPos1,
                       'selected' : True,
                       'circle'   : circ,
                       'pathcode' : pathCode,
                       'isanchor' : isAnchor,
                       })
                                      
    transPos2 = [xdata2, ydata2, 1,]
    absPos2   = np.matmul(transPos2, invMat)
    circ      = mpatches.Circle(transPos2, ptRad, fc='b')
    ax.add_patch(circ)
    print '  [Path.CURVE4 [%5.2f, %5.2f]]' % (xdata2, ydata2)
    pathCode  = [Path.CURVE4, [xdata2, ydata2]]
    isAnchor  = False
    nPts += 1
    ptList.insert(-2, { 'xPos'     : transPos2[0],
                        'yPos'     : transPos2[1],
                        'absPos'   : absPos2,
                        'transPos' : transPos2,
                        'selected' : True,
                        'circle'   : circ,
                        'pathcode' : pathCode,
                        'isanchor' : isAnchor,
                       })

    transPos3 = [xdata3, ydata3, 1,]
    absPos3   = np.matmul(transPos3, invMat)
    circ      = mpatches.Circle(transPos3, ptRad/2, fc='g')
    ax.add_patch(circ)
    print '  [Path.CURVE4 [%5.2f, %5.2f]]' % (xdata3, ydata3)
    pathCode  = [Path.CURVE4, [xdata3, ydata1]]
    isAnchor  = True
    nAnchors += 1
    ptList.insert(-2, {'xPos'     : transPos3[0],
                       'yPos'     : transPos3[1],
                       'absPos'   : absPos1,
                       'transPos' : transPos1,
                       'selected' : True,
                       'circle'   : circ,
                       'pathcode' : pathCode,
                       'isanchor' : isAnchor,
                       })


############[ event handlers ]##############

def on_press(event):
    global selectedPt, buttonState
    global nPts, nAnchors, ptList
#    print 'in on_press nPts %d nAnchors %d' % (nPts, nAnchors)
    if event.inaxes is not ax:
        return
    inAPoint = False
    for pt in ptList:
        if pt['circle']:
            contains, attrd = pt['circle'].contains(event)
            if contains:
                inAPoint = True
                print 'inAPoint'
                if pt['isanchor']:
                    fc='g'
                else:
                    fc='b'
                if pt['selected']:
                    pt['selected'] = False
                    pt['circle'].set_fc(fc)
                else:
                    pt['selected'] = True
                    selectedPt = pt
                    pt['circle'].set_fc('red')
                break
    fig.canvas.draw()
    buttonState = True
        
    if not inAPoint:
        (xdata, ydata) = event.xdata, event.ydata
        transPos = [xdata, ydata, 1,]
        absPos = np.matmul(transPos, invMat)
        circ = mpatches.Circle(transPos, ptRad)
                        
        if curveType == 'line':
            if nPts == 0:
                appendPtMOVETO([xdata, ydata])  
            elif nPts == 1:
                    appendPtLINETO([xdata, ydata])  
            elif nPts == 2:
                    appendPtLINECLOSE([xdata, ydata])
            elif nPts > 2:
                    appendPtLINEAPPEND([xdata, ydata])
                    
        elif curveType == 'curve3':
            if nPts == 0:
                appendPtMOVETO([xdata, ydata])  
            elif nPts == 1:
                appendPtCURVE3TO([ptList[0]['xPos'], ptList[0]['yPos']],
                                 [xdata,             ydata]) 
            elif nPts == 2:
                appendPtCURVE3CLOSE([ptList[1+nAnchors]['xPos'],    ptList[1+nAnchors]['yPos']],
                                    [xdata,                         ydata])
            elif nPts > 2:
                (xdata0, ydata0) = (ptList[0]['xPos'], ptList[0]['yPos'])
                appendPtCURVE3APPEND([ptList[nPts+nAnchors-2]['xPos'], ptList[nPts+nAnchors-2]['yPos']],
                                    [xdata,                        ydata])
                            
        elif curveType == 'curve4':
            if nPts == 0:
                appendPtMOVETO([xdata, ydata])  
            elif nPts == 1:
                appendPtCURVE4TO([ptList[0]['xPos'], ptList[0]['yPos']],
                                 [xdata,             ydata]) 
            elif nPts == 2:
                (xdata0, ydata0) = (ptList[0]['xPos'], ptList[0]['yPos'])
                appendPtCURVE4CLOSE([ptList[1+nAnchors]['xPos'], ptList[1+nAnchors]['yPos']],
                                 [xdata,             ydata]) 
            elif nPts > 2:
#                (xdata0, ydata0) = (ptList[0]['xPos'], ptList[0]['yPos'])
                listLen = len(ptList)
                rmList = (listLen-1, listLen-3, listLen-4)
                for ix in rmList:
                    print '  remove %1d'%ix
                    ptList[ix]['circle'].remove()
                ptList = ptList[:-4]
                nPts -= 2
                nAnchors -= 2
                lastPt = nPts + nAnchors - 1
                appendPtCURVE4CLOSE([ptList[lastPt]['xPos'], ptList[lastPt]['yPos']],
                                     [xdata,            ydata])

        generatePolyPath()  
        fig.canvas.draw()
        selectedPt = ptList[-1]    
#        print 'out on_press nPts %d nAnchors %d' % (nPts, nAnchors)

def on_release(event):
    global buttonState, selectedPt
#    print 'in on_release'
    for pt in ptList:
        if pt['isanchor']:
            fc='g'
        else:
            fc='b'
        if pt['selected']:
            if shiftKeyPressed:
                (xdata, ydata) = (event.xdata, event.ydata)
                transPos = [xdata, ydata, 1,]
                pt['absPos'] = np.matmul(transPos, invMat)
                pt['circle'].center = transPos[:2]
            buttonState = False
            pt['selected'] = False
            selectedPt = None
            pt['circle'].set_fc(fc)
            fig.canvas.draw()
    buttonState = False

    
def on_motion(event):
    global xdata, ydata, selectedPt, ptList
    global fig
    
#    print 'in on_motion'
    xdata = event.xdata
    ydata = event.ydata
        
    if buttonState and shiftKeyPressed:
        absPos = [xdata, ydata, 1]
        transPos = np.matmul(absPos, transMat)
        selectedPt['circle'].center = transPos[:2]
        selectedPt['xPos'] = transPos[0]
        selectedPt['yPos'] = transPos[1]
        selectedPt['pathcode'][1] = [xdata, ydata]
        generatePolyPath()            
        fig.canvas.draw()

# Figure and axes
def openBlobWindow(imgPSF):
    global fig, ax
    fig = plt.figure(figsize=(6, 6))
    ax = fig.add_axes([.1, .1, .8, .8])
    ax.set_xlim([-1., 1.])
    ax.set_ylim([-1., 1.])
    ax.grid()
    curveTypeText = fig.text(.48, .95, curveType, size=18)
    fig.canvas.set_window_title('Blob')
    fig.canvas.mpl_connect('button_release_event', on_release)
    fig.canvas.mpl_connect('button_press_event',   on_press)
    fig.canvas.mpl_connect('motion_notify_event',  on_motion)
    
    plt.show()
    
    imgPSF = numpy.fromstring ( fig.canvas.tostring_argb(), dtype=numpy.float() )
        
    
def returnBlobImage():
    
    global fig, ax
    print '*** in blob ***'
    fig.canvas.draw()
    imgPil = Image.frombytes('RGB', fig.canvas.get_width_height(),
                             fig.canvas.tostring_rgb())
    imgNp = np.array(imgPil.convert('L'))/255.
    return imgNp
    
if __name__ == '__main__':
    openBlobWindow()
    
