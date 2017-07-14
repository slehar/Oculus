"""
Demo of an interactive PathPatch object.
See 
https://matplotlib.org/examples/shapes_and_collections/path_patch_demo.html
"""
import numpy as np
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
from matplotlib.path import Path


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

# Figure and axes
plt.close('all')
fig = plt.figure(figsize=(6, 6))
ax = fig.add_axes([.1, .1, .8, .8])
ax.set_xlim([-1., 1.])
ax.set_ylim([-1., 1.])
ax.grid()
#ax.axis('equal')

curveTypeText = fig.text(.48, .95, curveType, size=18)

# Keypress 'q' to quit
def keypress(event):
    global ptList, data, mute, altKeyPressed, shiftKeyPressed
    global curveType
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
    fig.canvas.draw()
fig.canvas.mpl_connect('key_press_event', keypress)

# Key release
def keyrelease(event):
    global ptList, data, mute, altKeyPressed, shiftKeyPressed
        
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
#        if pt['isanchor']:
#            pathcode = [path.MOVETO(pt)]
#    print path_data
    for pth in path_data:
        print '  [%2d [%5.2f, %5.2f]]' % (pth[0], pth[1][0], pth[1][1])
    print ax.patches
    
#    for pat in ax.patches:
#        pat.remove()
    
    for pa in patchList:    # Remove previous polygon before adding new path
        pa.remove()
    patchList = []
    codes, verts = zip(*path_data)
    path = Path(verts, codes)
    patch = mpatches.PathPatch(path, facecolor='none', visible=True)
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
    (xdata0, ydata0) = (ptList[0]['xPos'], ptList[0]['yPos'])             
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
    print 'in appendPtAPPENDCLOSE [%5.2f, %5.2f]' % (xdata, ydata)
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

def appendPtCURVE3TO(fromLoc, toLoc):
    global nPts, nAnchors
    (xdata, ydata)    = fromLoc
    (xdata2, ydata2)  = toLoc
    print 'toLoc = %r' % toLoc
    (dx, dy) = (xdata2 - xdata, ydata2 - ydata)
    (xdata1, ydata1) = ((xdata + xdata2)/2, ydata+.2)
                                          
    transPos1 = [xdata1, ydata1, 1,]
    absPos1   = np.matmul(transPos1, invMat)
    circ      = mpatches.Circle(transPos1, ptRad/2, fc='g')
    ax.add_patch(circ)
    print 'in appendPtCURVE3TO [%5.2f, %5.2f]' % (xdata1, ydata1)
    pathCode  = [Path.CURVE3, [xdata1, ydata1]]
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
                                      
    transPos2 = [xdata2, ydata2, 1,]
    absPos2   = np.matmul(transPos2, invMat)
    circ      = mpatches.Circle(transPos2, ptRad, fc='b')
    ax.add_patch(circ)
    print 'in appendPtCURVE3TO [%5.2f, %5.2f]' % (xdata2, ydata2)
    pathCode  = [Path.CURVE3, [xdata2, ydata2]]
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
    (xdata, ydata)    = fromLoc
    (xdata2, ydata2)  = toLoc
    (dx, dy) = (xdata2 - xdata, ydata2 - ydata)
    (xdata1, ydata1) = ((xdata + xdata2)/2, ydata+.2)
                                          
    transPos = [xdata, ydata, 1,]
    absPos   = np.matmul(transPos, invMat)
    print '  [Path.MOVETO [%5.2f, %5.2f]]' % (xdata, ydata)
    pathCode     = [Path.MOVETO, [xdata, ydata]]
    isAnchor     = False
    nPts += 1
    ptList.append({'xPos'     : transPos[0],
                   'yPos'     : transPos[1],
                   'absPos'   : absPos,
                   'transPos' : transPos,
                   'selected' : True,
                   'circle'   : None,
                   'pathcode' : pathCode,
                   'isanchor' : isAnchor,
                   })

    transPos1 = [xdata1, ydata1, 1,]
    absPos1   = np.matmul(transPos1, invMat)
    circ      = mpatches.Circle(transPos1, ptRad/2, fc='g')
    ax.add_patch(circ)
    print 'in appendPtCURVE3CLOSE [%5.2f, %5.2f] [%5.2f, %5.2f]' % (
                    xdata, ydata, xdata1, ydata2)
    print '  [Path.CURVE3 [%5.2f, %5.2f]]' % (xdata1, ydata1)
    pathCode  = [Path.CURVE3, [xdata1, ydata1]]
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
                                      
    transPos2 = [xdata2, ydata2, 1,]
    absPos2   = np.matmul(transPos2, invMat)
    circ      = mpatches.Circle(transPos2, ptRad, fc='b')
    ax.add_patch(circ)
    print '  [Path.CURVE3 [%5.2f, %5.2f]]' % (xdata2, ydata2)
    pathCode  = [Path.CURVE3, [xdata2, ydata2]]
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

    transPos     = [ptList[0]['xPos'], ptList[0]['yPos'], 1,]
    absPos       = np.matmul(transPos, invMat)
    (xdata0, ydata0) = (ptList[0]['xPos'], ptList[0]['yPos'])
    print '  [Path.CLOSEPOLY [%5.2f, %5.2f]]' % (xdata0, ydata0)
    pathCode = [Path.CLOSEPOLY, [xdata0, ydata0]]
    circle = None
    isAnchor = False
    nPts += 1
    ptList.append({'xPos'     : transPos[0],
                   'yPos'     : transPos[1],
                   'absPos'   : absPos,
                   'transPos' : transPos,
                   'selected' : False,
                   'circle'   : circ,
                   'pathcode' : pathCode,
                   'isanchor' : isAnchor,
                   })
                   
def appendPtCURVE3APPEND(fromLoc, toLoc):
    global nPts, nAnchors
    (xdata, ydata)    = fromLoc
    (xdata2, ydata2)  = toLoc
    (dx, dy) = (xdata2 - xdata, ydata2 - ydata)
    (xdata1, ydata1) = ((xdata + xdata2)/2. + .2, (ydata + ydata2)/2. + .2)

    print 'in appendPtCURVE3APPEND [%5.2f, %5.2f] [%5.2f, %5.2f]' % (
                    xdata, ydata, xdata2, ydata2)


    transPos1 = [xdata1, ydata1, 1,]
    absPos1   = np.matmul(transPos1, invMat)
    circ      = mpatches.Circle(transPos1, ptRad/2, fc='g')
    ax.add_patch(circ)
    print '  [Path.CURVE3 [%5.2f, %5.2f]]' % (xdata1, ydata1)
    pathCode  = [Path.CURVE3, [xdata1, ydata1]]
    isAnchor  = True
    nAnchors += 1
    ptList.insert(-1, {'xPos'     : transPos1[0],
                   'yPos'     : transPos1[1],
                   'absPos'   : absPos1,
                   'transPos' : transPos1,
                   'selected' : True,
                   'circle'   : circ,
                   'pathcode' : pathCode,
                   'isanchor' : isAnchor,
                   })
                                      
#    transPos2 = [xdata2, ydata2, 1,]
#    absPos2   = np.matmul(transPos2, invMat)
#    circ      = mpatches.Circle(transPos2, ptRad, fc='b')
#    ax.add_patch(circ)
#    print 'in appendPtCURVE3APPEND [%5.2f, %5.2f]' % (xdata2, ydata2)
#    pathCode  = [Path.CURVE3, [xdata2, ydata2]]
#    isAnchor  = False
#    ptList.append( {'xPos'     : transPos2[0],
#                    'yPos'     : transPos2[1],
#                    'absPos'   : absPos2,
#                    'transPos' : transPos2,
#                    'selected' : True,
#                    'circle'   : circ,
#                    'pathcode' : pathCode,
#                    'isanchor' : isAnchor,
#                   })

#    transPos     = [xdata, ydata, 1,]
#    absPos       = np.matmul(transPos, invMat)
#    print '  [Path.CURVE3 [%5.2f, %5.2f]]' % (xdata1, ydata1)
#    pathCode = [Path.CLOSEPOLY, path_data[0][1]]
#    isAnchor = False
#    ptList.insert(-0, {'xPos'     : transPos[0],
#                       'yPos'     : transPos[1],
#                       'absPos'   : absPos,
#                       'transPos' : transPos,
#                       'selected' : False,
#                       'circle'   : circ,
#                       'pathcode' : pathCode,
#                       'isanchor' : isAnchor,
#                   })

def appendPtCURVE4TO(xyLoc):
    global nPts
    xdata,  ydata  = xyLoc
    xdata1, ydata1 =  [xdata + .2, ydata + .2]
    xdata2, ydata2 =  [xdata + .4, ydata     ]
    xdata3, ydata3 =  [xdata + .2, ydata - .2]
    
    transPos = [xdata, ydata, 1,]
    absPos   = np.matmul(transPos, invMat)
    circ     = mpatches.Circle(transPos, ptRad)
    ax.add_patch(circ)
    print 'in appendPtCURVE4 [%5.2f, %5.2f]' % (xdata, ydata)
    pathCode = [Path.MOVETO, [xdata, ydata]]
    isAnchor = False
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
                   
    transPos1 = [xdata1, ydata1, 1,]
    absPos1   = np.matmul(transPos1, invMat)
    circ      = mpatches.Circle(transPos1, ptRad/2, fc='g')
    ax.add_patch(circ)
    print 'in appendPtCURVE4 [%5.2f, %5.2f]' % (xdata1, ydata1)
    pathCode = [Path.CURVE4, [xdata1, ydata1]]
    isAnchor = True
    ptList.append( {'xPos'     : transPos1[0],
                    'yPos'     : transPos1[1],
                    'absPos'   : absPos,
                    'transPos' : transPos,
                    'selected' : True,
                    'circle'   : circ,
                    'pathcode' : pathCode,
                    'isanchor' : isAnchor,
                   })
                   
    transPos2 = [xdata2, ydata2, 1,]
    absPos2   = np.matmul(transPos2, invMat)
    circ      = mpatches.Circle(transPos2, ptRad/2, fc='g')
    ax.add_patch(circ)
    print 'in appendPtCURVE4 [%5.2f, %5.2f]' % (xdata2, ydata2)
    pathCode  = [Path.CURVE4, [xdata2, ydata2]]
    isAnchor  = True
    ptList.append( {'xPos'     : transPos2[0],
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
    circ      = mpatches.Circle(transPos3, ptRad, fc='b')
    ax.add_patch(circ)
    print 'in appendPtCURVE4 [%5.2f, %5.2f]' % (xdata3, ydata3)
    pathCode  = [Path.CURVE4, [xdata3, ydata3]]
    isAnchor  = False
    nPts += 1
    ptList.append( {'xPos'     : transPos3[0],
                    'yPos'     : transPos3[1],
                    'absPos'   : absPos3,
                    'transPos' : transPos3,
                    'selected' : True,
                    'circle'   : circ,
                    'pathcode' : pathCode,
                    'isanchor' : isAnchor,
                   })
    
def appendPtCURVE4CLOSE(xyLoc):
    global nPts
    xdata,  ydata  = xyLoc
    xdata1, ydata1 =  [xdata + .2, ydata + .2]
    xdata2, ydata2 =  [xdata + .4, ydata     ]
    xdata3, ydata3 =  [xdata + .2, ydata - .2]
    
    transPos = [xdata, ydata, 1,]
    absPos   = np.matmul(transPos, invMat)
    circ     = mpatches.Circle(transPos, ptRad)
    ax.add_patch(circ)
    print 'in appendPtCURVE4 [%5.2f, %5.2f]' % (xdata, ydata)
    pathCode = [Path.MOVETO, [xdata, ydata]]
    isAnchor = False
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
                   
    transPos1 = [xdata1, ydata1, 1,]
    absPos1   = np.matmul(transPos1, invMat)
    circ      = mpatches.Circle(transPos1, ptRad/2, fc='g')
    ax.add_patch(circ)
    print 'in appendPtCURVE4 [%5.2f, %5.2f]' % (xdata1, ydata1)
    pathCode = [Path.CURVE4, [xdata1, ydata1]]
    isAnchor = True
    ptList.append( {'xPos'     : transPos1[0],
                    'yPos'     : transPos1[1],
                    'absPos'   : absPos,
                    'transPos' : transPos,
                    'selected' : True,
                    'circle'   : circ,
                    'pathcode' : pathCode,
                    'isanchor' : isAnchor,
                   })
                   
    transPos2 = [xdata2, ydata2, 1,]
    absPos2   = np.matmul(transPos2, invMat)
    circ      = mpatches.Circle(transPos2, ptRad/2, fc='g')
    ax.add_patch(circ)
    print 'in appendPtCURVE4 [%5.2f, %5.2f]' % (xdata2, ydata2)
    pathCode  = [Path.CURVE4, [xdata2, ydata2]]
    isAnchor  = True
    ptList.append( {'xPos'     : transPos2[0],
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
    circ      = mpatches.Circle(transPos3, ptRad, fc='b')
    ax.add_patch(circ)
    print 'in appendPtCURVE4 [%5.2f, %5.2f]' % (xdata3, ydata3)
    pathCode  = [Path.CURVE4, [xdata3, ydata3]]
    isAnchor  = False
    nPts += 1
    ptList.append( {'xPos'     : transPos3[0],
                    'yPos'     : transPos3[1],
                    'absPos'   : absPos3,
                    'transPos' : transPos3,
                    'selected' : True,
                    'circle'   : circ,
                    'pathcode' : pathCode,
                    'isanchor' : isAnchor,
                   })
    
def appendPtCURVE4APPEND(xyLoc):
    global nPts
    xdata,  ydata  = xyLoc
    xdata2, ydata2 =  [ptList[0]['xLoc']], [ptList[0]['Loc']]
    
    xdata1, ydata1 =  [xdata + .2, ydata + .2]
    xdata3, ydata3 =  [xdata + .2, ydata - .2]
    
    transPos = [xdata, ydata, 1,]
    absPos   = np.matmul(transPos, invMat)
    circ     = mpatches.Circle(transPos, ptRad)
    ax.add_patch(circ)
    print 'in appendPtCURVE4 [%5.2f, %5.2f]' % (xdata, ydata)
    pathCode = [Path.MOVETO, [xdata, ydata]]
    isAnchor = False
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
                   
    transPos1 = [xdata1, ydata1, 1,]
    absPos1   = np.matmul(transPos1, invMat)
    circ      = mpatches.Circle(transPos1, ptRad/2, fc='g')
    ax.add_patch(circ)
    print 'in appendPtCURVE4 [%5.2f, %5.2f]' % (xdata1, ydata1)
    pathCode = [Path.CURVE4, [xdata1, ydata1]]
    isAnchor = True
    ptList.append( {'xPos'     : transPos1[0],
                    'yPos'     : transPos1[1],
                    'absPos'   : absPos,
                    'transPos' : transPos,
                    'selected' : True,
                    'circle'   : circ,
                    'pathcode' : pathCode,
                    'isanchor' : isAnchor,
                   })
                   
    transPos2 = [xdata2, ydata2, 1,]
    absPos2   = np.matmul(transPos2, invMat)
    circ      = mpatches.Circle(transPos2, ptRad/2, fc='g')
    ax.add_patch(circ)
    print 'in appendPtCURVE4 [%5.2f, %5.2f]' % (xdata2, ydata2)
    pathCode  = [Path.CURVE4, [xdata2, ydata2]]
    isAnchor  = True
    ptList.append( {'xPos'     : transPos2[0],
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
    circ      = mpatches.Circle(transPos3, ptRad, fc='b')
    ax.add_patch(circ)
    print 'in appendPtCURVE4 [%5.2f, %5.2f]' % (xdata3, ydata3)
    pathCode  = [Path.CURVE4, [xdata3, ydata3]]
    isAnchor  = False
    nPts += 1
    ptList.append( {'xPos'     : transPos3[0],
                    'yPos'     : transPos3[1],
                    'absPos'   : absPos3,
                    'transPos' : transPos3,
                    'selected' : True,
                    'circle'   : circ,
                    'pathcode' : pathCode,
                    'isanchor' : isAnchor,
                   })
    
############[ event handlers ]##############

def on_press(event):
    global selectedPt, buttonState, patchList
    global nPts, nAnchors
    print 'in on_press'
    if event.inaxes is not ax:
        return
    inAPoint = False
    for pt in ptList:
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
                                    [xdata,                ydata])
            elif nPts > 2:
                appendPtCURVE3APPEND([ptList[-1+nAnchors]['xPos'], ptList[-1+nAnchors]['yPos']],
                                    [xdata,                ydata])
                            
        elif curveType == 'curve4':
            if lnPts == 0:
                appendPtMOVETO([xdata, ydata])  
            elif nPts == 1:
                appendPtCURVE4TO([xdata, ydata],
                                 [ptList[0]['xPos'], ptList[0]['yPos']]) 
            elif nPts == 2:
                appendPtCURVE4CLOSE([xdata, ydata])
            elif nPts > 2:
                appendPtCURVE4APPEND([xdata, ydata])
                            
        generatePolyPath()  
        fig.canvas.draw()
        selectedPt = ptList[-1]    

def on_release(event):
    global buttonState, selectedPt
    print 'in on_release'
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

fig.canvas.mpl_connect('button_release_event', on_release)
fig.canvas.mpl_connect('button_press_event',   on_press)
fig.canvas.mpl_connect('motion_notify_event',  on_motion)

plt.show()

# Gef fig manager to raise window in top left corner (10,10)
figmgr=plt.get_current_fig_manager()
figmgr.canvas.manager.window.raise_()
geom=figmgr.window.geometry()
(xLoc,yLoc,dxWidth,dyHeight)=geom.getRect()
figmgr.window.setGeometry(10,10,dxWidth,dyHeight)
