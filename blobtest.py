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
path_data = []
patchList = []
inAPoint = False
selectedPt = None
scale = 1.
ptRad = 0.02
buttonState = False
altKeyPressed = False

# Figure and axes
plt.close('all')
fig = plt.figure(figsize=(6, 6))
ax = fig.add_axes([.1, .1, .8, .8])
ax.set_xlim([-1., 1.])
ax.set_ylim([-1., 1.])
ax.grid()
#ax.axis('equal')

# Keypress 'q' to quit
def keypress(event):
    global ptList, data, mute, altKeyPressed
        
    if event.key == 'alt':
        altKeyPressed = True
        print 'alt key pressed'
    if event.key == 'q':
        plt.close()        
fig.canvas.mpl_connect('key_press_event', keypress)

# Key release
def keyrelease(event):
    global ptList, data, mute, altKeyPressed
        
    if event.key == 'alt':
        altKeyPressed = False
        print 'alt key UN-pressed'
        
fig.canvas.mpl_connect('key_release_event', keyrelease)

# callback functions, handling events


transMat = np.array([[scale,   0,     0],
                     [0,       scale, 0],
                     [0,       0,     1]])

invMat = np.linalg.inv(transMat)

def generatePolyPath(): 
    global path_data, patchList
    print 'in generatePolyPath'
    
    path_data = []
    for pt in ptList:
        path_data.append(pt['pathcode'])
    print path_data
    
    for pa in patchList:    # Remove previous polygon before adding new path
        pa.remove()
    patchList = []
    codes, verts = zip(*path_data)
    path = Path(verts, codes)
    patch = mpatches.PathPatch(path, facecolor='r', alpha=0.5)
    patchList.append(patch)
    ax.add_patch(patch)
    

def extendPolyPath(xdata, ydata): 
    global path_data, patchList
#    print '\nin extendPolyPath [%5.2f, %5.2f] len=%d'%(xdata,ydata,len(path_data))

    nPts = len(ptList)
    if nPts == 1:             # If first point just MOVETO
        path_data.append([Path.MOVETO, [xdata, ydata]])
#        print '  Path.MOVETO (%5.2f, %5.2f)'%(xdata, ydata)
           
    elif nPts == 2:           # If second point draw LINETO
        path_data.append([Path.LINETO, [xdata, ydata]])
#        print '  Path.LINETO (%5.2f, %5.2f)'%(xdata, ydata)
            
    elif nPts == 3:           # If third point draw LINETO and CLOSEPOLY
        path_data.append([Path.LINETO, [xdata, ydata]])
#        print '  Path.LINETO (%5.2f, %5.2f)'%(xdata, ydata)
        path_data.append([Path.CLOSEPOLY,  [ptList[0]['xPos'],ptList[0]['yPos']]])
#        print '  Path.CLOSEPOLY (%5.2f %5.2f)'%(ptList[0]['xPos'],ptList[0]['yPos'])
        for pa in patchList:    # Remove previous polygon before adding new path
            pa.remove()
  
    elif nPts >= 4:           # All further points insert LINETO before CLOSEPOLY        
        path_data[-1] = [Path.LINETO, [xdata, ydata]]
#        print '  Path.LINETO (%5.2f, %5.2f)'%(xdata, ydata)
        path_data.append([Path.CLOSEPOLY,  [ptList[0]['xPos'],ptList[0]['yPos']]])
#        print '  Path.CLOSEPOLY (%5.2f %5.2f)'%(ptList[0]['xPos'],ptList[0]['yPos'])
        for pa in patchList:    # Remove previous polygon before adding new path
            pa.remove()

    patchList = []
    codes, verts = zip(*path_data)
    path = Path(verts, codes)
    patch = mpatches.PathPatch(path, facecolor='r', alpha=0.5)
    patchList.append(patch)
    ax.add_patch(patch)


def editPolyPath(xdata, ydata):
    global ptList, path_data, patchList
#    print '\nin editPolyPath [%5.2f, %5.2f] path_len=%d'%(xdata,ydata,len(path_data))

    nPts = len(ptList)
#    print 'nPts = %d' % nPts
    for pt in ptList:
        if pt['selected']:
#            print '  SELECTED'
            if nPts == 0:             # If no points just return
#                print '  nPts = %d'%nPts
                return
                   
            elif nPts == 1:           # If first point
#                print '  nPts = %d'%nPts
                if altKeyPressed:
                    [oldX, oldY] = path_data[0][1]
                    path_data = []
                    path_data.append([Path.MOVETO, [oldX, oldY]])
                    path_data.append([Path.CURVE3, [xdata, ydata]])
                    path_data.append([Path.CURVE3, [0.,    0.]])
#                    path_data.append([Path.CURVE3, [xdata, ydata]])
#                    print '  Path.MOVETO (%5.2f, %5.2f)'%(oldX, oldY)
#                    print '  Path.CURVE3 (%5.2f, %5.2f)'%(xdata, ydata)
                    print path_data
#                    for pa in patchList:    # Remove previous polygon before adding new path
#                        pa.remove()

                else:
                    path_data[0] = [Path.MOVETO, [xdata, ydata]]
                    print '  Path.MOVETO (%5.2f, %5.2f)'%(xdata, ydata)
                    
            elif nPts == 2:           # If third point draw LINETO and CLOSEPOLY
#                print '  nPts = %d'%nPts
                for ix, pt in enumerate(ptList):
#                    print '  ix = %d sel=%r '%(ix, pt['selected'])
                    if pt['selected']:
                        path_data[ix][1] = [xdata, ydata]
#                        print '  ix %d  Path.MOVETO (%5.2f, %5.2f)'%(ix, xdata, ydata)
          
            elif nPts >= 3:           # All further points insert LINETO before CLOSEPOLY        
#                print '  nPts = %d'%nPts
                for ix, pt in enumerate(ptList):
                    if pt['selected']:
                        path_data.pop()
                        path_data[ix][1] = [xdata, ydata]
#                        print '  ix %d  Path.MOVETO (%5.2f, %5.2f)'%(ix, xdata, ydata)
                        path_data.append((Path.CLOSEPOLY,  (ptList[0]['xPos'],ptList[0]['yPos'])))
#                        print '  Path.CLOSEPOLY (%5.2f %5.2f)'%(ptList[0]['xPos'],ptList[0]['yPos'])

    for pa in patchList:    # Remove previous polygon before adding new path
        pa.remove()

    patchList = []
    codes, verts = zip(*path_data)
    path = Path(verts, codes)
    patch = mpatches.PathPatch(path, facecolor='r', alpha=0.5)
    patchList.append(patch)
    ax.add_patch(patch)

def on_press(event):
    global selectedPt, buttonState, path_data, patchList
    if event.inaxes is not ax:
        return
    inAPoint = False
    for pt in ptList:
        contains, attrd = pt['circle'].contains(event)
        if contains:
            inAPoint = True
            if pt['selected']:
                pt['selected'] = False
                pt['circle'].set_fc('blue')
            else:
                pt['selected'] = True
                selectedPt = pt
                pt['circle'].set_fc('red')
            break
    fig.canvas.draw()
    buttonState = True
        
    if not inAPoint:
        xdata = event.xdata
        ydata = event.ydata
        transPos = [xdata, ydata, 1,]
        absPos = np.matmul(transPos, invMat)
        circ = mpatches.Circle(transPos, ptRad)
        ax.add_patch(circ)
        if len(ptList) == 0:
            pathCode = [Path.MOVETO, [xdata, ydata]]
            ptList.append({'xPos'     : transPos[0],
                           'yPos'     : transPos[1],
                           'absPos'   : absPos,
                           'transPos' : transPos,
                           'selected' : True,
                           'circle'   : circ,
                           'pathcode' : pathCode,
                           })
    
        elif len(ptList) == 1:
            pathCode = [Path.LINETO, [xdata, ydata]]
            ptList.append({'xPos'     : transPos[0],
                           'yPos'     : transPos[1],
                           'absPos'   : absPos,
                           'transPos' : transPos,
                           'selected' : True,
                           'circle'   : circ,
                           'pathcode' : pathCode,
                           })
        elif len(ptList) == 2:
            pathCode = [Path.LINETO, [xdata, ydata]]
            ptList.append({'xPos'     : transPos[0],
                           'yPos'     : transPos[1],
                           'absPos'   : absPos,
                           'transPos' : transPos,
                           'selected' : True,
                           'circle'   : circ,
                           'pathcode' : pathCode,
                           })
            pathCode = [Path.CLOSEPOLY, path_data[0][1]]
            ptList.append({'xPos'     : transPos[0],
                           'yPos'     : transPos[1],
                           'absPos'   : absPos,
                           'transPos' : transPos,
                           'selected' : False,
                           'circle'   : circ,
                           'pathcode' : pathCode,
                           })
        generatePolyPath()  
#        extendPolyPath(xdata, ydata)

        fig.canvas.draw()
        selectedPt = ptList[-1]
    

def on_release(event):
    global buttonState, selectedPt
    for pt in ptList:
        if pt['selected']:
            xdata = event.xdata
            ydata = event.ydata
            transPos = [xdata, ydata, 1,]
            pt['absPos'] = np.matmul(transPos, invMat)
            pt['circle'].center = transPos[:2]
            buttonState = False
            pt['selected'] = False
            selectedPt = None
            pt['circle'].set_fc('blue')
            fig.canvas.draw()
    buttonState = False

    
def on_motion(event):
    global xdata, ydata, selectedPt, ptList
    
    xdata = event.xdata
    ydata = event.ydata
        

    if buttonState:
        
        if altKeyPressed:
            print '*** ALT MOTION ***'
            selectedPt['anchor1'] = (xdata, ydata)
            print 'anchor1 = (%5.2f, %5.2f)' % (xdata, ydata)
            editPolyPath(xdata, ydata)
        else:
            
            absPos = [xdata, ydata, 1]
            transPos = np.matmul(absPos, transMat)
            selectedPt['circle'].center = transPos[:2]
            selectedPt['xPos'] = transPos[0]
            selectedPt['yPos'] = transPos[1]
            selectedPt['pathcode'][1] = [xdata, ydata]
            generatePolyPath()
#            editPolyPath(xdata, ydata)
            fig.canvas.draw()


fig.canvas.mpl_connect('button_release_event', on_release)
fig.canvas.mpl_connect('button_press_event',   on_press)
fig.canvas.mpl_connect('motion_notify_event',  on_motion)

plt.show()

# Gef fig manager to raise window in top left corner (10,10)
figmgr=plt.get_current_fig_manager()
#figmgr.canvas.manager.window.raise_()
#geom=figmgr.window.geometry()
#(xLoc,yLoc,dxWidth,dyHeight)=geom.getRect()
#figmgr.window.setGeometry(10,10,dxWidth,dyHeight)
