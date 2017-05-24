"""
Demo of a PathPatch object.
See 
https://matplotlib.org/examples/shapes_and_collections/path_patch_demo.html
"""
import numpy as np
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
from matplotlib.path import Path


# global variables
ptList = []
inAPoint = False
selectedPt = None
scale = 1.
ptRad = 0.001
buttonState = False
path_data = []
patchList = []

# Figure and axes
plt.close('all')
fig = plt.figure(figsize=(6, 6))
ax = fig.add_axes([.1, .1, .8, .8])



ax.grid()
ax.axis('equal')

# Keypress 'q' to quit
def keypress(event):
    global ptList, data, mute
    if event.key == 'q':
        plt.close()        
fig.canvas.mpl_connect('key_press_event',       keypress)


# callback functions, handling events


transMat = np.array([[scale,   0,     0],
                     [0,       scale, 0],
                     [0,       0,     1]])

invMat = np.linalg.inv(transMat)

def updatePolyPath(xdata, ydata): # Update polygon when points move
    global path_data, patchList
    
    if len(path_data) == 0:             # If first point just MOVETO
        path_data.append((Path.MOVETO, (xdata, ydata)))
    elif len(path_data) == 1:           # If second point draw LINETO 
        path_data.append((Path.LINETO, (xdata, ydata)))
    elif len(path_data) == 2:           # If third point draw LINETO and CLOSEPOLY
        path_data.append((Path.LINETO, (xdata, ydata)))
        path_data.append((Path.CLOSEPOLY, (ptList[0]['xPos'],ptList[0]['yPos'])))
    elif len(path_data) >= 3:           # All further points insert LINETO before CLOSEPOLY
        path_data = path_data[:-1]
        path_data.append((Path.LINETO, (xdata, ydata)))
        path_data.append((Path.CLOSEPOLY, (ptList[0]['xPos'],ptList[0]['yPos'])))
        for pa in patchList:    # Remove previous polygon before adding new point
            pa.remove()
            print '\n%r'%pa
        patchList = []
      
    codes, verts = zip(*path_data)
    path = Path(verts, codes)
    patch = mpatches.PathPatch(path, facecolor='r', alpha=0.5)
    print '\n%r'%patch
    patchList.append(patch)
    ax.add_patch(patch)
    print '\n%r'%path_data


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
        ptList.append({'xPos':transPos[0],
                       'yPos': transPos[1],
                       'absPos':absPos,
                       'transPos':transPos,
                       'selected':True,
                       'circle':circ})

    
#        if len(path_data) == 0:             # If first point just MOVETO
#            path_data.append((Path.MOVETO, (xdata, ydata)))
#        elif len(path_data) == 1:           # If second point draw LINETO 
#            path_data.append((Path.LINETO, (xdata, ydata)))
#        elif len(path_data) == 2:           # If third point draw LINETO and CLOSEPOLY
#            path_data.append((Path.LINETO, (xdata, ydata)))
#            path_data.append((Path.CLOSEPOLY, (ptList[0]['xPos'],ptList[0]['yPos'])))
#        elif len(path_data) >= 3:           # All further points insert LINETO before CLOSEPOLY
#            path_data = path_data[:-1]
#            path_data.append((Path.LINETO, (xdata, ydata)))
#            path_data.append((Path.CLOSEPOLY, (ptList[0]['xPos'],ptList[0]['yPos'])))
#            for pa in patchList:    # Remove previous polygon before adding new point
#                pa.remove()
#                print '\n%r'%pa
#            patchList = []
#  
#        codes, verts = zip(*path_data)
#        path = Path(verts, codes)
#        patch = mpatches.PathPatch(path, facecolor='r', alpha=0.5)
#        print '\n%r'%patch
#        patchList.append(patch)
#        ax.add_patch(patch)
#        print '\n%r'%path_data

        updatePolyPath(xdata, ydata)

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
    
    if buttonState:
        xdata = event.xdata
        ydata = event.ydata
        absPos = [xdata, ydata, 1]
        transPos = np.matmul(absPos, transMat)
        selectedPt['circle'].center = transPos[:2]
        selectedPt['xPos'] = transPos[0]
        selectedPt['yPos'] = transPos[1]
        fig.canvas.draw()
        
fig.canvas.mpl_connect('button_release_event',on_release)
fig.canvas.mpl_connect('button_press_event',  on_press)
fig.canvas.mpl_connect('motion_notify_event', on_motion)

plt.show()

# Gef fig manager to raise window in top left corner (10,10)
figmgr=plt.get_current_fig_manager()
#figmgr.canvas.manager.window.raise_()
#geom=figmgr.window.geometry()
#(xLoc,yLoc,dxWidth,dyHeight)=geom.getRect()
#figmgr.window.setGeometry(10,10,dxWidth,dyHeight)
