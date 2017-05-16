"""
Demo of a PathPatch object.
"""
import numpy as np
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt


# global variables
ptList = []
inAPoint = False
selectedPt = None
scale = 1.
ptRad = 0.001
buttonState = False

#fig, ax = plt.subplots()
fig = plt.figure(figsize=(6, 6))
ax = fig.add_axes([.1, .1, .8, .8])

#Path = mpath.Path
#path_data = [
#    (Path.MOVETO, (1.58, -2.57)),
#    (Path.CURVE4, (0.35, -1.1)),
#    (Path.CURVE4, (-1.75, 2.0)),
#    (Path.CURVE4, (0.375, 2.0)),
#    (Path.LINETO, (0.85, 1.15)),
#    (Path.CURVE4, (2.2, 3.2)),
#    (Path.CURVE4, (3, 0.05)),
#    (Path.CURVE4, (2.0, -0.5)),
#    (Path.CLOSEPOLY, (1.58, -2.57)),
#    ]
#codes, verts = zip(*path_data)
#path = mpath.Path(verts, codes)
#patch = mpatches.PathPatch(path, facecolor='r', alpha=0.5)
#ax.add_patch(patch)

# plot control points and connecting lines
#x, y = zip(*path.vertices)
#line, = ax.plot(x, y, 'go-')

ax.grid()
ax.axis('equal')

# callback functions, handling events


transMat = np.array([[scale,   0,     0],
                     [0,       scale, 0],
                     [0,       0,     1]])

invMat = np.linalg.inv(transMat)


def on_press(event):
    global selectedPt, buttonState
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
            fig.canvas.draw()
            break
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
figmgr.canvas.manager.window.raise_()
geom=figmgr.window.geometry()
(xLoc,yLoc,dxWidth,dyHeight)=geom.getRect()
figmgr.window.setGeometry(10,10,dxWidth,dyHeight)
