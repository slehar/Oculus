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
ptRad = 0.02
buttonState = False
path_data = []
patchList = []

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
    global ptList, data, mute
    if event.key == 'q':
        plt.close()        
fig.canvas.mpl_connect('key_press_event',       keypress)


# callback functions, handling events


transMat = np.array([[scale,   0,     0],
                     [0,       scale, 0],
                     [0,       0,     1]])

invMat = np.linalg.inv(transMat)

def extendPolyPath(xdata, ydata): 
    global path_data, patchList
    print '\nin extendPolyPath [%5.2f, %5.2f] len=%d'%(xdata,ydata,len(path_data))

    nPts = len(ptList)
    if nPts == 1:             # If first point just MOVETO
        path_data.append([Path.MOVETO, [xdata, ydata]])
        print '  Path.MOVETO (%5.2f, %5.2f)'%(xdata, ydata)
           
    elif nPts == 2:           # If second point draw LINETO
        path_data.append([Path.LINETO, [xdata, ydata]])
        print '  Path.LINETO (%5.2f, %5.2f)'%(xdata, ydata)
            
    elif nPts == 3:           # If third point draw LINETO and CLOSEPOLY
        path_data.append([Path.LINETO, [xdata, ydata]])
        print '  Path.LINETO (%5.2f, %5.2f)'%(xdata, ydata)
        path_data.append([Path.CLOSEPOLY,  [ptList[0]['xPos'],ptList[0]['yPos']]])
        print '  Path.CLOSEPOLY (%5.2f %5.2f)'%(ptList[0]['xPos'],ptList[0]['yPos'])
        for pa in patchList:    # Remove previous polygon before adding new path
            pa.remove()
  
    elif nPts >= 4:           # All further points insert LINETO before CLOSEPOLY        
        path_data[-1] = [Path.LINETO, [xdata, ydata]]
        print '  Path.LINETO (%5.2f, %5.2f)'%(xdata, ydata)
        path_data.append([Path.CLOSEPOLY,  [ptList[0]['xPos'],ptList[0]['yPos']]])
        print '  Path.CLOSEPOLY (%5.2f %5.2f)'%(ptList[0]['xPos'],ptList[0]['yPos'])
        for pa in patchList:    # Remove previous polygon before adding new path
            pa.remove()

    patchList = []
    codes, verts = zip(*path_data)
    path = Path(verts, codes)
    patch = mpatches.PathPatch(path, facecolor='r', alpha=0.5)
#    print '\n%r'%patch
    patchList.append(patch)
    ax.add_patch(patch)
#    print '\n%r'%path_data




def editPolyPath(xdata, ydata):
    global ptList, path_data, patchList
    print '\nin editPolyPath [%5.2f, %5.2f] path_len=%d'%(xdata,ydata,len(path_data))

    nPts = len(ptList)
    print 'nPts = %d' % nPts
    for pt in ptList:
        if pt['selected']:
            print '  SELECTED'
            if nPts == 0:             # If no points just return
                print '  nPts = %d'%nPts
                return
                   
            elif nPts == 1:           # If first point replace with MOVETO
                print '  nPts = %d'%nPts
                path_data[0] = [Path.MOVETO, [xdata, ydata]]
                print '  Path.MOVETO (%5.2f, %5.2f)'%(xdata, ydata)
                    
            elif nPts == 2:           # If third point draw LINETO and CLOSEPOLY
                print '  nPts = %d'%nPts
                for ix, pt in enumerate(ptList):
                    print '  ix = %d sel=%r '%(ix, pt['selected'])
                    if pt['selected']:
                        print 'XXXXXXXXXXXXXXXXX!!!'
                        path_data[ix][1] = [xdata, ydata]
                        print '  ix %d  Path.MOVETO (%5.2f, %5.2f)'%(ix, xdata, ydata)
          
            elif nPts >= 3:           # All further points insert LINETO before CLOSEPOLY        
                print '  nPts = %d'%nPts
                for ix, pt in enumerate(ptList):
                    if pt['selected']:
                        path_data[ix][1] = [xdata, ydata]
                        print 'ix %d  Path.MOVETO (%5.2f, %5.2f)'%(ix, xdata, ydata)
                for pa in patchList:    # Remove previous polygon before adding new path
                    pa.remove()
          
        
#        path_data.append((Path.CLOSEPOLY,  (ptList[0]['xPos'],ptList[0]['yPos'])))
#        print '  Path.CLOSEPOLY (%5.2f %5.2f)'%(ptList[0]['xPos'],ptList[0]['yPos'])

    for pa in patchList:    # Remove previous polygon before adding new path
        pa.remove()

    patchList = []
    codes, verts = zip(*path_data)
    path = Path(verts, codes)
    patch = mpatches.PathPatch(path, facecolor='r', alpha=0.5)
#    print '\n%r'%patch
    patchList.append(patch)
    ax.add_patch(patch)
#    print '\n%r'%path_data


'''
def updatePolyPath(xdata, ydata, replace=False): # Update polygon when points move
    global path_data, patchList

    print '\nin updatePolyPath [%5.2f, %5.2f] replace=%r, len=%d'%(xdata,ydata,replace, len(path_data))

#    path_len = len(path_data) 
    nPts = len(ptList)
    if nPts == 0:             # If first point just MOVETO
#        print 'append len=0'
        path_data.append((Path.MOVETO, (xdata, ydata)))
        print '  Path.MOVETO (%5.2f, %5.2f)'%(xdata, ydata)
           
    elif nPts == 1:           # If second point draw LINETO
        if replace:
#            print 'replace len=1'
            path_data[:] = []
            path_data = [(Path.MOVETO, (xdata, ydata))]
            print '  Path.MOVETO (%5.2f, %5.2f)'%(xdata, ydata)
        else:
#            print 'append len=1'
            path_data.append((Path.LINETO, (xdata, ydata)))
            print '  Path.LINETO (%5.2f, %5.2f)'%(xdata, ydata)
            
    elif nPts == 2:           # If third point draw LINETO and CLOSEPOLY
        if replace:
#            print 'replace len=2'
            path_data[:] = []
            
            
            if ptList[0] is selectedPt:
                print '  Selected Pt: 0'
                path_data.append((Path.MOVETO,     (xdata,             ydata)))
                path_data.append((Path.LINETO,     (ptList[1]['xPos'], ptList[1]['yPos'])))
                print '  Path.MOVETO (%5.2f, %5.2f)'%(xdata, ydata)
                print '  Path.LINETO (%5.2f %5.2f)'%(xdata,ydata)
            elif ptList[1] is selectedPt:
                print '  Selected Pt: 1'
                path_data.append((Path.MOVETO,     (ptList[0]['xPos'], ptList[0]['yPos'])))
                path_data.append((Path.LINETO,     (xdata,             ydata)))
                print '  Path.MOVETO (%5.2f, %5.2f)'%(ptList[0]['xPos'], ptList[0]['yPos'])
                print '  Path.LINETO (%5.2f %5.2f)'%(xdata,ydata)
                        
        else:
#            print 'append len=2'
            path_data.append((Path.LINETO, (xdata, ydata)))
            path_data.append((Path.CLOSEPOLY, (ptList[0]['xPos'],ptList[0]['yPos'])))
            print '  Path.LINETO (%5.2f, %5.2f)'%(xdata, ydata)
            print '  Path.CLOSEPOLY (%5.2f, %5.2f)'%(ptList[0]['xPos'],ptList[0]['yPos'])
  
    elif nPts >= 3:           # All further points insert LINETO before CLOSEPOLY
        if replace:
#            print 'replace len=>3'
            path_data[:] = []
            
            
            for vertex in range(path_len - 1):
                print 'vertex : %d'%vertex
                if ptList[vertex] is selectedPt:
                    print '  Selected Pt: %d'%vertex
                    if vertex == 0:
                        path_data.append((Path.MOVETO,     (xdata,             ydata)))
                        print '  Path.MOVETO (%5.2f %5.2f)'%(xdata,ydata)
                    else:
                        path_data.append((Path.LINETO,      (ptList[vertex]['xPos'], ptList[vertex]['yPos'])))
                        print '  Path.LINETO (%5.2f %5.2f)'%(ptList[vertex]['xPos'], ptList[vertex]['yPos'])
                else:
                    if vertex == 0:
                        path_data.append((Path.MOVETO,     (xdata,             ydata)))
                        print '  Path.MOVETO (%5.2f %5.2f)'%(xdata,ydata)
                    else:
                        path_data.append((Path.LINETO,      (ptList[vertex]['xPos'], ptList[vertex]['yPos'])))
                        print '  Path.LINETO (%5.2f %5.2f)'%(ptList[vertex]['xPos'], ptList[vertex]['yPos'])
            
            path_data.append((Path.CLOSEPOLY,  (ptList[0]['xPos'],ptList[0]['yPos'])))
            print '  Path.CLOSEPOLY (%5.2f %5.2f)'%(ptList[0]['xPos'],ptList[0]['yPos'])
        else:
#            print 'append len>=3'
            path_data = path_data[:-1]
            path_data.append((Path.LINETO, (xdata, ydata)))
            path_data.append((Path.CLOSEPOLY, (ptList[0]['xPos'],ptList[0]['yPos'])))
            print '  Path.LINETO (%5.2f %5.2f))'%(xdata, ydata)
            print '  Path.CLOSEPOLY (%5.2f %5.2f)'%(ptList[0]['xPos'],ptList[0]['yPos'])



    for pa in patchList:    # Remove previous polygon before adding new path
        pa.remove()
#            print '\n%r'%pa
    patchList = []
    codes, verts = zip(*path_data)
    path = Path(verts, codes)
    patch = mpatches.PathPatch(path, facecolor='r', alpha=0.5)
#    print '\n%r'%patch
    patchList.append(patch)
    ax.add_patch(patch)
#    print '\n%r'%path_data
'''


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
#        updatePolyPath(xdata, ydata, replace=False)
        nPts = len(ptList)
        if nPts == 0:
            polycode = ''
        elif nPts == 1:
            polycode = Path.MOVETO
        elif nPts >= 2:
            polycode = Path.LINETO
        ptList.append({'xPos'     : transPos[0],
                       'yPos'     : transPos[1],
                       'absPos'   : absPos,
                       'transPos' : transPos,
                       'selected' : True,
                       'circle'   : circ,
                       'polycode' : polycode,
                       })

        extendPolyPath(xdata, ydata)           ###############<<<<

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
#        updatePolyPath(xdata, ydata, replace=True)
        editPolyPath(xdata, ydata)                  ######<<<<<<<<
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
