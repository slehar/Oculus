"""
Demo of a PathPatch object.
"""
import matplotlib.path as mpath
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
from matplotlib.widgets import RadioButtons

#fig, ax = plt.subplots()
fig = plt.figure(figsize=(6, 6))
ax = fig.add_axes([.1, .1, .8, .8])


rax=plt.axes([.01,.1,.1,.5])
rax.set_xticks([])
rax.set_yticks([])

#this is list containing the data points, 3-6 pairs
pointlist = []
Radio=RadioButtons(rax, ('3','4','5','6'))
def howmanypoints(label):
    global pointlist
    print label
    if label == '3':
        print 'three points'
#        pointlist.append([[0,0],[2,-1],[1,2.8]])
        pointlist.append([0,  0])
        pointlist.append([2, -1])
        pointlist.append([1, 2.8])
        print pointlist
        Path = mpath.Path
        path_data = [
                (Path.MOVETO,    pointlist[0]),
                (Path.LINETO,    pointlist[1]),
                (Path.LINETO,    pointlist[2]),
                (Path.CLOSEPOLY, pointlist[0])]
        print 'path data:'
        print path_data
        codes, verts = zip(*path_data)
        path = mpath.Path(verts, codes)
        patch = mpatches.PathPatch(path, facecolor='r', alpha=0.5)
        ax.add_patch(patch)
        ax.set_xlim(-4,4)
        ax.set_ylim(-4,4)
        
Radio.on_clicked(howmanypoints)


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

def on_press(event):
    print 'button_press'
    
fig.canvas.mpl_connect('button_press_event', on_press)
plt.show()
