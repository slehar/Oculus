''' utils module

General utilities:

    axes = get_axes(fig, rect, title=None)

    slider = get_slider(fig, rect, label, valmin, valmax, valinit

'''

from   matplotlib.widgets import Slider

# Get 'axes' (rectangular regions) for each image in the figure
def get_axes(fig, rect, title=None):
    ax = fig.add_axes(rect)
    ax.axes.set_xticks([])
    ax.axes.set_yticks([])
    if title:
        ax.set_title(title)
    return ax

# Get slider for a given figure
def get_slider(fig, rect, label, valmin, valmax, valinit):
    ax = get_axes(fig, rect)
    return Slider(ax, label, valmin, valmax, valinit)

# Re-size to even dimensions
def resize_even_dim(imgNp):
    ySize, xSize = imgNp.shape
    if ySize % 2 == 0:
        print 'y dimension is even'
    else:
        print 'y dimension is odd, so drop a row'
        ySize = ySize - 1
        imgNp = imgNp[0:ySize]

    if xSize % 2 == 0:
        print 'x dimension is even'
    else:
        print 'x dimension is odd, so drop a column'
        xSize = xSize - 1
        imgNp = imgNp[:, :xSize]
    return imgNp
