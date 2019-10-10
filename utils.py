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
