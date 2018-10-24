""" Plotting routines"""
from __future__ import print_function, absolute_import, division, unicode_literals
import numpy as np
from matplotlib import pyplot as plt

###Display a white light image
def show_whitelight(whiteim):

    ydim, xdim = whiteim.shape

    plt.imshow(np.log10(np.abs(whiteim)), cmap=plt.get_cmap('viridis'))
    plt.xticks(np.arange(0, xdim, 5))
    plt.yticks(np.arange(0, ydim, 5))

    plt.grid(color='r', linestyle='-', linewidth=2, alpha=0.3)
    plt.gca().invert_yaxis()
    plt.show()

def show_narrowband(image):

    plt.imshow(nbimage, cmap=plt.get_cmap('viridis'))
    plt.colorbar()
    plt.show()
