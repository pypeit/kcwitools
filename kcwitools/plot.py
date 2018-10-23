""" Plotting routines"""
from __future__ import print_function, absolute_import, division, unicode_literals

###Display a white light image
def show_whitelight(infil):
    wave, flux, hdr = open_kcwi_cube(infil)
    minwave = 3600.
    maxwave = 5500.
    wavedim, ydim, xdim = flux.shape
    whiteim = np.zeros((ydim, xdim))
    slices = np.where((wave >= minwave) & (wave <= maxwave))
    slices = slices[0]

    for i in range(slices.size):
        whiteim[:, :] += flux[slices[i], :, :]

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
