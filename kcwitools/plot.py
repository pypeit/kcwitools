""" Plotting routines"""
from __future__ import print_function, absolute_import, division, unicode_literals
import numpy as np
from matplotlib import pyplot as plt
from astropy.wcs import WCS



###Display a white light image
def show_whitelight(whiteim, header=False):
    ydim, xdim = whiteim.shape
    fig = plt.figure()
    if header == False:
        ax = fig.add_subplot(111)
        d=ax.imshow(np.log10(np.abs(whiteim)), cmap=plt.get_cmap('viridis'), origin='lower')
    else:
        hdr = tweak_header(header)
        ax = fig.add_subplot(111, projection=WCS(hdr))
        d = ax.imshow(np.log10(np.abs(whiteim)), cmap=plt.get_cmap('viridis'), origin='lower')
        ax.set_xlabel('Right Ascension (J2000)')
        ax.set_ylabel('Declination (J2000)')
    plt.xticks(np.arange(0, xdim, 5))
    plt.yticks(np.arange(0, ydim, 5))
    plt.grid(color='gray', ls='dashed')
    fig.colorbar(d)
    plt.show()


def show_narrowband(nbimage, header=False):
    fig = plt.figure()
    if header == False:
        ax = fig.add_subplot(111)
        d=ax.imshow(nbimage, cmap=plt.get_cmap('viridis'), origin='lower')
    else:
        hdr = tweak_header(header)
        ax = fig.add_subplot(111, projection=WCS(hdr))
        d = ax.imshow(nbimage, cmap=plt.get_cmap('viridis'), origin='lower')
        plt.grid(color='gray', ls='dashed')
        ax.set_xlabel('Right Ascension (J2000)')
        ax.set_ylabel('Declination (J2000)')
    fig.colorbar(d)
    plt.show()


def tweak_header(header):
    '''
    Header tweaks to make the 3D kcwi header compatible with astropy image header object
    '''
    header['NAXIS']=2
    header['WCSDIM']=2
    if 'CRVAL3' in header.keys():
        header.remove('CRVAL3')
    if 'CRPIX3' in header.keys():
        header.remove('CRPIX3')
    if 'CD3_3' in header.keys():
        header.remove('CD3_3')
    if 'NAXIS3' in header.keys():
        header.remove('NAXIS3')
    if 'CUNIT3' in header.keys():
        header.remove('CUNIT3')
    if 'CTYPE3' in header.keys():
        header.remove('CTYPE3')
    if 'CNAME3' in header.keys():
        header.remove('CNAME3')

    return header
