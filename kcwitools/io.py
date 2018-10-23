""" I/O routines"""
from __future__ import print_function, absolute_import, division, unicode_literals

import numpy as np
import os
from astropy.io import fits

def open_kcwi_cube(infil):
    """
    Open KCWI cube and return a flux cube and wavelength array

    Args:
        infil:

    Returns:
        wave : ndarray
        flux : ndarray
        hdu_hdr

    """
    # Open
    hdu = fits.open(infil)
    hdu_hdr = hdu['PRIMARY'].header

    # Flux
    flux = hdu['PRIMARY'].data

    # Wavelengths
    crval3 = hdu_hdr['CRVAL3']
    crpix3 = hdu_hdr['CRPIX3']
    cd3_3 = hdu_hdr['CD3_3']
    wavedim, ydim, xdim = flux.shape
    wave = crval3 + (crpix3 + np.arange(0, wavedim, 1.0)) * cd3_3
    # Return
    return wave, flux, hdu_hdr



