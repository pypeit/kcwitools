""" I/O routines"""
from __future__ import print_function, absolute_import, division, unicode_literals

import numpy as np
import os
from astropy.io import fits

def open_kcwi_cube(infil):
    """
    Open KCWI cube and return a flux cube and the header (which contains the wavelengths)

    Args:
        infil:

    Returns:
        flux : ndarray
        hdu_hdr

    """
    # Open
    hdu = fits.open(infil)
    hdu_hdr = hdu['PRIMARY'].header

    # Flux
    flux = hdu['PRIMARY'].data

    # Return
    return hdu_hdr, flux


def write_cube(hdr, cube, outfile):
    """
    Write data cube to disk

    Args:
        hdr:
        cube:
        outfile:

    Returns:

    """
    hdu_coadd = fits.PrimaryHDU(cube, header=hdr)
    hdu_coadd.writeto(outfile, overwrite=True)
