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


def write_image(hdr, image, outfile):
    """
    Write image to disk, e.g. whitelight

    Args:
        hdr:
        cube:
        outfile:

    Returns:

    """
    # Fuss with header
    header = hdr.copy()
    header.remove('CRVAL3')
    header.remove('CRPIX3')
    header.remove('CD3_3')
    header.remove('CDELT3')
    header.remove('NAXIS3')
    header['NAXIS']=2
    #

    hdu_coadd = fits.PrimaryHDU(image, header=header)
    hdu_coadd.writeto(outfile, overwrite=True)
    print("Wrote image to {}".format(outfile))
