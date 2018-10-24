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

def write_spec():
    prihdu = fits.PrimaryHDU(spec)
    hdu_out = fits.HDUList([prihdu])
    prihdu.name = 'FLUX'

    sighdu = fits.ImageHDU(err_spec)
    sighdu.name = 'ERROR'
    hdu_out.append(sighdu)

    wvhdu = fits.ImageHDU(wave)
    wvhdu.name = 'WAVELENGTH'
    hdu_out.append(wvhdu)

    # WATCH OUT:  purposely flipping xcen and ycen here to name things according to the image diplayed
    if filebase != False:
        hdu_out.writeto(filebase + '_' + np.str(ycen) + '_' + np.str(xcen) + '_5x5.fits', overwrite=True)

