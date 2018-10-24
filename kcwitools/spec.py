""" Spectra routines"""
from __future__ import print_function, absolute_import, division, unicode_literals

import numpy as np
import pdb

from linetools.spectra.xspectrum1d import XSpectrum1D

def extract_square(xcen, ycen, wave, flux, var, squaresize=5, outfile=None):
    """

    Args:
        xcen: int
        ycen: int
        wave:
        flux:
        var:
        squaresize: int, optional
          Size of the square to extract
          Should be an odd number
        outfile: str, optional
          FITS file name

    Returns:
        xspec: XSpectrum1D

    """
    #wave, flux, hdr = open_kcwi_cube(fil)
    #wave, var, hdr = open_kcwi_cube(varfil)

    halfbox = (squaresize-1)//2

    # create an error array assuming poisson noise for flux
    sub_flux = flux[:, ycen-halfbox:ycen+halfbox+1, xcen-halfbox:xcen+halfbox+1]
    sub_var = var[:, ycen-halfbox:ycen+halfbox+1, xcen-halfbox:xcen+halfbox+1]

    spec = np.sum(sub_flux, axis=(1,2))
    err_spec = np.sqrt(np.sum(sub_var, axis=(1,2)))

    # Object
    xspec = XSpectrum1D.from_tuple((wave, spec, err_spec))

    if outfile is not None:
        xspec.write_to_fits(outfile)

    return xspec


def extract_1D_spec(xcen, ycen, fil, varfil, delx, dely):
    wave, flux, hdr = open_kcwi_cube(fil)
    wave, var, hdr = open_kcwi_cube(varfil)

    # assumes sky subtracted data, and returns sqrt(var) var
    wavedim, ydim, xdim = flux.shape

    ###WATCH OUT FOR THIS:  NOW SWITCHING X WITH Y BECAUSE OF HOW CUBES READ IN and displayed
    # in the white light image

    xtmp = ycen
    ycen = xcen
    xcen = xtmp

    delxtmp = dely
    dely = delx
    delx = delxtmp

    # Intitiate Spectrum
    spec = np.zeros(wavedim)
    err_spec = np.zeros(wavedim)

    # create an error array assuming poisson noise for flux

    for i in range(wavedim):
        slice1 = flux[i, xcen - delx:xcen + delx, ycen - dely:ycen + dely]
        varslice = var[i, xcen - delx:xcen + delx, ycen - dely:ycen + dely]
        spec[i] = np.sum(slice1)
        err_spec[i] = np.sqrt(np.abs(np.sum(varslice)))

    return spec, err_spec
