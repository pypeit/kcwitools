""" Spectra routines"""
from __future__ import print_function, absolute_import, division, unicode_literals

import numpy as np
import pdb
import warnings

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

    spec = np.nansum(sub_flux, axis=(1,2))
    err_spec = np.sqrt(np.sum(sub_var, axis=(1,2)))

    # Replace NAN
    bad = np.isnan(err_spec)
    if np.any(bad):
        warnings.warn("You have NAN in your variance cube.  Be warned. Replacing with 0.")
        err_spec[bad] = 0.


    # Object
    xspec = XSpectrum1D.from_tuple((wave, spec, err_spec))

    if outfile is not None:
        xspec.write_to_fits(outfile)

    return xspec

def extract_rectangle(xcen, ycen, wave, flux, var, deltax=5,deltay=5, outfile=None):
    """

    Args:
        xcen: int
        ycen: int
        wave:
        flux:
        var:
        deltax: int, optional
          Size of the deltax to extract
        deltay: int, optional
          Size of the deltay t extract
        outfile: str, optional
          FITS file name

    Returns:
        xspec: XSpectrum1D

    """


    # create an error array assuming poisson noise for flux
    sub_flux = flux[:, ycen-deltay:ycen+deltay, xcen-deltax:xcen+deltax]
    sub_var = var[:, deltay:ycen+deltay, xcen-deltax:xcen+deltax]

    spec = np.nansum(sub_flux, axis=(1,2))
    err_spec = np.sqrt(np.sum(sub_var, axis=(1,2)))

    # Replace NAN
    bad = np.isnan(err_spec)
    if np.any(bad):
        warnings.warn("You have NAN in your variance cube.  Be warned. Replacing with 0.")
        err_spec[bad] = 0.


    # Object
    xspec = XSpectrum1D.from_tuple((wave, spec, err_spec))

    if outfile is not None:
        xspec.write_to_fits(outfile)

    return xspec

