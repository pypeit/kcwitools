""" Image generation and manipulation routines"""
from __future__ import print_function, absolute_import, division, unicode_literals
import numpy as np
import pdb

from kcwitools import io
from kcwitools import utils
from astropy.io import fits

###trim off the crap parts of the KCWI cube medium slicer
###This is *aggressive* in the y direction because of dispersion
def kcwi_cube_trim_Medium(infil):
    hdr, flux = io.open_kcwi_cube(infil)
    trimflux = flux[:, 16:80, 6:28]  # note that python reads arrays weird.  Trim *down in y* then x.

    a, b = infil.split(".fits")
    outfil = a + '_trimmed.fits'

    hdu_out = fits.PrimaryHDU(trimflux, header=hdr)
    hdu_out.writeto(outfil, overwrite=True)


###trim off crap parts of the KCWI large slicer
def kcwi_cube_trim_Large(infil):
    hdr, flux = io.open_kcwi_cube(infil)
    trimflux = flux[:, 15:81, 2:26]  # note that python reads arrays weird.  Trim *down in y* then x.

    a, b = infil.split(".fits")
    outfil = a + '_trimmed.fits'

    hdu_out = fits.PrimaryHDU(trimflux, header=hdr)
    hdu_out.writeto(outfil, overwrite=True)


def build_whitelight(hdr, flux, minwave=3600., maxwave=5500., outfile=None):
    """
    Generate a white light image over a range of wavelengths

    Args:
        hdr: Header
        flux: ndarray
        minwave: float, optional
        maxwave: float, optional
        outfile: str, optional
          If provided, write image to disk

    Returns:
        whiteim: ndarray

    """
    # Wavelength
    wave = utils.build_wave(hdr)
    # Slices
    slices = np.where((wave >= minwave) & (wave <= maxwave))[0]
    # Do it
    whiteim = np.sum(flux[slices,:,:], axis=0)
    #for i in range(slices.size):
    #    whiteim[:, :] += flux[slices[i], :, :]

    if outfile is not None:
        io.write_image(hdr, whiteim, outfile)

    # Return
    return whiteim


def build_narrowband(hdr, flux, line, z=None, del_wave=2.0, sub_offimage=False, outfile=None):
    """
    Generate a narrow band image at a given wavelength

    Args:
        hdr:
        flux:
        line: float
          Wavelength of interest
        z: float, optional
          Convert to rest wavelength based on this redshift
        del_wave: float, optional
          Width of NB in rest-frame
          Default = 2 Ang
        sub_offimage: bool, optional
          Subtract off an off-band image?
        outfile: str, optional

    Returns:
        nbimage: ndarray

    """
    # Rest wave
    if z is None:
        z = 0.
    # Wavelengths
    wave = utils.build_wave(hdr)
    wrest = wave/(1+z)

    # Work out the slices for the on-band image:
    slices = np.abs(wrest - line) < del_wave

    # define windows on either side of the slice (watch out for a,b and <= vs <)
    offimage_width = 2. * del_wave

    slice_low = np.where((wrest >= (line - del_wave - offimage_width)) &
                         (wrest < (line - del_wave)))[0]
    slice_high = np.where((wrest >= (line + del_wave)) &
                          (wrest < (line + del_wave + offimage_width)))[0]

    # Make the on-band and two off-band images
    nbimage = np.average(flux[slices, :, :], 0)
    high = np.average(flux[slice_high, :, :], 0)
    low = np.average(flux[slice_low, :, :], 0)

    # Average the two off-band images and subtract.
    if sub_offimage:
        offimage = (low + high) / 2.0
        nbimage -= offimage

    # Output?
    if outfile is not None:
        io.write_image(hdr, nbimage, outfile)

    # Return
    return nbimage

def cube_skysub(fil, skyy1, skyy2, skyx1, skyx2, outfil):
    wave, flux, hdr = open_kcwi_cube(fil)
    wavedim, ydim, xdim = flux.shape
    sky = np.zeros(wavedim)

    for i in range(wavedim):
        sky[i] = np.median(flux[i, skyy1:skyy2, skyx1:skyx2])

    # Now Sky subtract the whole cube
    fluxmsky = np.zeros((wavedim, ydim, xdim))
    for i in range(wavedim):
        fluxmsky[i, :, :] = flux[i, :, :] - sky[i]

    hdu_out = fits.PrimaryHDU(fluxmsky, header=hdr)
    hdu_out.writeto(outfil, overwrite=True)


def var_skysub(varfil, skyy1, skyy2, skyx1, skyx2, outfil):
    wave, var, hdr = open_kcwi_cube(varfil)
    wavedim, ydim, xdim = var.shape
    sky = np.zeros(wavedim)

    for i in range(wavedim):
        sky[i] = np.median(var[i, skyy1:skyy2, skyx1:skyx2])

    # Now Sky subtract the whole cube
    varwsky = np.zeros((wavedim, ydim, xdim))
    for i in range(wavedim):
        varwsky[i, :, :] = var[i, :, :] + sky[i]

    hdu_out = fits.PrimaryHDU(varwsky, header=hdr)
    hdu_out.writeto(outfil, overwrite=True)

