""" Image generation and manipulation routines"""
from __future__ import print_function, absolute_import, division, unicode_literals
import numpy as np
from astropy.wcs import WCS
from astropy.wcs.utils import proj_plane_pixel_scales
from kcwitools import plot as kp
from copy import deepcopy
import pdb
from scipy.stats import sigmaclip

from kcwitools import io
from kcwitools import utils
from astropy.io import fits

###trim off the crap parts of the KCWI cube medium slicer
###This is *aggressive* in the y direction because of dispersion
def kcwi_cube_trim_Medium(infil):
    hdr, flux = io.open_kcwi_cube(infil)
    trimflux = flux[:, 18:81, 6:28]  # note that python reads arrays weird.  Trim *down in y* then x.

    a, b = infil.split(".fits")
    outfil = a + '_trimmed.fits'

    hdu_out = fits.PrimaryHDU(trimflux, header=hdr)
    hdu_out.writeto(outfil, overwrite=True)


###trim off crap parts of larger kcwi cube
def kcwi_cube_trim_large(infil):
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
        whiteim: ndarray in units of 10^-16 (erg/s/cm2/arcsec2)

    """
    # Wavelength
    wave = utils.build_wave(hdr)
    # Slices
    slices = np.where((wave >= minwave) & (wave <= maxwave))[0]
    # Do it
    header =kp.tweak_header(deepcopy(hdr))
    #Modify and sum data to be in units of 10^-16 (erg/s/cm2/arcsec2)
    plate_scale = proj_plane_pixel_scales(WCS(header))
    plate_scale*= 3600. #arcsec
    pixArea = plate_scale[0]*plate_scale[1] #area in arcsec^2
    dLambda = wave[2]-wave[1]#(maxwave-minwave)


    whiteim = np.sum(flux[slices,:,:], axis=0)

    # convert units
    whiteim /= pixArea
    whiteim *= dLambda  

    #clipped = sigmaclip(whiteim).clipped
    #med = np.median(clipped)
    #whiteim -= med


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
        nbimage: ndarray in units of 10^-16 (erg/s/cm2/arcsec2)

    """
    # Rest wave
    if z is None:
        z = 0.
    # Wavelengths
    wave = utils.build_wave(hdr)
    wrest = wave/(1+z)

    header =kp.tweak_header(deepcopy(hdr))
    #Modify and sum data to be in units of 10^-16 (erg/s/cm2/arcsec2)
    plate_scale = proj_plane_pixel_scales(WCS(header))
    plate_scale*= 3600. #arcsec
    pixArea = plate_scale[0]*plate_scale[1] #area in arcsec^2
    dLambda = wave[2]-wave[1]#del_wave delta_lambda_per_pixel


    # Work out the slices for the on-band image:
    slices = np.abs(wrest - line) < del_wave

    # define windows on either side of the slice (watch out for a,b and <= vs <)
    offimage_width = 2. * del_wave

    slice_low = np.where((wrest >= (line - del_wave - offimage_width)) &
                         (wrest < (line - del_wave)))[0]
    slice_high = np.where((wrest >= (line + del_wave)) &
                          (wrest < (line + del_wave + offimage_width)))[0]

    # Make the on-band and two off-band images
    nbimage = np.sum(flux[slices, :, :], 0)
    high = np.sum(flux[slice_high, :, :], 0)
    low = np.sum(flux[slice_low, :, :], 0)
    
    # convert units
    # If flambda is in units of erg s^-1 cm^-2 A^-1 
    # get SB in units of erg s^-1 cm^-2 arcsec^-2
    nbimage /= pixArea
    nbimage *= dLambda  

    high /= pixArea
    high *= dLambda  
    low /= pixArea
    low *= dLambda  

    # Average the two off-band images and subtract.
    if sub_offimage:
        offimage = (low + high) / 2.0
        nbimage -= offimage

    # Output?
    if outfile is not None:
        io.write_image(hdr, nbimage, outfile)

    # Return
    return nbimage


def cube_skysub(fil, skyy1, skyy2, skyx1, skyx2, outfil=None):
    hdr, flux = io.open_kcwi_cube(fil)
    # Wavelengths
    wave= utils.build_wave(hdr)

    wavedim, ydim, xdim = flux.shape
    sky = np.zeros(wavedim)

    for i in range(wavedim):
        sky[i] = np.nanmedian(flux[i, skyy1:skyy2, skyx1:skyx2])

    # Now Sky subtract the whole cube
    fluxmsky = np.zeros((wavedim, ydim, xdim))
    for i in range(wavedim):
        fluxmsky[i, :, :] = flux[i, :, :] - sky[i]


    if outfil is not None:
        hdu_out = fits.PrimaryHDU(fluxmsky, header=hdr)
        hdu_out.writeto(outfil, overwrite=True)

    #Returns
    return fluxmsky


def var_skysub(varfil, skyy1, skyy2, skyx1, skyx2, outfil=None):
    hdr, var = io.open_kcwi_cube(varfil)
    # Wavelengths
    wave= utils.build_wave(hdr)
    wavedim, ydim, xdim = var.shape

    sky = np.zeros(wavedim)

    for i in range(wavedim):
        sky[i] = np.nanmedian(var[i, skyy1:skyy2, skyx1:skyx2])

    # Now Sky subtract the whole cube
    varwsky = np.zeros((wavedim, ydim, xdim))
    for i in range(wavedim):
        varwsky[i, :, :] = var[i, :, :] + sky[i]

    if outfil is not None:
        hdu_out = fits.PrimaryHDU(varwsky, header=hdr)
        hdu_out.writeto(outfil, overwrite=True)

    #Returns
    return varwsky

