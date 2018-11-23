""" utility routines"""
from __future__ import print_function, absolute_import, division, unicode_literals
import numpy as np
import pdb

from astropy import units
from astropy.coordinates import SkyCoord
from astropy.time import Time

from pypeit.core import wave as pyp_wave


def build_wave(hdu_hdr, air2vac=True, helio=True):
    """

    Args:
        hdu_hdr:

    Returns:
        wave: ndarray

    """
    crval3 = hdu_hdr['CRVAL3']
    crpix3 = hdu_hdr['CRPIX3']
    cd3_3 = hdu_hdr['CD3_3']
    wavedim = hdu_hdr['NAXIS3']
    # Do it
    wave = crval3 + (crpix3 + np.arange(0, wavedim, 1.0)) * cd3_3
    # Air2vac?
    if air2vac:
        print("Correcting wavelengths to vacuum")
        wave = airtovac(wave)
    # Helio-centric?
    if helio and 'MJD' in hdu_hdr.keys():
        print("Correcting wavelengths to helio-centric")
        vel_corr = helio_correction(hdu_hdr)
        wave *= vel_corr
    # Return
    return wave


def helio_correction(hdr):
    # Keck
    longitude=155.47833
    latitude=19.82833
    altitude=4160.0
    loc = (longitude * units.deg, latitude * units.deg, altitude * units.m,)
    # Grab coord
    radec = SkyCoord(ra=hdr['CRVAL1'], dec=hdr['CRVAL2'], unit='deg')
    tval = Time(hdr['MJD'], scale='tt', format='mjd')
    obstime = Time(tval.value, format=tval.format, scale='utc', location=loc)
    # Calculate
    vel = pyp_wave.geomotion_velocity(obstime, radec, frame='heliocentric')
    vel_corr = np.sqrt((1. + vel/299792.458) / (1. - vel/299792.458))
    # Return
    return vel_corr


def airtovac(wave):
    """ Convert air-based wavelengths to vacuum

    Parameters:
    ----------
    wave: ndarray
      Wavelengths

    Returns:
    ----------
    wavelenght: ndarray
      Wavelength array corrected to vacuum wavelengths
    """
    # Assume AA
    wavelength = wave

    # Standard conversion format
    sigma_sq = (1.e4/wavelength)**2. #wavenumber squared
    factor = 1 + (5.792105e-2/(238.0185-sigma_sq)) + (1.67918e-3/(57.362-sigma_sq))
    factor = factor*(wavelength>=2000.) + 1.*(wavelength<2000.) #only modify above 2000A

    # Convert
    wavelength = wavelength*factor

    return wavelength
