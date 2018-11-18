""" RA/DEC and WCS fussing for KCWI """
from __future__ import print_function, absolute_import, division, unicode_literals
import numpy as np
import pdb

from astropy import units

from kcwitools import io
from kcwitools import utils


def offsets(coord1, coord2, verbose=True):
    """
    Input a pair of RA/DEC and calculate the RA/DEC offsets between them

    Taken from xastropy

    Parameters:
    ----------
    irad1 : SkyCoord
      RA/DEC of source 1 (origin)
    irad2 : SkyCoord
      RA/DEC of source 2 (destination)
    verbose: bool, optional

    Returns:
    -------
    offsets, PA : (Quantity, Quantity), Quntity
       Tuple of offsets (arcsec)
       Position Angle (degrees)
    """

    # Angular separation
    sep = coord1.separation(coord2).to('arcsec')

    # PA
    PA = coord1.position_angle(coord2)

    # RA/DEC
    dec_off = np.cos(PA) * sep # arcsec
    ra_off = np.sin(PA) * sep # arcsec (East is *higher* RA)

    # Print
    if verbose:
        print('RA Offset from 1 to 2 is {:g}'.format(ra_off))
        print('DEC Offset from 1 to 2 is {:g}'.format(dec_off))
        print('PA = {:g}'.format(PA.degree*units.degree))

    # Return
    return (ra_off, dec_off), PA.degree * units.degree


def offset_radec(hdu_hdr, ra_offset, dec_offset, coord):
    """
    Impose offset on the header given offsets in *arcseconds*

    Args:
        hdu_hdr: Header
        ra_offset: float (deg)
          Arcseconds
        dec_offset:  float (deg)
          Arcseconds
        coord: SkyCoord
          Needed to translatate arcseconds to degrees of RA

    Returns:
        new_hdr : Header
    """
    ra_off = ra_offset.to('deg').value/np.cos(coord.dec).value
    dec_off = dec_offset.to('deg').value
    #
    new_hdr = hdu_hdr.copy()
    # Offset
    new_hdr['CRVAL1'] += ra_off
    new_hdr['CRVAL2'] += dec_off
    # Return
    return new_hdr


def set_radec(hdu_hdr, ra, dec, x, y):
    """
    Modify the Header to force RA/DEC value at a given x,y position

    Args:
        hdu_hdr: Header
        ra: float (deg)
        dec: float (deg)
        x: float
          x pixel position in the image
        y: float
          y pixel position in the image

    Returns:
        new_hdr: Header

    """
    new_hdr = hdu_hdr.copy()
    # Offset
    new_hdr['CRVAL1'] = ra
    new_hdr['CRVAL2'] = dec
    new_hdr['CRPIX1'] = x
    new_hdr['CRPIX2'] = y
    return new_hdr
