""" utility routines"""
from __future__ import print_function, absolute_import, division, unicode_literals
import numpy as np

from astropy import units
from astropy.io import fits
from astropy.wcs import WCS
from astropy.nddata import Cutout2D
from copy import deepcopy
from kcwitools import io as ki
from kcwitools import plot as kp


def build_wave(hdu_hdr, air2vac=True):
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
        wave = airtovac(wave)
    # Return
    return wave

###Add the CD3_3 and CDELT3 keywords to a cube made by montage
def fix_kcwi_cube_montage(mosaicfil):
    hdu = fits.open(mosaicfil)
    flux = hdu['PRIMARY'].data
    hdu_hdr = hdu['PRIMARY'].header
    crval3 = hdu_hdr['CRVAL3']
    crpix3 = hdu_hdr['CRPIX3']
    cdelt3 = hdu_hdr['CDELT3']

    hdu_hdr['CDELT3'] = 1.0
    hdu_hdr['CD3_3'] = 1.0

    hdu_out = fits.PrimaryHDU(flux, header=hdu_hdr)
    hdu_out.writeto(mosaicfil, overwrite=True)


###Add the CD3_3 and CDELT3 keywords to a cube made by montage for higher resolution
def fix_kcwi_cube_montage_2(mosaicfil):
    hdu = fits.open(mosaicfil)
    flux = hdu['PRIMARY'].data
    hdu_hdr = hdu['PRIMARY'].header
    crval3 = hdu_hdr['CRVAL3']
    crpix3 = hdu_hdr['CRPIX3']
    cdelt3 = hdu_hdr['CDELT3']

    hdu_hdr['CDELT3'] = 0.5
    hdu_hdr['CD3_3'] = 0.5

    hdu_out = fits.PrimaryHDU(flux, header=hdu_hdr)
    hdu_out.writeto(mosaicfil, overwrite=True)


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


def trim_kcwi_cube_with_header(fluxfile,varfile,size = (66, 25), position = (16.5, 47.5)):
     """ Trims KCWI datacube to remove bad pixels and reconstructs the header

    Parameters:
    ----------
    fluxfile  : name of the flux file
    varfile   : name of the variance file 
    size= (y,x) size of the cutout image in pixels
    position = (y,x) central position of the cutout in pixels

    Returns:
    ----------
    newflux, newvar :  New trimmed flux and variance cubes
    newhdr_flx,newhrd_var: new updated headers
    """

    hdu_hdr,fluxmsky = ki.open_kcwi_cube(fluxfile)
    hdrv,var = ki.open_kcwi_cube(varfile)
    wave = build_wave(hdu_hdr)
    
    wavedim, ydim, xdim = fluxmsky.shape
    
    header =kp.tweak_header(hdu_hdr)
    headerv =kp.tweak_header(hdrv)
    

    
    q=np.isnan(fluxmsky)
    fluxmsky[q]=0.
    qq=np.isnan(var)
    fluxmsky[q]=0.
    #position = (16.5, 47.5)
    #size = (66, 25)     # pixels


    newflux=np.zeros((wavedim,size[0],size[1]))
    newvar=np.zeros((wavedim,size[0],size[1]))
    
    
    
    for i in range(0,wavedim):
        cutout = Cutout2D(fluxmsky[i,:,:], position, size, wcs=WCS(header))
        cutout_v = Cutout2D(var[i,:,:], position, size, wcs=WCS(hdrv))
        #pdb.set_trace()
        newflux[i,:,:]=cutout.data
        newvar[i,:,:]=cutout_v.data
    
    
    
    # Construct new header
    
    temp=cutout.wcs.to_header()
    temp1=cutout_v.wcs.to_header()

    
    newhdr_flx=deepcopy(hdu_hdr)
    newhdr_flx['NAXIS1']=size[1]
    newhdr_flx['NAXIS2']=size[0]
    newhdr_flx['CRPIX1']=temp['CRPIX1']
    newhdr_flx['CRPIX2']=temp['CRPIX2']
    
    
    newhrd_var=deepcopy(hdrv)
    newhrd_var['NAXIS1']=size[1]
    newhrd_var['NAXIS2']=size[0]
    newhrd_var['CRPIX1']=temp1['CRPIX1']
    newhrd_var['CRPIX2']=temp1['CRPIX2']
    
    
    
    return newflux,newvar,newhdr_flx,newhrd_var



