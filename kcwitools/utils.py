""" utility routines"""
from __future__ import print_function, absolute_import, division, unicode_literals
import numpy as np

from astropy import units
import subprocess
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
    #wave = crval3 + (crpix3 + np.arange(0, wavedim, 1.0)) * cd3_3
    wave = crval3 + cd3_3 * (np.arange(wavedim) + 1. - crpix3)
    # Air2vac?
    if air2vac:
        wave = airtovac(wave)
    # Return
    return wave




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
    crval3 = hdu_hdr['CRVAL3']
    crpix3 = hdu_hdr['CRPIX3']
    cd3_3 = hdu_hdr['CD3_3']
    wavedim = hdu_hdr['NAXIS3']


    
    wavedim, ydim, xdim = fluxmsky.shape
    
    header =kp.tweak_header(deepcopy(hdu_hdr))
    headerv =kp.tweak_header(deepcopy(hdrv))
    

    
    q=np.isnan(fluxmsky)
    fluxmsky[q]=0.
    qq=np.isnan(var)
    fluxmsky[q]=0.
    #position = (16.5, 47.5)
    #size = (66, 25)     # pixels


    newflux=np.zeros((wavedim,size[0],size[1]))
    newvar=np.zeros((wavedim,size[0],size[1]))

    cutout=Cutout2D(fluxmsky[0,:,:], position, size, wcs=WCS(header))
    cutout_v=Cutout2D(var[0,:,:], position, size, wcs=WCS(headerv))
    
    
    index=cutout.bbox_original
    index1=cutout_v.bbox_original
    
    newflux=fluxmsky[:,index[0][0]:index[0][1],index[1][0]:index[1][1]]
    newvar=var[:,index1[0][0]:index1[0][1],index1[1][0]:index1[1][1]]

    
    
    
    #for i in range(0,wavedim):
    #    cutout = Cutout2D(fluxmsky[i,:,:], position, size, wcs=WCS(header))
    #    cutout_v = Cutout2D(var[i,:,:], position, size, wcs=WCS(headerv))
    #    #pdb.set_trace()
    #    newflux[i,:,:]=cutout.data
    #    newvar[i,:,:]=cutout_v.data
    
    
    
    # Construct new header
    
    temp=cutout.wcs.to_header()
    temp1=cutout_v.wcs.to_header()

    
    newhdr_flx=deepcopy(hdu_hdr)
    newhdr_flx['NAXIS1']=size[1]
    newhdr_flx['NAXIS2']=size[0]
    newhdr_flx['CRPIX1']=temp['CRPIX1']
    newhdr_flx['CRPIX2']=temp['CRPIX2']
    newhdr_flx['CRVAL3']=crval3
    newhdr_flx['CRPIX3']=crpix3
    newhdr_flx['CD3_3']=cd3_3
    newhdr_flx['NAXIS3']=wavedim

    
    
    newhrd_var=deepcopy(hdrv)
    newhrd_var['NAXIS1']=size[1]
    newhrd_var['NAXIS2']=size[0]
    newhrd_var['CRPIX1']=temp1['CRPIX1']
    newhrd_var['CRPIX2']=temp1['CRPIX2']
    newhrd_var['CRVAL3']=crval3
    newhrd_var['CRPIX3']=crpix3
    newhrd_var['CD3_3']=cd3_3
    newhrd_var['NAXIS3']=wavedim

        
    return newflux,newvar,newhdr_flx,newhrd_var




def build_moment_map(hdr, flux, line, z=None, del_wave=2.0):
    """ Simple 1st moment of an emission line to compute velocity

    Parameters:
    ----------
    hdr    : header of the file needed
    flux   : flux datacube matrix
    line   : rest frame wavelength of line center
    z      : line redshift [defult 0]
    del_wave: delta wavelength window within which 1st moment is computed

    Returns:
    ----------
    nbimage: new image with 1st moment map
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
    dLambda = del_wave


    # Work out the slices for the on-band image:
    slices = np.abs(wrest - line) < del_wave
    
    wave_grid=wrest[slices]
    flux_cube=flux[slices, :, :]
    
    nbimage=np.zeros((np.shape(flux_cube)[1],np.shape(flux_cube)[2]))
    z=  (wave_grid/line) -1.#;
    C = 299792.458#;  % speed of light [km/sec]
    Beta = ((z+1.)**2. - 1.)/(1. + (z+1.)**2.);
    del_v = Beta*C;



    for i in range(0,(np.shape(flux_cube)[1])):
        for j in range(0,(np.shape(flux_cube)[2])):
            # Make the on-band and two off-band images
            nbimage[i,j] = np.trapz(del_v*flux_cube[:,i,j],x=del_v,axis=0)
    

    return (nbimage)


