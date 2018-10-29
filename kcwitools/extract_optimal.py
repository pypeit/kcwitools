""" Routines for extracting optimal 1D spectrum using Horne 1986 algorithm"""
from __future__ import print_function, absolute_import, division, unicode_literals
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
import pdb

from kcwitools import optimal_extract_utility as r
from kcwitools import utils
from kcwitools import defs
from kcwitools import image
from linetools.spectra.xspectrum1d import XSpectrum1D


def extract_optimal_1D(fluxcube, varcube,header,corners,sigma=5,fit_order = 200,outfile=None,verbose=False,plot=False):
    """

    Args:
        fluxcube         : 3D kcwi fluxcube
        varcube          : 3D KCWI variance cube
        header           : Header information       
        corners          : A tuple defining x1,x2,y1,y2: Four corners of an aribaty rectangular pseudo-slit for extraction
        sigma            : Sigma clipping values [default 5]
        fit_order        : polinomial fit order for spatial profile [default 200]
        
        outfile          : str, optional
                          FITS file name
        verbose [False]  : optional, prints information on terminal
        plot [False]     : optional plots diagnostic stuff
    
    Returns:
        xspec            : XSpectrum1D object [optimal extraction]
        xspec_std        : XSpectrum1D object [boxcar extraction]

    """ 

    # CHECK CHECK CHECK
    #THESE VALUES NEED TO BE DOUBLE CHECKED. WHICH ONE SHOULD I USE?


    val= defs.defs(header)

    gain=val['gain']
    rdnoise=val['rdnoise']  #  Don't know this value. Please instruct what should I read in here?

    #if val['BIASSUB'] == 'T':  # If read noise 



    # Wavelength
    wavegrid = utils.build_wave(header)

    # Expand the four corners of the pseudo-slit
    x1=corners[0]
    x2=corners[1]
    y1=corners[2]
    y2=corners[3]



    # Create slices within wavelength ranges that are good for KCWI
    # From now on only deal with these smaller slices
    #slices=np.where((wavegrid >=3500.) & (wavegrid <= 5500.))
    slices=(wavegrid >=3500.) & (wavegrid <= 5500.)

    flux=fluxcube[slices,:,:]
    wave=wavegrid[slices]
    var_origin=(varcube[slices,:,:])



    # Get the required shapes
    wavedim, ydim, xdim = flux.shape


    #Changing header Naxis3 to reflect the new sliced cube
    header['NAXIS3'] = int(sum(slices))


    # Plot and show the slit position on the white light image
    if plot is True:
        whitelight=image.build_whitelight(header,flux)
        plt.figure()
        plt.imshow(whitelight,cmap=plt.get_cmap('viridis'),origin="lower")
        plt.colorbar()

        # Showing the slit position
        plt.plot([y1,y1,y2,y2,y1],[x1,x2,x2,x1,x1],'g')
        plt.show()

    # create a mask
    indx1=mask_index(xdim,y1,xdim,0)
    indx2=mask_index(xdim,y2,xdim,0)
    indy1=mask_index(ydim,x1,ydim,0)
    indy2=mask_index(ydim,x2,ydim,0)

    #Making a 3d masked cube, rather than just a masked slice
    #True=Good, False=Bad

    mask = np.zeros((wavedim, ydim, xdim))
    mask[:,indy1:indy2,indx1:indx2]=1

    # Create a boolean 3D True False mask
    aperture_mask=np.isin(mask, 1)
    # This is the masked cube that can be used for extraction
    cube=flux*mask
    var_origin[~aperture_mask]=np.inf

    #Get the indices of the aperture spaxels at each wavelength
    aperture_indices = np.where((aperture_mask[900,:,:] == True))
    #pdb.set_trace()
    

    N=len(aperture_indices[0])


    # Do a boxcar extraction
    spec=np.sum(np.sum(cube, axis=2), axis=1)
    var=var_origin+rdnoise**2


    cube_bad=(~np.isfinite(cube))
    cube[cube_bad]=0
    cube[cube<0]=0

    varbad=(~np.isfinite(var))
    var[varbad]=np.inf

    sum_spec=np.sum(np.sum(aperture_mask*cube, axis=2), axis=1)
    sum_err= np.sqrt(np.abs(np.sum(np.sum(var,axis=2),axis=1)))



    # Now starting Optimal Extraction
    counter=1
    c=1
    #Use the cube of pixel weights, if given

    while counter !=0:
        if verbose is True:
            print("Wavelength Iteration {}".format(c))
            print("\n\nCONSTRUCTING SPATIAL PROFILE\n")
        mcube=r.spatial_profile( cube, var, aperture_mask, wave, fit_order, aperture_indices, plot=False, verbose=False)
        
        if verbose is True:
            print("\n\nREVISING VARIANCE ESTIMATES\n")

        var=r.variance(mcube, spec, rdnoise, gain)
        varbad=(~np.isfinite(var))
        var[varbad]=np.inf
    
        if verbose is True:
            print("\n\nMASKING COSMICS/BAD PIXELS\n")
        aperture_mask, counter=r.sigma_clip3D(cube, var, mcube, spec, aperture_mask, sigma, aperture_indices, verbose=False)
        var[~aperture_mask]=np.inf
        cube[~aperture_mask]=0

        if verbose is True:
            print("\n\nEXTRACTING OPTIMAL SPECTRUM\n")
        spec=np.sum(np.sum(aperture_mask*cube*mcube/var, axis=2), axis=1)/np.sum(np.sum(aperture_mask*mcube*mcube/var, axis=2), axis=1)
        vspec=np.sum(np.sum(aperture_mask*mcube, axis=2), axis=1)/np.sum(np.sum(aperture_mask*mcube*mcube/var , axis=2), axis=1)
        c+=1

    # Object
    xspec = XSpectrum1D.from_tuple((wave, spec, np.sqrt(vspec)))

    # Object
    xspec_std = XSpectrum1D.from_tuple((wave, sum_spec, sum_err))
    

    if plot is True:
        plt.figure()
        plt.title("Optimal Extract (Blue) and Normal (red)")
        plt.plot()
        plt.plot(wave, spec, c="b")
        plt.plot(wave,sum_spec,c="r")
        plt.show()
    # Object
    
    if outfile is not None:
        xspec.write_to_fits(outfile)
    
    return xspec,xspec_std



def mask_index(xdim,x1,xmax,xmin):
    """ 
    Creates a quick indexs for a mask
       Args:
          xdim  :  dimension of the vector
          x1    :  x value for which index is sought
          xmax  :  maximum value of x axis 
          xmin  :  minimum value of x axis 

       Returns :
            ind  : Final indexes for x1    
    """ 
    ind= int(np.floor( (xdim+1) * (x1 - xmin)/ (xmax - xmin)  ))
    return ind
