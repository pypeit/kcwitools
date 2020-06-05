from __future__ import print_function, absolute_import, division, unicode_literals
from astropy.modeling import models, fitting
from linetools.spectra.xspectrum1d import XSpectrum1D
import numpy as np
import matplotlib.pyplot as plt
import warnings

def extract_weighted_spectrum(flux,variance,wave,verbose=False,weights='Gaussian',porder=9):
    """
    Args:
        flux             : 3D kcwi fluxcube in the aperture of extraction
        variance          : 3D KCWI variance cube in the aperture of extraction
        wave             : wavelength of the datacube that is fed in     
        verbose          : [default False]: if set True shows the fitted light profile
        weights          : if = 'Gaussian',fits a Gaussian profile
                              = 'poly', fits an n order polynomial where n is set by porder
                              = 'Data', Weights by the raw white light image
        porder           :  polynomial order [default 9]
        
     
    Returns:
        xspec            : XSpectrum1D object [optimal extraction]

    """
    img=flux.sum(axis=0)
    img[img<0] =0
    amp_init = flux.max()
    ydim,xdim=img.shape
    stdev_init_x = 0.33 * ydim
    stdev_init_y = 0.33 * xdim
    theta=0
    
    if weights=='gaussian':
        g_init = models.Gaussian2D(amp_init, 5, 5, stdev_init_x, stdev_init_y,theta=10)
    elif weights =='poly':
        g_init =models.Polynomial2D(degree=porder)# 
    else:
        g_init = models.Gaussian2D(amp_init, 5, 5, stdev_init_x, stdev_init_y,theta=10)


    yi, xi = np.indices(img.shape)

    fit_g = fitting.LevMarLSQFitter()
    with warnings.catch_warnings():
        # Ignore model linearity warning from the fitter
        warnings.simplefilter('ignore')
        p = fit_g(p_init, xi, yi, img)
        
    if verbose == True:
        # Plot the data with the best-fit model
        print('Plotting fitted model profile')
        plt.figure(figsize=(8, 2.5))
        plt.subplot(1, 3, 1)
        plt.imshow(img, origin='lower', interpolation='nearest')#, vmin=0, vmax=1.2e3)\
        plt.title("Data")
        plt.subplot(1, 3, 2)
        plt.imshow(p(xi, yi), origin='lower', interpolation='nearest')#, vmin=0, vmax=1.2e3)#, vmin=-1e4,
           #vmax=5e4)
        plt.title("Model")
        plt.subplot(1, 3, 3)
        plt.imshow(img - p(xi, yi), origin='lower', interpolation='nearest')#, vmin=-10, vmax=10)#, vmin=-1e4,
         #  vmax=5e4)
        plt.title("Residual")
    
    if weights =='Data':
        weights=img
    else:
        weights = p(xi, yi)

        
        
    w = wave
    n = len(w)
    fl = np.zeros(n)
    sig = np.zeros(n)
    
    for wv_ii in range(n):
        # n_spaxels = np.sum(mask)
        weights = weights / np.sum(weights)
        fl[wv_ii] = np.nansum(flux[wv_ii] * weights)  # * n_spaxels
        sig[wv_ii] = np.sqrt(np.nansum(variance[wv_ii] * (weights ** 2)))  # * n_spaxels
        
    # renormalize
    fl_sum = np.nansum(flux,axis=(1,2))
    norm = np.sum(fl_sum) / np.sum(fl)
    fl = fl * norm
    sig = sig * norm
    
     # Object
    xspec = XSpectrum1D.from_tuple((wave, fl, sig))

    return xspec
