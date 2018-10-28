from __future__ import print_function

import numpy as np
import matplotlib.pyplot as plt
from numpy.polynomial.chebyshev import chebfit,chebval



# Adopted from some online stuff that uses these for optimal extraction

def wavelegth_fitter( fluxcube, varcube, aperture_mask, wavelength, order, aperture_indices, plot=False, verbose=False):

	"""
	Fit polynomials in the wavelength direction of each pixel in an aperture. 
	
	Input: 
		 fluxcube 			:    KCWI flux cube 
		 varceub  			:	 KCWI variance cube
		 aperture_mask  	:    3D mask for the pseudo-slit. Good values=True, Bad ones =False
		 wavelength     	:    wavelength grid of the fluxcube
		 order          	:    Chebyshev polynomial order for fitting
		 aperture_indices 	:    Indices of the cube corresponding to pseudo-slit aperture, in an Nx2 array. 

		 plot [False]		:  If true plots the polynomials (both before and after normalising spatially)    
		 verbose [False]    :  True prints stuff in terminal


	"""


	wavedim,ydim,xdim=fluxcube.shape
	wavelength_fits=np.zeros([wavedim, ydim, xdim])


	#Lamdas must be scaled to between -1 and 1 for Chebyshev fitting
	min_val=wavelength.min()
	max_val=wavelength.max()
	scaled_l=np.array([2*(lamda-min_val)/(max_val-min_val) -1 for lamda in wavelength])

	for n, (i, j) in enumerate(zip(aperture_indices[0], aperture_indices[1])):

		"""Take each spaxel and fit a polynomial along the wavelength direction, each pixel weighted by one over its variance"""
		if verbose==True:
			print("Fitting Pixel {}".format(n))

		weights=1.0/varcube[:, i, j]

		#Fit the Chebyshev coefficients
		coefficients=chebfit(scaled_l, fluxcube[:, i, j], order, w=weights)

		#Make a polynomial from these
		polynomial=chebval(scaled_l, coefficients)
		#Ensure all values are positve
		polynomial[polynomial<0]=0


		wavelength_fits[:, i, j]=polynomial

	if plot==True:
		for i, j in zip(aperture_indices[0], aperture_indices[1]):
			plt.plot(wavelength, wavelength_fits[:, :])

		plt.show()

	if verbose==True:
		print("Normalising Spatially")
	#Normalise spatially
	spatial_norm=np.sum(np.sum(wavelength_fits, axis=2), axis=1)

	# This is the P(x, lambda) value from Horne 1986 in 3D
	mcube=wavelength_fits/spatial_norm[:, np.newaxis, np.newaxis]



	if plot==True:
		for i, j in zip(aperture_indices[0], aperture_indices[1]):

			plt.plot(wavelength, mcube[:,i, j])

		plt.show()

	return mcube



def extract1d(datacube,verbose=False):
	if verbose==True:
		print("Boxcar extraction\n\n")
	#boxcar extraction
	spectrum=np.sum(np.sum(datacube, axis=2), axis=1)
	return spectrum


def variance(mcube, spec, rdnoise, gain,verbose=False):
	"""Step 6 in Horne 1986"""
	if verbose==True:
		print("Updating the Variance Cube")

	variance=(np.abs(spec[:, np.newaxis, np.newaxis]*mcube)/gain)+rdnoise**2

	return variance


def sigma_clip3D(datacube, varcube, mcube, spec, aperture_mask, sigma, aperture_indices, verbose=False):

	"""
	Perform a sigma clipping on residuals
	"""
	if verbose == True:
		print(" Sigma clipping outliers")

	#Forming the residuals
	tab=datacube-spec[:, np.newaxis, np.newaxis]*mcube
	residuals=(aperture_mask)*np.power(tab, 2.)/varcube
	

	#Find the outliers
	outliers=np.where(residuals>sigma**2.)
	noutliers=len(residuals[outliers])

	if verbose == True:
		print("There are {} outliers".format(noutliers))
	
	#Flip the True values to False
	aperture_mask[outliers]=~aperture_mask[outliers]


	return (aperture_mask, noutliers)
