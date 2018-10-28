from __future__ import print_function

import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
from kcwitools import optimal_extract_utility as r
from linetools.spectra.xspectrum1d import XSpectrum1D


def extract_optimal_1D(fluxcube, varcube,header,wavegrid,x1,x2,y1,y2,outfile=None,verbose=False,plot=False):
	"""

	Args:
		fluxcube: 3D kcwi fluxcube
		varcube: 3D KCWI variance cube
		header:  header information
		wavegrid:  1D kcwi wavelength grid
		
		x1,x2,y1,y2: Four corners of an aribaty rectangular pseudo-slit for extraction
		outfile: str, optional
		  FITS file name
		verbose [False] : optional, prints information on terminal
		plot [False] : optional plots diagnostic stuff
	
	Returns:
		xspec: XSpectrum1D

	"""	

	# CHECK CHECK CHECK
	#THESE VALUES NEED TO BE DOUBLE CHECKED. WHICH ONE SHOULD I USE?
	gain=header['GAIN2']
	rdnoise=3.5  #  Don't know this value. Please instruct what should I read in here?


	# Hardcoding some values
	sigma=5  # Sigma Clipping Value
	fit_order=200   # Chebyshev polinomial order


	# Create slices within wavelength ranges that are good for KCWI
	# From now on only deal with these smaller slices
	slices=np.where((wavegrid >=3500.) & (wavegrid <= 5500.))
	flux=fluxcube[slices[0],:,:]
	wave=wavegrid[slices[0]]
	var_origin=(varcube[slices[0],:,:])



	# Get the required shapes
	wavedim, ydim, xdim = flux.shape


	# Plot and show the slit position on the white light image
	if plot==True:
		whitelight=np.sum(flux, axis=0)
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


	N=len(aperture_indices[0])


	# Do a boxcar extraction
	spec=r.extract1d(cube)
	var=var_origin+rdnoise**2


	cube_bad=np.where(~np.isfinite(cube))
	cube[cube_bad]=0
	cube[cube<0]=0

	varbad=np.where(~np.isfinite(var))
	var[varbad]=np.inf

	sum_spec=np.sum(np.sum(aperture_mask*cube, axis=2), axis=1)
	sum_err= np.sqrt(np.abs(np.sum(np.sum(var,axis=2),axis=1)))



	# Now starting Optimal Extraction
	counter=1.0
	c=1
	#Use the cube of pixel weights, if given

	while counter !=0:
		if verbose == True:
			print("Wavelength Iteration {}".format(c))
			print("\n\nCONSTRUCTING SPATIAL PROFILE\n")
		mcube=r.wavelegth_fitter( cube, var, aperture_mask, wave, fit_order, aperture_indices, plot=False, verbose=False)
		
		if verbose == True:
			print("\n\nREVISING VARIANCE ESTIMATES\n")

		var=r.variance(mcube, spec, rdnoise, gain)
		varbad=np.where(~np.isfinite(var))
		var[varbad]=np.inf
	
		if verbose==True:
			print("\n\nMASKING COSMICS/BAD PIXELS\n")
		aperture_mask, counter=r.sigma_clip3D(cube, var, mcube, spec, aperture_mask, sigma, aperture_indices, verbose=False)
		var[~aperture_mask]=np.inf
		cube[~aperture_mask]=0

		if verbose==True:
			print("\n\nEXTRACTING OPTIMAL SPECTRUM\n")
		spec=np.sum(np.sum(aperture_mask*cube*mcube/var, axis=2), axis=1)/np.sum(np.sum(aperture_mask*mcube*mcube/var, axis=2), axis=1)
		vspec=np.sum(np.sum(aperture_mask*mcube, axis=2), axis=1)/np.sum(np.sum(aperture_mask*mcube*mcube/var , axis=2), axis=1)
		c+=1

	# Object
	xspec = XSpectrum1D.from_tuple((wave, spec, np.sqrt(vspec)))

	# Object
	xspec_std = XSpectrum1D.from_tuple((wave, sum_spec, sum_err))
	

	if plot==True:
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
	ind= int(np.floor( (xdim+1) * (x1 - xmin)/ (xmax - xmin)  ))
	return ind
