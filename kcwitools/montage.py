""" montage routines"""
from __future__ import print_function, absolute_import, division, unicode_literals
import numpy as np

from astropy.io import fits
from astropy import units
import subprocess

def run_montage(infils,outfil="Montage.fits",grating='BL',clean=False):
    """ take a list of (ideally trimmed) KCWI cubes and run montage on them
    Args:
    ----------
    infils: list
      list of kcwi cubes to montage
    outfil: string
      name of montage output file, defaults to Montage.fits'

    
    """

    #create the driectories
    subprocess.call(["mkdir","Input"])
    subprocess.call(["mkdir","Projection"])

    #copy
    for fil in infils:
        subprocess.call(["cp",fil,"Input"])

    #first part of montage
    subprocess.call(["mImgtbl","-c","Input/","cubes.tbl"],shell=True)
    subprocess.call(["mMakeHdr","cubes.tbl","cubes.hdr"],shell=True)

    #second art of montage
    for fil in infils:
        tmp=fil.split(".")
        subprocess.call(["mProjectCube","Input/"+fil,"projection/"+tmp[0]+"_proj.fits","cubes.hdr"],shell=True)

    #final part of montage
    subprocess.call(["mImgtbl","-c","projection/","cubes-proj.tbl"],shell=True)
    subprocess.call(["mAddCube","-p projection/","cubes-proj.tbl","cubes.hdr",outfil],shell=True)

    #remove all montage created bits aside from output file (default is off)
    if(clean):
        subprocess.call(["rm","-rf","Input"])
        subprocess.call(["rm","-rf","Projection"])
        subprocess.call(["rm","cubes.tbl","cubes.hdr","cubes-proj.tbl"])

    #fix the header after montage is done
    if grating == 'BL':
        fix_kcwi_cube_montage_BL(outfil)
    if grating == 'BM':
        fix_kcwi_cube_montage_BM(outfil)
        

###Add the CD3_3 and CDELT3 keywords to a BL grating cube cube made by montage
def fix_kcwi_cube_montage_BL(mosaicfil):
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
def fix_kcwi_cube_montage_BM(mosaicfil):
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


    
