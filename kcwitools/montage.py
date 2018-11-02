""" montage routines"""
from __future__ import print_function, absolute_import, division, unicode_literals
import numpy as np

from astropy.io import fits
from astropy import units
import subprocess

def run_montage(infils,outdir="./",outfil="Montage.fits",grating='BL',clean=False):
    """ take a list of (ideally trimmed) KCWI cubes and run montage on them
    Args:
    ----------
    infils: list
      list of kcwi cubes to montage
    outfil: string
      name of montage output file, defaults to Montage.fits'

    
    """

    #create the driectories
    subprocess.run(["mkdir",outdir+"Input"])
    subprocess.run(["mkdir",outdir+"Projection"])

    #copy
    for fil in infils:
        subprocess.run(["cp",fil,outdir+"Input"],shell=True)

    #first part of montage
    subprocess.call(["mImgtbl","-c","Input/",outdir+"cubes.tbl"],shell=True)
    subprocess.call(["mMakeHdr",outdir+"cubes.tbl",outdir+"cubes.hdr"],shell=True)

    #second art of montage
    for fil in infils:
        tmp=fil.split(".")
        subprocess.call(["mProjectCube",outdir+"Input/"+fil,outdir+"projection/"+tmp[0]+"_proj.fits",outdir+"cubes.hdr"],shell=True)

    #final part of montage
    subprocess.call(["mImgtbl","-c",outdir+"projection/",outdir+"cubes-proj.tbl"],shell=True)
    subprocess.call(["mAddCube","-p",outdir+"projection/",outdir+"cubes-proj.tbl",outdir+"cubes.hdr",outdir+outfil],shell=True)

    #remove all montage created bits aside from output file (default is off)
    if(clean):
        subprocess.call(["rm","-rf",outdir+"Input"])
        subprocess.call(["rm","-rf",outdir+"Projection"])
        subprocess.call(["rm",outdir+"cubes.tbl",outdir+"cubes.hdr",outdir+"cubes-proj.tbl"])

    #fix the header after montage is done
    if grating == 'BL':
        fix_kcwi_cube_montage_BL(outdir+outfil)
    if grating == 'BM':
        fix_kcwi_cube_montage_BM(outdir+outfil)
        

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


    
