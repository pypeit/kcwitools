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

    #in case the user forgets the trailing /
    outdir = outdir+'/'
    
    #create the driectories
    subprocess.Popen(["mkdir",outdir+"Input"])
    subprocess.Popen(["mkdir",outdir+"Projection"])

    #copy
    for fil in infils:
        subprocess.Popen(["cp",fil,outdir+"Input"]).wait()

    inputs = []
    for files in os.walk(outdir+'Input/'):
        for file in files:
            if file.endswith('.fits'):
                inputs.append(file)
                
    #first part of montage
    subprocess.Popen(["mImgtbl","-c",outdir+"Input/",outdir+"cubes.tbl"]).wait()
    subprocess.Popen(["mMakeHdr",outdir+"cubes.tbl",outdir+"cubes.hdr"]).wait()

    #second art of montage
    for fil in inputs:
        tmp=fil.split(".")
        subprocess.Popen(["mProjectCube",outdir+"Input/"+fil,outdir+"projection/"+tmp[0]+"_proj.fits",outdir+"cubes.hdr"]).wait()

    #final part of montage
    subprocess.Popen(["mImgtbl","-c",outdir+"projection/",outdir+"cubes-proj.tbl"]).wait()
    subprocess.Popen(["mAddCube","-p",outdir+"projection/",outdir+"cubes-proj.tbl",outdir+"cubes.hdr",outdir+outfil]).wait()

    #remove all montage created bits aside from output file (default is off)
    if(clean):
        subprocess.Popen(["rm","-rf",outdir+"Input"])
        subprocess.Popen(["rm","-rf",outdir+"Projection"])
        subprocess.Popen(["rm",outdir+"cubes.tbl",outdir+"cubes.hdr",outdir+"cubes-proj.tbl"])

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


    
