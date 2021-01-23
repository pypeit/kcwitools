""" montage routines"""
from __future__ import print_function, absolute_import, division, unicode_literals
import numpy as np

from astropy.io import fits
from astropy import units
from astropy.table import Table

from kcwitools import image as im

import subprocess
import os

def run_montage(infils,outdir='Montage/',outfil="Montage.fits",northup=False,
                trimBL=False,trimBM=False,grating='BL',clean=False,wl=3500,wh=5500):
    """ take a list of KCWI cubes and run montage on them
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
    subprocess.Popen(["mkdir",outdir])
    subprocess.Popen(["mkdir",outdir+"Input"])
    subprocess.Popen(["mkdir",outdir+"Projection"])

    #copy
    for fil in infils:
        #check to see if you are trimming, and a BM slicer
        if(trimBM):
            trim_medium_premontage(fil,wl,wh)
            a, b = fil.split(".fits")
            trim = a + '.c.fits'
            subprocess.Popen(["cp",trim,outdir+"Input"]).wait()
        elif(trimBL):
            trim_large_premontage(fil,wl,wh)
            a, b = fil.split(".fits")
            trim = a + '.c.fits'
            fix_kcwi_cube_montage_BL(trim)
            subprocess.Popen(["cp",trim,outdir+"Input"]).wait()
        else:
            subprocess.Popen(["cp",fil,outdir+"Input"]).wait()

    inputs = []
    for root, dirs,files in os.walk(outdir+'Input/'):
        for file in files:
            inputs.append(file)
                
    #first part of montage
    subprocess.Popen(["mImgtbl","-c",outdir+"Input/",outdir+"cubes.tbl"]).wait()
    subprocess.Popen(["mMakeHdr",outdir+"cubes.tbl",outdir+"cubes.hdr"]).wait()

    ##are you forcing pixels to be north up?
    if northup:
        tab=Table.read(outdir+'cubes.hdr',format='ascii',data_end=19)
        tab[17][2]='0.000000000'
        tab.add_row(['END',' ',' '])
        tab.write(outdir+'cubes.hdr',format='ascii',overwrite=True)
        
    #second part of montage
    for fil in inputs:
        tmp=fil.split(".")
        subprocess.Popen(["mProjectCube",outdir+"Input/"+fil,outdir+"Projection/"+tmp[0]+"_proj.fits",outdir+"cubes.hdr"]).wait()

    #final part of montage
    subprocess.Popen(["mImgtbl","-c",outdir+"Projection/",outdir+"cubes-proj.tbl"]).wait()
    subprocess.Popen(["mAddCube","-p",outdir+"Projection/",outdir+"cubes-proj.tbl",outdir+"cubes.hdr",outdir+outfil]).wait()

    #remove all montage created bits aside from output file (default is off)
    if(clean):
        subprocess.Popen(["rm","-rf",outdir+"Input"])
        subprocess.Popen(["rm","-rf",outdir+"Projection"])
        subprocess.Popen(["rm",outdir+"cubes.tbl",outdir+"cubes.hdr",outdir+"cubes-proj.tbl"])
        subprocess.Popen(["rm",outdir+"*.c.fits"])
        subprocess.Popen(["rm",outdir+"*.area.fits"])

        
    #fix the header after montage is done
    if grating == 'BL':
        fix_kcwi_cube_montage_BL(outdir+outfil)
    if grating == 'BM':
        fix_kcwi_cube_montage_BM(outdir+outfil)
        

###Add the CD3_3 and CDELT3 keywords to a BL grating cube
def fix_kcwi_cube_montage_BL(mosaicfil):
    hdu = fits.open(mosaicfil)
    flux = hdu['PRIMARY'].data
    hdu_hdr = hdu['PRIMARY'].header
    crval3 = hdu_hdr['CRVAL3']
    crpix3 = hdu_hdr['CRPIX3']
    
    #cdelt3 = hdu_hdr['CDELT3']
    #Make it so that there is no offset (shift zeropint and zero the offset)
    hdu_hdr['CRVAL3'] = crval3 - crpix3
    hdu_hdr['CRPIX3'] = 0.0
    hdu_hdr['CDELT3'] = 1.0
    hdu_hdr['CD3_3'] = 1.0
    

    hdu_out = fits.PrimaryHDU(flux, header=hdu_hdr)
    hdu_out.writeto(mosaicfil, overwrite=True)


###Add the CD3_3 and CDELT3 keywords to a BM grating cube
def fix_kcwi_cube_montage_BM(mosaicfil):
    hdu = fits.open(mosaicfil)
    flux = hdu['PRIMARY'].data
    hdu_hdr = hdu['PRIMARY'].header
    crval3 = hdu_hdr['CRVAL3']
    crpix3 = hdu_hdr['CRPIX3']
    #cdelt3 = hdu_hdr['CDELT3']
    hdu_hdr['CRVAL3'] = crval3 - crpix3
    hdu_hdr['CRPIX3'] = 0.0
    hdu_hdr['CDELT3'] = 0.5
    hdu_hdr['CD3_3'] = 0.5

    hdu_out = fits.PrimaryHDU(flux, header=hdu_hdr)
    hdu_out.writeto(mosaicfil, overwrite=True)

###Trim the large slicer for montage
def trim_large_premontage(infil,wl,wh):
    subprocess.Popen(["cwi_crop","-cube",infil,"-wcrop",str(wl)+":"+str(wh),"-ycrop","17:81","-xcrop","3:26"]).wait()
    
###Trim the medium slicer for montage
def trim_medium_premontage(infil,wl,wh):
    subprocess.Popen(["cwi_crop","-cube",infil,"-wcrop",str(wl)+":"+str(wh),"-ycrop","17:81","-xcrop","6:29"]).wait()
    
