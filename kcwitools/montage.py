""" montage routines"""
from __future__ import print_function, absolute_import, division, unicode_literals

import numpy as np
from astropy.io import fits
from astropy import units
from astropy.table import Table
from kcwitools import image as im
import subprocess
import os
import shutil
from MontagePy.main import *
from MontagePy.archive import *

from cwitools.scripts import cwi_crop 


#if it's from the python pipeline, create Flux and Variance arrays that are
#separate

def pyfil_strip(infil,flux=True):   
    pyfil = fits.open(infil)  #read
  
    fhdr = pyfil[0].header    #pull flux header and data
    fdata = pyfil[0].data

    vhdr = pyfil[0].header    #pull variance header and data
    #RB hack: the *icubes.fits files don't have a valid variance header. Using flux header instead.
    vdata = pyfil[2].data
    
    if(flux):
        a, b = infil.split(".fits")
        outfil = a+'_flux.fits'
        fits.writeto(outfil,fdata,fhdr)
    else:      
        a, b = infil.split(".fits")
        outfil = a+'_var.fits'
        fits.writeto(outfil,vdata,vhdr)
    return

def fix_kcwi_cube_pre_montage(infil):  #to make sure cdelt3 is populated
    hdu = fits.open(infil)
    flux = hdu['PRIMARY'].data
    hdu_hdr = hdu['PRIMARY'].header
    crval3 = hdu_hdr['CRVAL3']
    crpix3 = hdu_hdr['CRPIX3']
    cd3_3 = hdu_hdr['CD3_3']

    #Make it so that there is no offset (shift zeropint and zero the offset)
    hdu_hdr['NAXIS'] = 3
    hdu_hdr['CRVAL3'] = crval3 - crpix3
    hdu_hdr['CRPIX3'] = 0.0
    hdu_hdr['CDELT3'] = cd3_3
    

    hdu_out = fits.PrimaryHDU(flux, header=hdu_hdr)
    hdu_out.writeto(infil, overwrite=True)

###Trim the large slicer for montage
def trim_premontage(infil,wc1=3500,wc2=5500,xc1=3,xc2=26,yc1=15,
                    yc2=78):
    cwi_crop(infil,xcrop=(xc1,xc2),ycrop=(yc1,yc2),wcrop=(wc1,wc2))

###Trim the large slicer for montage
def trim_premontage_medium(infil,wc1=3500,wc2=5500,xc1=6,xc2=29,yc1=17,
                         yc2=81):
    cwi_crop(infil,xcrop=(xc1,xc2),ycrop=(yc1,yc2),wcrop=(wc1,wc2))

#now for the main routine
def run_montage(infils,outfil="Montage.fits",northup=False,trim=True,
                pycube=True,flux=True,wc1=3500,wc2=5500,xc1=3,xc2=26,
                yc1=15,yc2=78):
    
    """ take a list of KCWI cubes and run montage on them
    Args:
    ----------
    infils: list
      list of kcwi cubes to montage
    outfil: string
      name of montage output file, defaults to Montage.fits'
    northup:  Force North = up in the reprojection (default = no)
    trim:  Remove junk from the cube edges (default = yes, and assumes
           large slicer)
    pycube:  Set True if file from the new python pipeline, false if from the
             IDL pipleline (defaults to true)
    flux: Set true for a flux file, false if a variance file
    """

    #start by setting up the directories, remove the Montage directory if it exists.
    workdir = "Montage"

    try:
        home
    except:
        home = os.getcwd()
    print("Startup folder: " + home)

    #try:
    #    shutil.rmtree(workdir)
    #except:
    #    print("                Can't delete work tree; probably doesn't exist yet", flush=True)


    print("Work directory: " + workdir, flush=True)
    #--------------------
    #make the directories in the nested structure
    #Check if workdir exists if not do the following commands

    #Check if workdir exists if not create it

    if os.path.isdir(workdir)==False:
        os.makedirs(workdir)  
        os.chdir(workdir)
        os.makedirs("raw")
        os.makedirs("projected")

    #need to come back, to see the files
    os.chdir(home)
    #--------------------
    
    #now fix,strip (if needed), and trim and fill up the raw montage directories
    for fil in infils:
        print(fil)
        #fix_kcwi_cube_pre_montage(fil)
        #check to see if a python pipeline cube
        if (pycube):
            #flux or variance array?
            if (flux):
                pyfil_strip(fil)
                a, b = fil.split(".fits")

                if(trim):
                    trim_premontage(a+'_flux.fits',wc1=wc1,wc2=wc2,
                                        xc1=6,xc2=29,yc1=17,yc2=81)
                    trim = a + '_flux.c.fits'
                    fix_kcwi_cube_pre_montage(trim)
                    shutil.move(trim,workdir+'/raw/')
                else:
                    fix_kcwi_cube_pre_montage(a+'_flux.fits')
                    shutil.move(a+'_flux.fits',workdir+'/raw/')
            else:
                pyfil_strip(fil,flux=False)
                a, b = fil.split(".fits")

                if(trim):
                    trim_premontage(a+'_var.fits',wc1=wc1,wc2=wc2,
                                        xc1=6,xc2=29,yc1=17,yc2=81)
                    trim = a + '_var.c.fits'
                    fix_kcwi_cube_pre_montage(trim)
                    shutil.move(trim,workdir+'/raw/')
                else:
                    fix_kcwi_cube_pre_montage(a+'_var.fits')
                    shutil.move(a+'_var.fits',workdir+'/raw/')
        else:
           if(trim):
               trim_premontage(fil,wc1=wc1,wc2=wc2,
                               xc1=6,xc2=29,yc1=17,yc2=81)
               a, b = fil.split(".fits")
               trim = a + '.c.fits'
               fix_kcwi_cube_pre_montage(trim)
               shutil.move(trim,workdir+'/raw/')
           else:
               fix_kcwi_cube_pre_montage(fil)
               shutil.move(fil,workdir+'/raw/')
               
    #now get to work on the montage routines
    os.chdir(workdir)

    #make an image table
    rtn = mImgtbl("raw", "rimages.tbl")
    print("mImgtbl (raw):    " + str(rtn), flush=True)

    #make a notional cube header
    rtn = mMakeHdr("rimages.tbl","rcubes.hdr")
    print(rtn)

    #make a list of input files
    inputs = []
    for root, dirs,files in os.walk('raw/'):
        for file in files:
            inputs.append(file)

    #do the projection
    for fil in inputs:
        a,b = fil.split(".fits")
        
        rtn = mProjectCube("raw/"+fil, "projected/"+a+"_proj.fits",
                            "rcubes.hdr")
        print("mProjectCube:           " + str(rtn), flush=True)

    #make a table of projected images
    mImgtbl("projected", "pimages.tbl")

    #add up the cubes!
    rtn = mAddCube("projected", "pimages.tbl", "rcubes.hdr", outfil)
    print("mAddCube:    " + str(rtn), flush=True)

    os.chdir(home)
    print("all done!!!")


