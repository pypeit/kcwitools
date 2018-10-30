""" montage routines"""
from __future__ import print_function, absolute_import, division, unicode_literals
import numpy as np

from astropy import units


def run_montage(infils,outfil="Montage.fits",clean=True):
    """ take a list of (ideally trimmed) KCWI cubes and run montage on them
    Args:
    ----------
    infils: list
      list of kcwi cubes to montage
    outfil: string
      name of montage output file, defaults to Montage.fits'

    
    """

    subprocess.call(["mkdir","Input"])
    subprocess.call(["mkdir","Projection"])

    for fil in infils:
        subprocess.call(["cp",fil,"Input"])

    subprocess.call(["mImgtbl","-c","Input/","cubes.tbl"])

    if(clean):
        subprocess.call(["rm","-rf","Input"])
        subprocess.call(["rm","-rf","Projection"])


#mImgtbl -c Input/ cubes.tbl
#mMakeHdr cubes.tbl cubes.hdr
#mProjectCube Input/kb180905_00062_icubes_trimmed.fits projection/kb180905_00062_icubes_trimmed_proj.fits cubes.hdr
#mProjectCube Input/kb180905_00063_icubes_trimmed.fits projection/kb180905_00063_icubes_trimmed_proj.fits cubes.hdr
#mProjectCube Input/kb180905_00064_icubes_trimmed.fits projection/kb180905_00064_icubes_trimmed_proj.fits cubes.hdr
#mProjectCube Input/kb180905_00065_icubes_trimmed.fits projection/kb180905_00065_icubes_trimmed_proj.fits cubes.hdr
#mImgtbl -c projection/ cubes-proj.tbl
#mAddCube -p projection/ cubes-proj.tbl cubes.hdr J2222_0918.fits


###Add the CD3_3 and CDELT3 keywords to a cube made by montage
def fix_kcwi_cube_montage(mosaicfil):
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
def fix_kcwi_cube_montage_2(mosaicfil):
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


    
