""" utility routines"""
from __future__ import print_function, absolute_import, division, unicode_literals

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






