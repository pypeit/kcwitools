""" Spectra routines"""
from __future__ import print_function, absolute_import, division, unicode_literals

def extract_1D_spec_5box(xcen, ycen, fil, varfil, filebase=False):
    wave, flux, hdr = open_kcwi_cube(fil)
    wave, var, hdr = open_kcwi_cube(varfil)

    # assumes sky subtracted data, and returns sqrt(var) var
    wavedim, ydim, xdim = flux.shape

    ###WATCH OUT FOR THIS:  NOW SWITCHING X WITH Y BECAUSE OF HOW CUBES READ IN and displayed
    # in the white light image

    xtmp = ycen
    ycen = xcen
    xcen = xtmp

    # Intitiate Spectrum
    spec = np.zeros(wavedim)
    err_spec = np.zeros(wavedim)

    # create an error array assuming poisson noise for flux

    for i in range(wavedim):
        slice1 = flux[i, xcen - 2:xcen + 3, ycen - 2:ycen + 3]
        varslice = var[i, xcen - 2:xcen + 3, ycen - 2:ycen + 3]
        spec[i] = np.sum(slice1)
        err_spec[i] = np.sqrt(np.abs(np.sum(varslice)))

    prihdu = fits.PrimaryHDU(spec)
    hdu_out = fits.HDUList([prihdu])
    prihdu.name = 'FLUX'

    sighdu = fits.ImageHDU(err_spec)
    sighdu.name = 'ERROR'
    hdu_out.append(sighdu)

    wvhdu = fits.ImageHDU(wave)
    wvhdu.name = 'WAVELENGTH'
    hdu_out.append(wvhdu)

    # WATCH OUT:  purposely flipping xcen and ycen here to name things according to the image diplayed
    if filebase != False:
        hdu_out.writeto(filebase + '_' + np.str(ycen) + '_' + np.str(xcen) + '_5x5.fits', overwrite=True)

    return spec, err_spec


def extract_1D_spec(xcen, ycen, fil, varfil, delx, dely):
    wave, flux, hdr = open_kcwi_cube(fil)
    wave, var, hdr = open_kcwi_cube(varfil)

    # assumes sky subtracted data, and returns sqrt(var) var
    wavedim, ydim, xdim = flux.shape

    ###WATCH OUT FOR THIS:  NOW SWITCHING X WITH Y BECAUSE OF HOW CUBES READ IN and displayed
    # in the white light image

    xtmp = ycen
    ycen = xcen
    xcen = xtmp

    delxtmp = dely
    dely = delx
    delx = delxtmp

    # Intitiate Spectrum
    spec = np.zeros(wavedim)
    err_spec = np.zeros(wavedim)

    # create an error array assuming poisson noise for flux

    for i in range(wavedim):
        slice1 = flux[i, xcen - delx:xcen + delx, ycen - dely:ycen + dely]
        varslice = var[i, xcen - delx:xcen + delx, ycen - dely:ycen + dely]
        spec[i] = np.sum(slice1)
        err_spec[i] = np.sqrt(np.abs(np.sum(varslice)))

    return spec, err_spec
