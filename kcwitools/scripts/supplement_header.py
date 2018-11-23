#!/usr/bin/env python

""" Supplement a KCWI data cube header
"""

import pdb

def parser(options=None):

    import argparse

    parser = argparse.ArgumentParser(description='Fuss with RA/DEC for a KCWI image (v1)')
    parser.add_argument("orig_cube", type=str, help="Name of original KCWI data cube file")
    parser.add_argument("raw_file", type=str, help="Name of KCWI raw data file for the header")
    parser.add_argument("new_cube", type=str, help="Name for new KCWI data cube file (can be same)")
    #parser.add_argument("--square", type=int, help="Use a square aperture with input size (odd integer)")
    #parser.add_argument("--outfile", help="Write the narrow band image to this data file")

    if options is None:
        args = parser.parse_args()
    else:
        args = parser.parse_args(options)
    return args


def main(args):
    """ Run
    """
    import numpy as np
    from astropy.io import fits

    from kcwitools import io as kcwi_io

    # Open original data cube
    hdr, flux = kcwi_io.open_kcwi_cube(args.orig_cube)

    # Open raw file
    hdul = fits.open(args.raw_file)
    head_orig = hdul[0].header

    # Keys
    new_hdr = hdr.copy()
    for key in ['MJD', 'BGRATNAM', 'FRAMENO', 'ROTPOSN', 'AIRMASS', 'DATE-OBS']:
        new_hdr[key] = head_orig[key]

    # Write
    kcwi_io.write_cube(new_hdr, flux, args.new_cube)


