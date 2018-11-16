#!/usr/bin/env python

""" Impose RA/DEC offsets
"""

import pdb

def parser(options=None):

    import argparse

    parser = argparse.ArgumentParser(description='Fuss with RA/DEC for a KCWI image (v1)')
    parser.add_argument("good_coord", type=str, help="Correct coordinate, e.g. J112233.21+112233.3")
    parser.add_argument("kcwi_coord", type=str, help="KCWI current coord, e.g. J112244.12+112244.5")
    parser.add_argument("orig_cube", type=str, help="Name of original KCWI data cube file")
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
    from linetools import utils as ltu

    from kcwitools import io as kcwi_io
    from kcwitools import radec as kcwi_r

    # Coords
    coord1 = ltu.radec_to_coord(args.good_coord)
    coord2 = ltu.radec_to_coord(args.kcwi_coord)

    # Offset
    offsets, PA = kcwi_r.offsets(coord2, coord1)

    # Open original
    hdr, flux = kcwi_io.open_kcwi_cube(args.orig_cube)

    # Generate the new hdr
    ra_off = offsets[0].to('deg').value/np.cos(coord1.dec.value)
    dec_off = offsets[1].to('deg').value
    new_hdr = kcwi_r.offset_radec(hdr, ra_off, dec_off)

    # Write
    kcwi_io.write_cube(new_hdr, flux, args.new_cube)


