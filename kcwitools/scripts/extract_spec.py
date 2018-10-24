#!/usr/bin/env python

""" Extracts a 1D spectrum
"""

import pdb

def parser(options=None):

    import argparse

    parser = argparse.ArgumentParser(description='Make/Show/Write a narrow band image (v1)')
    parser.add_argument("flux", type=str, help="Data flux file (cube)")
    parser.add_argument("var", type=str, help="Data variance file (cube)")
    parser.add_argument("x", type=int, help="x centroid in the cube; integer")
    parser.add_argument("y", type=int, help="y centroid in the cube; integer")
    parser.add_argument("--square", type=int, help="Use a square aperture with input size (odd integer)")
    parser.add_argument("--outfile", help="Write the narrow band image to this data file")

    if options is None:
        args = parser.parse_args()
    else:
        args = parser.parse_args(options)
    return args


def main(args):
    """ Run
    """
    from kcwitools import io as kcwi_io
    from kcwitools import spec as kcwi_s
    from kcwitools import utils as kcwi_u

    # Load
    hdr, flux = kcwi_io.open_kcwi_cube(args.flux)
    _, var = kcwi_io.open_kcwi_cube(args.var)
    wave = kcwi_u.build_wave(hdr)

    # Spectrum

    # Square?
    xspec = None
    if args.square is not None:
        xspec = kcwi_s.extract_square(args.x, args.y, wave, flux, var, squaresize=args.square,
                                      outfile=args.outfile)
    # Show
    if xspec is not None:
        xspec.plot(xspec=True)
