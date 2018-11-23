#!/usr/bin/env python

""" Loads and plots a requested spectrum
"""

import pdb

def parser(options=None):

    import argparse

    parser = argparse.ArgumentParser(description='Make/Show/Write a white light image (v1)')
    parser.add_argument("file", type=str, help="Data file (cube)")
    #parser.add_argument("--meta", default=True, help="Show meta data? [default: True]", action="store_true")
    #parser.add_argument("-s", "--select", default=0, type=int, help="Index of spectrum to plot (when multiple exist)")
    parser.add_argument("--outfile", help="Write the image to this data file")
    parser.add_argument("--wcs", default=False, help="If true plot image in WCS coordinate", action='store_true')


    if options is None:
        args = parser.parse_args()
    else:
        args = parser.parse_args(options)
    return args


def main(args):
    """ Run
    """
    from kcwitools import io as kcwi_io
    from kcwitools import utils as kcwi_u
    from kcwitools import image as kcwi_img
    from kcwitools import plot as kcwi_p

    # Load
    hdr, flux = kcwi_io.open_kcwi_cube(args.file)

    # White light
    whiteim = kcwi_img.build_whitelight(hdr, flux, outfile=args.outfile)

    # Show
    if args.wcs is not False:
        head = hdr
        kcwi_p.show_whitelight(whiteim, header=head)
    else:
        head = None
        # Show
        kcwi_p.show_whitelight(whiteim, header=head)

