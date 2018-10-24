#!/usr/bin/env python

""" Loads and plots a requested spectrum
"""

import pdb

def parser(options=None):

    import argparse

    parser = argparse.ArgumentParser(description='Make/Show/Write a narrow band image (v1)')
    parser.add_argument("file", type=str, help="Data file (cube)")
    parser.add_argument("line", type=float, help="Line (wavelength; Ang) to center on")
    parser.add_argument("--sub_off", default=False, help="Subtract off an off-band image?", action="store_true")
    parser.add_argument("--z", default=0., type=float, help="Redshift of the emission")
    parser.add_argument("--outfile", help="Write the narrow band image to this data file")
    parser.add_argument("--del_wave", default=2., type=float, help="Width of Narrow Band in rest-frame")


    if options is None:
        args = parser.parse_args()
    else:
        args = parser.parse_args(options)
    return args


def main(args):
    """ Run
    """
    from kcwitools import io as kcwi_io
    from kcwitools import image as kcwi_img
    from kcwitools import plot as kcwi_p

    # Load
    hdr, flux = kcwi_io.open_kcwi_cube(args.file)

    # Narrow band
    nb = kcwi_img.build_narrowband(hdr, flux, args.line, z=args.z,del_wave=args.del_wave,
                                   sub_offimage=args.sub_off,
                                   outfile=args.outfile)

    # Show
    kcwi_p.show_narrowband(nb)


