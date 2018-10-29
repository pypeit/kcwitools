""" Routine to extract header information from KCWI datacube"""
from __future__ import print_function, absolute_import, division, unicode_literals
from astropy.io import fits



def defs(header):
    """
    Args: 
        header  :  KCWI datacube header

    Returns: 
        out  :  a dictionary with header information


    """
    # CHECK CHECK CHECK
    #THESE VALUES NEED TO BE DOUBLE CHECKED. WHICH ONE SHOULD I USE?
    out= {}
    out['gain']=header['GAIN2']
    out['rdnoise']=3.5   #  CHECK CHECK CHECK

    return out
