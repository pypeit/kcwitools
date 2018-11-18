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
    #Check if keyword exists else put them to default values
    if 'GAIN2' in header.keys():
        out['gain']=header['GAIN2']
    else:
        out['gain']=1.
    if 'BIASRN2' in header.keys():
        out['rdnoise']=header['BIASRN2']   #  CHECK CHECK CHECK
        out['BIASSUB']=header['BIASSUB']  # Has bias been subtracted or not T/F
    else:
        out['rdnoise']=0.   #  CHECK CHECK CHECK
        out['BIASSUB']='T'  # Has bias been subtracted or not T/F


    return out
