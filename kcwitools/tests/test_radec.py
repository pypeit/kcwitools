# Module to run tests on radec
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import numpy as np
import pytest

from astropy import units
from astropy.io import fits

from linetools import utils as ltu

from kcwitools import radec

coord1 = ltu.radec_to_coord('J214425.25-405400.81')
coord2 = ltu.radec_to_coord('J214425.31-405402.6')

def test_offsets():
    # Run
    offsets, PA = radec.offsets(coord1, coord2)
    # Test
    assert np.isclose(offsets[0].value, 0.68026069)
    assert PA.unit == units.deg

def test_offset_radec():
    # Dummy header
    hdr = fits.Header()
    hdr['CRVAL1'] = 0.
    hdr['CRVAL2'] = 0.
    #
    offsets, PA = radec.offsets(coord1, coord2)
    # Run
    new_hdr = radec.offset_radec(hdr, offsets[0], offsets[1], coord2)
    # Test
    assert np.isclose(new_hdr['CRVAL1'], 0.000250)

