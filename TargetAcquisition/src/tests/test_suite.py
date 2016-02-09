import unittest
import tautils
import numpy as np
from astropy.io import fits

class TestNIRSpec(unittest.TestCase):
    pass

class TestPierreF(TestNIRSpec):
    def setUp(self):
        star = '../../PierreFerruitData/TA_cutouts/postageout_star_     101 quad_       3 quad_star        1.fits'
        self.master_img = fits.getdata(star, 0)

    def test_readimage(self):
        self.assertIsInstance(tautils.readimage(self.master_img), np.ndarray)

if __name__ == '__main__':
    unittest.main()
