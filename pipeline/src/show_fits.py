import numpy as np
import os
import string
from astropy.io import fits


#  Read the fits file
path4files = '/Users/pena/Documents/PyCharmProjects/nirspec/pipeline/build6/test_data/FS_data/'
#hdulist = fits.open(path4files+'jwtest1008001_01101_00001_NRS2_uncal_mod.fits')
#hdulist = fits.open(path4files+'jwtest1008001_01101_00001_NRS2_uncal_mod_dq_init.fits')
hdulist = fits.open(path4files+'jwtest1008001_01101_00001_NRS2_uncal_mod_dq_init_saturation.fits')

# print on screen what extensions are in the file
print ('\n FILE INFORMATION: \n')
hdulist.info()
print

# get and print header
print ('\n FILE HEADER: \n')
hdr = hdulist[0].header
print repr(hdr)

# close the fits file
hdulist.close()


