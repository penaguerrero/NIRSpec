#! /usr/bin/env python

# This script converts the University of Arizona NIRCam reference files into SSB format.

import sys, os,re,math
import pyfits,optparse,scipy
import numpy as np
from jwst_lib import models

from jwst_lib import stpipe
jwstver=stpipe.__svn_revision__
print "*********Hello, you are currently using the JWST pipeline version****************:", jwstver


# read in a MIRI reference file as a data model,
# make a few modifications and output it back in the correct format

# Get the input and output file names from the user
input_name = raw_input ('Enter the input file name: ')
output_name = raw_input ('Enter the output file name: ')

# Open the input FITS file and load the header and data arrays
# MIRI straylight mask has [0] header and [MASK] extension
#input=models.StrayLightModel(input_name)
# read in a science file and output it with NaNs, not a mask file
input=models.ImageModel(input_name)

# make changes to data
# add a box of NaN values to check if pipeline handles them correctly.
input.data[310:330,380:400] = float('NaN')


# Write the output file
input.save(output_name)