from __future__ import print_function
from astropy.io import fits
from configobj import ConfigObj
import argparse
import string
import os
import matplotlib.pyplot as plt
import numpy as np


"""
This script reads the header of a fits file and gives info of its extensions, as well as cheching the saturation
step of the pipeline. The script is based on the satcheck2.pro provided by Cheryl Pavlovsky.

Example usage:
    In a terminal type the following (make sure you call the code from the right directory):
    > python satcheck.py scifile_name.fits input_file4saturation_step.fits

    If the header is to be saved into a text file and/or plot the data then type -s and/or -p,respectively:
    > python modify_fits.py file_name.fits -s -p

"""

# Header
__author__ = "M. A. Pena-Guerrero"
__version__ = "1.0"


# Type the name of the person modifying the fits files
name_editor = 'Maria A. Pena-Guerrero'

# Get arguments to run script
parser = argparse.ArgumentParser(description='')
parser.add_argument("scifile_name",
                    action='store',
                    default=None,
                    help='Name of fits file, i.e. blah.fits')
#parser.add_argument("ref_file_name",
#                    action='store',
#                    default=None,
#                    help='Name of reference file, i.e. blah.fits')
parser.add_argument("input_file4saturation_step",
                    action='store',
                    default=None,
                    help='Name of input file for the saturation step, i.e. blah.fits')
parser.add_argument("-p",
                    dest="make_plot",
                    action='store_true',
                    default=False,
                    help='Plot the first data column versus the second.')
parser.add_argument("-s",
                    dest="save_txt",
                    action='store_true',
                    default=False,
                    help='Save the text file with the header information.')
args = parser.parse_args()

# Set the variables
scifile_name = args.scifile_name
#ref_file_name = args.ref_file_name
input_file4saturation_step = args.input_file4saturation_step
make_plot = args.make_plot
save_txt = args.save_txt

#  Read the science fits files
hdulist = fits.open(scifile_name)
# print on screen what extensions are in the file
print ('\n *** SCIENCE FILE INFORMATION: ')
hdulist.info()
# get and print header from file
sci_hdr = hdulist[0].header
#print ('\n *** SCIENCE FILE HEADER: \n')
#print repr(sci_hdr)
science = hdulist[1].data
pixdq = hdulist[2].data
grpdq = hdulist[3].data
# close the fits file
hdulist.close()

# Obtain values from keywords
scix = sci_hdr["SUBSIZE1"]
sciy = sci_hdr["SUBSIZE2"]
ngroups = sci_hdr["NGROUPS"]
tframe = sci_hdr["TFRAME"]
ref_file_name = sci_hdr["R_SATURA"].split("//")[1]   # Get the reference file used
subarray = sci_hdr["SUBARRAY"].split()[0]
fullf = "FULL"
allslitsf = "ALLSLITS"
print (" saturation map used:  ", ref_file_name)
print (" saturation corr done:  ", sci_hdr["S_SATURA"])

# Get the infro from the reference file used
pileline_path = string.split(os.getcwd(), sep='pipeline')[0]
build6_path = os.path.join(pileline_path, 'pipeline/build6')
test_data_path = os.path.join(build6_path, 'test_data')
FS_data_path = os.path.join(test_data_path, 'FS_data')
reference_files_path = os.path.join(build6_path, "reference_files")
ref_hdulist = fits.open(os.path.join(reference_files_path, ref_file_name))
print ("reference file name = ", os.path.join(reference_files_path, ref_file_name))
# print on screen what extensions are in the file
print ('\n *** REFERENCE FILE INFORMATION: ')
ref_hdulist.info()
# get and print header from file
ref_hdr = ref_hdulist[0].header
#print ('\n ***  REFERENCE  FILE HEADER: \n')
#print repr(ref_hdr)
satref = ref_hdulist[1].data
refdq = ref_hdulist[2].data
# close the fits file
ref_hdulist.close()

input_hdulist = fits.open(input_file4saturation_step)
# print on screen what extensions are in the file
print ('\n *** INPUT FILE FOR SATURATIOM STEP INFORMATION: ')
input_hdulist.info()
# get and print header from file
input_hdr = input_hdulist[0].header
#print ('\n ***  INPUT  FILE HEADER: \n')
#print repr(input_hdr)
incube = input_hdulist[1].data
indq = input_hdulist[2].data
# close the fits file
input_hdulist.close()

# Slice the data according to keyword
print('shape satref before slicing: ', np.shape(satref))
print('shape refdq before slicing: ', np.shape(refdq))
#print(satref)
# if subarray == fullf don't do anything
if subarray == allslitsf:
    satref = satref[896:1152, :]
    refdq = refdq[896:1152, :]
else:
    print ("Unknown array size.")
print('shape satref after slicing: ', np.shape(satref))
print('shape refdq after slicing: ', np.shape(refdq))

# Find and count where are the flagged pixels
print ('shape of science ', np.shape(science))
print ('shape of input file for saturation step ', np.shape(incube))
print ('ngroups = ', ngroups)
for frame_number in range(0, ngroups-1):
    diff = science[:, frame_number, :, :] - incube[:, frame_number, :, :]
    dd = np.where(diff != 0.0)
    count = len(dd[-1])
    print ("Sci pixels in group ", frame_number+1, " that changed value: ", count)

# Find where pixel values passed the threshold
print ('shape of grpdq: ', np.shape(grpdq))
badpix = 0
counter = 0
for hh in range(0, scix):
    for ii in range(0, sciy):
        threshold = satref[ii, hh]
        val = np.where(science[:, :, ii, hh] < threshold)
        c = len(val[-1])
        if c == 0:
            badpix += badpix
        else:
            flagged1 = np.where(grpdq[:, 0:c, ii, hh] >= 2.0)
            flagged2 = np.where(grpdq[:, c:ngroups, ii, hh] != 2.0)
            count = len(flagged1[-1])+len(flagged2[-1])
            counter = counter + count

print (' # of pixels/groups with inacurate flags: ', counter)
print (' # of bad (all groups saturated) pixels: ', badpix)

# Check for flag propagation from DQ to the sciene file PIXELDQ
print ('Size of refdq: ', np.shape(refdq))
print ('Size of indq: ', np.shape(indq))
newdq = refdq
diff = newdq - pixdq
flags = np.where(diff != 0.0)
flagcount = len(diff[-1])
print (' # of different flags in the pixeldq: ', flagcount)


print ('\n * Script satcheck.py finished * \n')