#! /usr/bin/env python

from astropy.io import fits
import argparse
import string

"""
This script reads the header of a fits file and gives info of its extensions.

Example usage:
    The code works from the terminal.
    To simply see the header on-screen type:
        > python hdr_info.py blah.fits
    To see and save the header into a text file type:
        > python hdr_info.py blah.fits -s

"""

# Header
__author__ = "M. A. Pena-Guerrero"
__version__ = "1.0"

# Get arguments to run script
parser = argparse.ArgumentParser(description='')
parser.add_argument("fits_file_name",
                    action='store',
                    default=None,
                    help='Name of fits file, i.e. blah.fits')
parser.add_argument("-s",
                    dest="save_txt",
                    action='store_true',
                    default=False,
                    help='Save the text file with the header information.')
args = parser.parse_args()

# Set the variables
fits_file_name = args.fits_file_name
save_txt = args.save_txt

#  Read the fits file
hdulist = fits.open(fits_file_name)

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

# save the info into a text file
if save_txt:
    # set the name of the text file
    text_file_name = string.replace(fits_file_name, '.fits', '_header.txt')
    tf = open(text_file_name, 'w')
    tf.write(repr(hdr))
    tf.close()
