from __future__ import print_function, division
import os
import collections
import numpy as np
import string
from astropy.io import fits

'''
This script contains functions frequently used in the test suite.
'''

def read_hdrfits(fits_file_name):
    '''
    This function reads the header fits file and returns a dictionary of the keywords with
    corresponding values. Keywords will be stored in the order they are read.
    Args:
        hdr_txt_file: full path with name of the header text file

    Returns:
        A dictionary of keywords with corresponding values
    '''
    #  Read the fits file
    hdulist = fits.open(fits_file_name)
    # print on screen what extensions are in the file
    #print ('\n FILE INFORMATION: \n')
    hdulist.info()
    # get and print header
    #print ('\n FILE HEADER: \n')
    hdr = hdulist[0].header
    #print (repr(hdr))
    # close the fits file
    hdulist.close()
    # set the name of the text file and save the header
    text_file_name = string.replace(fits_file_name, '.fits', '_header.txt')
    tf = open(text_file_name, 'w')
    tf.write(repr(hdr))
    tf.close()
    # read the text file
    keywd_dict = read_hdrtxt(text_file_name)
    # remove the text file
    os.system("rm "+text_file_name)
    return keywd_dict


def read_hdrtxt(hdr_txt_file):
    '''
    This function reads the header text file and returns a dictionary of the keywords with
    corresponding values. Keywords will be stored in the order they are read.
    Args:
        hdr_txt_file: full path with name of the header text file

    Returns:
        A dictionary of keywords with corresponding values
    '''
    keywd_dict = collections.OrderedDict()
    with open(hdr_txt_file, 'r') as htf:
        for line in htf.readlines():   # identify keywords by lines containing a =
            if '=' in line:
                line_list = line.split('=')
                keywd = line_list[0].split()[0]   # remove the white spaces from the keyword
                keywd_val = line_list[1].split()[0]  # remove the white spaces from the keyword value
                if "'" in keywd_val:
                    keywd_val = keywd_val.replace("'", "")   # remove the extra '
                keywd_dict[keywd] = keywd_val   # add dictionary entry
    return keywd_dict


