from __future__ import print_function, division
from astropy.io import fits

'''
This script contains functions frequently used in the test suite.
'''

def read_hdrfits(fits_file_name, info=False, show_hdr=False):
    '''
    This function reads the header fits file and returns a dictionary of the keywords with
    corresponding values. Keywords will be stored in the order they are read.
    Args:
        hdr_txt_file: full path with name of the header text file
        info: if True the function will show the contents and shapes of the file
        show_hdr: if True the function will print the header of the file

    Returns:
        hdrl: The header of the fits file
    '''
    #  Read the fits file
    hdulist = fits.open(fits_file_name)
    # print on screen what extensions are in the file
    if info:
        print ('\n FILE INFORMATION: \n')
        hdulist.info()
    # get and print header
    hdrl = hdulist[0].header
    if show_hdr:
        print ('\n FILE HEADER: \n')
        print (repr(hdrl))
    # close the fits file
    hdulist.close()
    return hdrl


