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
        fits_file_name: full path with name of the header text file
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


def get_filedata(fits_file_name, extension=None):
    """
    This function gets the data from the science extension of the fits file provided.
    Args:
        fits_file_name: name of the fits file used as input for the pipeline
        extension: integer, if no number is given as input, the default is extension 1

    Returns:
        fdata: the data extension of the file
    """
    if extension is None:
        ext = 1
    else:
        ext = extension
    fdata = fits.getdata(fits_file_name, ext)
    return fdata


def get_keywd_val(fits_file_name, keywd, ext=0):
    """
    This function obtains the value corresponding to the given keyword.
    Args:
        fits_file_name: name of the fits file used as input for the pipeline
        keywd: keyword for which to obtain the value
        ext: extension in which the kwyword lives, by default it is set to the
             primary extension.

    Returns:
        keywd_val: the value corresponding to the inputed keyword
    """
    keywd_val = fits.getval(fits_file_name, keywd, ext)
    return keywd_val


def get_sci_extensions(fits_file_name):
    """
    This functions obtains all the science extensions in the given file
    Args:
        fits_file_name: name of the fits file of interest

    Returns:
        sci_list: list of the numbers of the science extensions
    """
    hdulist = fits.open(fits_file_name)
    sci_list = []
    for ext, hdu in enumerate(hdulist):
        if hdu.name == "SCI":
            sci_list.append(ext)
    return sci_list
