from __future__ import print_function, division
import numpy as np
import os
from astropy.io import fits


"""
This script contains the auxiliary functions that the wcs FS, MOS, and IFU WCS scripts use.
"""


def find_nearest(arr, value):
    '''
    This function gives the content and the index in the array of the number that is closest to
    the value given.
    :param arr = 1-D numpy array
    :param value = float or integer
    :return: The array element closest to value and its index
    '''
    idx=(np.abs(arr-value)).argmin()
    return arr[idx], idx


def get_sci_extensions(fits_file_name):
    """
    This functions obtains all the science extensions in the given file
    Args:
        fits_file_name: name of the fits file of interest

    Returns:
        sci_list: list of the numbers of the science extensions
    """
    hdulist = fits.open(fits_file_name)
    hdulist.info()
    sci_list = []
    for ext, hdu in enumerate(hdulist):
        if hdu.name == "SCI":
            sci_list.append(ext)
    return sci_list


def do_idl_match(arrA, arrB):
    """
    This function does the same that the IDL match function does. It finds the elements common
    in both arrays and it returns two arrays with the index of those elements in each array.
    (The arrays do not need to be the same length)
    Args:
        arrA: numpy array
        arrB: numpy array

    Returns:
        subA: numpy array of index of arrA from elements also present in arrB
        subB: numpy array of index of arrB from elements also present in arrA
    """
    # Find the index corresponding to the intersection elements and return them as arrays
    subA, subB = [], []
    for i, ai in enumerate(arrA):
        if ai in arrB:
            subA.append(i)
    for i, bi in enumerate(arrB):
        if bi in arrA:
            subB.append(i)
    return np.array(subA), np.array(subB)


def do_idl_rebin(a, *args):
    '''
    * This function was copied from Example 2 in http://scipy-cookbook.readthedocs.io/items/Rebinning.html

    This acts identically to IDL's rebin command where all values in the original array are summed
    and divided amongst the entries in the new array. As in IDL, the new shape must be a factor of
    the old one. The ugly 'evList trick' builds and executes a python command of the form
    a.reshape(args[0],factor[0],).sum(1)/factor[0]
    a.reshape(args[0],factor[0],args[1],factor[1],).sum(1).sum(2)/factor[0]/factor[1]
    etc. This general form is extended to cover the number of required dimensions.

    This function rebins ndarray data into a smaller ndarray of the same rank whose dimensions
    are factors of the original dimensions. eg. An array with 6 columns and 4 rows
    can be reduced to have 6,3,2 or 1 columns and 4,2 or 1 rows.
    example usages:
     a=rand(6,4); b=rebin(a,3,2)
     a=rand(6); b=rebin(a,2)
    '''
    shape = a.shape
    lenShape = len(shape)
    factor = np.asarray(shape)/np.asarray(args)
    evList = ['a.reshape('] + \
             ['args[%d],factor[%d],'%(i,i) for i in range(lenShape)] + \
             [')'] + ['.sum(%d)'%(i+1) for i in range(lenShape)] + \
             ['/factor[%d]'%i for i in range(lenShape)]
    #print (''.join(evList))
    return eval(''.join(evList))


def get_esafile(esa_files_path, det, mode, configuration):
    """
    This function gets the ESA file corresponding to the input given.
    Args:
        esa_files_path: str, path where to find all the files with single configuration names
        det: string, detector e.g "NRS1"
        mode: string, either 'MOS', 'FS', or 'IFU'
        configuration: list of 2 strings, grating and filter of configuration

    Returns:
        esafile: str, full path of the ESA file corresponding to input given
    """
    # get all the fits into a searchable list
    fits_files_list = []
    for root, dirs, files in os.walk(os.path.join(esa_files_path, mode)):
        for file in files:
            if file.endswith('.fits'):
                fits_files_list.append(file)
    # Find the file of interest
    esafile = 'configuration not found'
    grat, filt = configuration
    for fitsfile in fits_files_list:
        if (grat in fitsfile) and (filt in fitsfile) and (det in fitsfile):
            esafile = fitsfile
    return esafile
