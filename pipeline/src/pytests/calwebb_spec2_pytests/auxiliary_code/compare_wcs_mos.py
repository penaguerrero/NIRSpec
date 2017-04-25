from __future__ import print_function, division
import numpy as np
import os

from .. import core_utils
from ..D_extract_2d import extract_2d_utils
from ..auxiliary_code import CV3_testdata_used4build7


"""
This script compares pipeline WCS info with ESA results for Multi-Object Spectroscopy (MOS) data.

slit = [shutter, ID]  --> e.g. [3, 19924]
det = 'NRS1' or 'NRS2'

"""

def compare_wcs(infile_hdr, infile_name, data_extension, esaroot):
    """
    This function does the WCS comparison from the world coordinates calculated using the
    compute_world_coordinates.py script with the ESA files. The function calls that script.

    Args:
        infile_hdr: list, header of the output fits file from the 2d_extract step
        infile_name: str, name of the output fits file from the 2d_extract step (with full path)
        data_extension: int, number of the data extension from the infile
        esaroot: str, full path of where to find all ESA intermediary products to make comparisons for the tests

    Returns:

    """

    # Run compute_world_coordinates.py in order to produce the necessary file
    os.system("compute_world_coordinates('"+infile_name+"')")

    # Read the resulting file
    cwc_fname = infile_name.replace(".fits", "_world_coordinates.fits")
    fdata = core_utils.get_filedata(cwc_fname, ext=data_extension)
    pwave = fdata[0,:,:]
    pdy = fdata[3,:,:]
    pskyx = fdata[1,:,:]
    pskyy = fdata[2,:,:]
    # get the origin of the subwindow and the grating from the extract_2d file header
    det = extract_2d_utils.find_DETECTOR(infile_hdr)
    # the science extensions are always the same in the extract_2d file
    sci_data_list, px0_list, py0_list = [], [], []
    grat = core_utils.get_keywd_val(infile_name, "GRATING", ext=0)
    filt = core_utils.get_keywd_val(infile_name, "FILTER", ext=0)
    print (grat, filt)
    sci_exts = core_utils.get_sci_extensions(infile_name)
    for se in sci_exts:
        sci_data_list.append(core_utils.get_filedata(infile_name, ext=se))
        px0_se = core_utils.get_keywd_val(infile_name, "SLTSTRT1", ext=se)+core_utils.get_keywd_val(infile_name, "SUBSTRT1", ext=0)-1
        py0_se = core_utils.get_keywd_val(infile_name, "SLTSTRT2", ext=se)+core_utils.get_keywd_val(infile_name, "SUBSTRT2", ext=0)-1
        px0_list.append(px0_se)
        py0_list.append(py0_se)
    n_p = np.shape(pwave)
    npx = n_p[0]
    npy = n_p[1]
    px = np.arange(npx)+px0_list
    py = np.arange(npy)+py0_list
    print  ("Pipeline subwindow corner pixel ID: ", px0_list, py0_list)
    # get the corresponding ESA file to the input file
    if det == "NRS1":
        file4detector = 0
    elif det == "NRS2":
        file4detector = 1
    CV3filename = CV3_testdata_used4build7.CV3_testdata_dict["MOS"][grat]["CV3filename"][file4detector]
    NID = CV3_testdata_used4build7.CV3_testdata_dict["MOS"][grat]["NID"]
    ESA_dir_name = CV3filename.split("_")[0].replace("NRS", "")+"_"+NID+"_JLAB88"
    # read in ESA data
    esafile = esaroot+"RegressionTestData_CV3_March2017_MOS/"+ESA_dir_name
    if det == "NRS1":
        eflux =
