"""
py.test module for unit testing the extract_2d step.
"""

import pytest
import os
import numpy as np
from jwst.extract_2d.extract_2d_step import Extract2dStep

from .. import core_utils
from . import extract_2d_utils


# Set up the fixtures needed for all of the tests, i.e. open up all of the FITS files

# Default names of pipeline output files
prev_file_name = "_assign_wcs_subtract_images"
out_file_name = "_imprint"

# This is the default name the pipeline will give an output file that did not have association (ASN) files
# for the background and imprint subtraction steps
no_ASN_files_in = "_assign_wcs"
no_ASN_files_out = "_extract_2d"

# fixture to read the config file
@pytest.fixture(scope="module")
def input_hdul(request, config):
    initiate_calwebb_spc2 = "calwebb_spec2_input_file"
    if config.has_option(initiate_calwebb_spc2, "input_file"):
        working_directory = config.get(initiate_calwebb_spc2, "working_directory")
        initial_input_file = config.get(initiate_calwebb_spc2, "input_file")
        initial_input_file_fullpath = os.path.join(working_directory, initial_input_file)
        next_input_file = initial_input_file_fullpath.replace(".fits", prev_file_name+".fits")
        if not os.path.isfile(next_input_file):
            next_input_file = initial_input_file_fullpath.replace(".fits", no_ASN_files_in+".fits")
        hdul = core_utils.read_hdrfits(next_input_file, info=True, show_hdr=True)
        return hdul, next_input_file
    else:
        pytest.skip("extract_2d needs an input_file")


# fixture to read the output file header
@pytest.fixture(scope="module")
def output_hdul(request, config):
    step = "extract_2d"
    initiate_calwebb_spc2 = "calwebb_spec2_input_file"
    working_directory = config.get(initiate_calwebb_spc2, "working_directory")
    initial_input_file = config.get(initiate_calwebb_spc2, "input_file")
    initial_input_file_fullpath = os.path.join(working_directory, initial_input_file)
    next_input_file = initial_input_file_fullpath.replace(".fits", prev_file_name+".fits")
    if os.path.isfile(next_input_file):
        output_file = next_input_file.replace(".fits", out_file_name+".fits")
    else:
        next_input_file = initial_input_file_fullpath.replace(".fits", no_ASN_files_in+".fits")
        output_file = next_input_file.replace(".fits", no_ASN_files_out+".fits")
    stp = Extract2dStep()
    run_calwebb_spec2 = config.get("run_calwebb_spec2_in_full", "run_calwebb_spec2")
    # if run_calwebb_spec2 is True calwebb_spec2 will be called, else individual steps will be ran
    if not run_calwebb_spec2:
        print ("Will run "+step+" step ...")
        #result = stp.call(next_input_file)
        #result.save(output_file)
    if config.has_option("steps", step):
        hdul = core_utils.read_hdrfits(output_file, info=True, show_hdr=True)
        return hdul, next_input_file
    else:
        pytest.skip("needs extract_2d output_file")


### THESE FUNCTIONS ARE TO VALIDATE BOTH THE WCS AND THE 2D_EXTRACT STEPS

# fixture to read the data of the fits file
@pytest.fixture(scope="module")
def get_fdata(output_hdul):
    if extract_2d_utils.check_FS_true(output_hdul):
        # Find what slit the data corresponds to
        ext, slit = extract_2d_utils.find_which_slit(output_hdul)
        if (slit is not None) or (slit != "NULL"):
            cwc_fname = output_hdul.next_input_file.replace(".fits", "_world_coordinates.fits")
            fdata = core_utils.get_filedata(cwc_fname, ext=ext)
            pwave = fdata[0,:,:]
            pdy = fdata[3,:,:]
            pskyx = fdata[1,:,:]
            pskyy = fdata[2,:,:]
            # get the origin of the subwindow and the grating from the extract_2d file header
            det = extract_2d_utils.find_DETECTOR(output_hdul)
            # the science extensions are always the same in the extract_2d file
            e2d_fname = output_hdul.next_input_file
            sci_data_list, px0_list, py0_list = [], [], []
            grat = core_utils.get_keywd_val(e2d_fname, "GRATING", ext=0)
            filt = core_utils.get_keywd_val(e2d_fname, "FILTER", ext=0)
            print (grat, filt)
            sci_exts = core_utils.get_sci_extensions(e2d_fname)
            for se in sci_exts:
                sci_data_list.append(core_utils.get_filedata(e2d_fname, ext=se))
                px0_se = core_utils.get_keywd_val(e2d_fname, "SLTSTRT1", ext=se)+core_utils.get_keywd_val(e2d_fname, "SUBSTRT1", ext=0)-1
                py0_se = core_utils.get_keywd_val(e2d_fname, "SLTSTRT2", ext=se)+core_utils.get_keywd_val(e2d_fname, "SUBSTRT2", ext=0)-1
                px0_list.append(px0_se)
                py0_list.append(py0_se)
            n_p = np.shape(pwave)
            npx = n_p[0]
            npy = n_p[1]
            px = np.arange(npx)+px0_list
            py = np.arange(npy)+py0_list
            print  ("Pipeline subwindow corner pixel ID: ", px0_list, py0_list)
            # read in ESA data

            #working_dir = os.path.dirname(output_hdul.output_file)

    elif extract_2d_utils.check_MOS_true(output_hdul):
        #fdata = core_utils.get_filedata(output_hdul.output_file, ext=None)
    else:
        pytest.skip("The fits file is not either FS nor MOS.")

# Unit tests

def test_s_ext2d_exists(output_hdul):
    assert extract_2d_utils.s_ext2d_exists(output_hdul)
