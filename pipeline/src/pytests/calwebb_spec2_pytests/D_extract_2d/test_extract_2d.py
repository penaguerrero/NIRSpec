from __future__ import print_function, division

"""
py.test module for unit testing the extract_2d step.
"""

import pytest
import os
from jwst.extract_2d.extract_2d_step import Extract2dStep

from .. import core_utils
from . import extract_2d_utils
from ..auxiliary_code import compare_wcs_fs


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
    # get information from the configuration file
    esaroot = config.get("esa_intermediary_products", "esaroot")
    initiate_calwebb_spc2 = "calwebb_spec2_input_file"
    working_directory = config.get(initiate_calwebb_spc2, "working_directory")
    initial_input_file = config.get(initiate_calwebb_spc2, "input_file")
    #data_directory = config.get(initiate_calwebb_spc2, "data_directory")
    #initial_input_file_fullpath = os.path.join(data_directory, initial_input_file)
    next_input_file = os.path.join(working_directory, initial_input_file.replace(".fits", prev_file_name+".fits"))
    if os.path.isfile(next_input_file):
        output_file = next_input_file.replace(".fits", out_file_name+".fits")
    else:
        next_input_file = os.path.join(working_directory, initial_input_file.replace(".fits", no_ASN_files_in+".fits"))
        output_file = next_input_file.replace(".fits", no_ASN_files_out+".fits")
    stp = Extract2dStep()
    run_calwebb_spec2 = config.get("run_calwebb_spec2_in_full", "run_calwebb_spec2")
    # if run_calwebb_spec2 is True calwebb_spec2 will be called, else individual steps will be ran
    if not run_calwebb_spec2:
        print ("Will run "+step+" step ...")
        #result = stp.call(next_input_file)
        #result.save(output_file)
    #if config.has_option("steps", step):
    if config.get("steps", step):
        hdul = core_utils.read_hdrfits(output_file, info=True, show_hdr=True)
        return initial_input_file, hdul, next_input_file, esaroot
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
            compare_wcs_fs.compare_wcs(output_hdul.hdul, cwc_fname, ext, output_hdul.esaroot)

    elif extract_2d_utils.check_MOS_true(output_hdul):

        #fdata = core_utils.get_filedata(output_hdul.output_file, ext=None)
    else:
        pytest.skip("The fits file is not either FS nor MOS.")

# Unit tests

def test_s_ext2d_exists(output_hdul):
    assert extract_2d_utils.s_ext2d_exists(output_hdul)
