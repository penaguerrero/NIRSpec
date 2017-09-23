from __future__ import print_function, division
"""
py.test module for unit testing the msa_flagging step.
"""

import pytest
import os
from jwst.imprint.imprint_step import ImprintStep

from .. import core_utils
from . import msa_flagging_utils


# Set up the fixtures needed for all of the tests, i.e. open up all of the FITS files

# Default names of pipeline output files
prev_file_name = "_assign_wcs_subtract_images_imprint"
out_file_name = "_msa_flag"

# This is the default name the pipeline will give an output file that did not have association (ASN) files
# for the background and imprint subtraction steps
no_ASN_files_in = "_assign_wcs"
no_ASN_files_out = "_msa_flag"

# fixture to read the config file
@pytest.fixture(scope="module")
def input_hdul(request, config):
    initiate_calwebb_spc2 = "calwebb_spec2_input_file"
    if config.has_option(initiate_calwebb_spc2, "input_file"):
        working_directory = config.get(initiate_calwebb_spc2, "working_directory")
        initial_input_file = config.get(initiate_calwebb_spc2, "input_file")
        initial_input_file_fullpath = os.path.join(working_directory, initial_input_file)
        step_input_file = initial_input_file_fullpath.replace(".fits", prev_file_name+".fits")
        if os.path.isfile(step_input_file):
            output_file = step_input_file.replace(".fits", out_file_name+".fits")
        else:
            step_input_file = os.path.join(working_directory, initial_input_file.replace(".fits", no_ASN_files_in+".fits"))
            output_file = step_input_file.replace(".fits", no_ASN_files_out+".fits")
        hdul = core_utils.read_hdrfits(step_input_file, info=True, show_hdr=True)
        return hdul, step_input_file
    else:
        pytest.skip("skipping extract_2d test, step needs an input_file")


# fixture to read the output file header
@pytest.fixture(scope="module")
def output_hdul(request, config):
    step = "imprint_subtract"
    initiate_calwebb_spc2 = "calwebb_spec2_input_file"
    working_directory = config.get(initiate_calwebb_spc2, "working_directory")
    initial_input_file = config.get(initiate_calwebb_spc2, "input_file")
    initial_input_file_fullpath = os.path.join(working_directory, initial_input_file)
    step_input_file = initial_input_file_fullpath.replace(".fits", prev_file_name+".fits")
    output_file = step_input_file.replace(".fits", out_file_name+".fits")
    stp = ImprintStep()
    run_calwebb_spec2 = config.get("run_calwebb_spec2_in_full", "run_calwebb_spec2")
    # if run_calwebb_spec2 is True calwebb_spec2 will be called, else individual steps will be ran
    if not run_calwebb_spec2:
        if config.get("steps", step):
            print ("Will run "+step+" step ...")
            if os.path.isfile(step_input_file):
                #result = stp.call(step_input_file)
                #result.save(output_file)
                hdul = core_utils.read_hdrfits(output_file, info=True, show_hdr=True)
                return hdul, output_file
        else:
            pytest.skip("Skiping "+step+" because the input file does not exist.")
    else:
        pytest.skip("Skiping "+step+". Step set to False in configuration file.")



# Unit tests

def test_msa_failed_open_exists(output_hdul):
    assert msa_flagging_utils.msa_failed_open_exists(output_hdul)
