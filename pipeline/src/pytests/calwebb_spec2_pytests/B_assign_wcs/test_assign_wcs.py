from __future__ import print_function, division
"""
py.test module for unit testing the assign_wcs step.
"""

import pytest
import os
from jwst.assign_wcs.assign_wcs_step import AssignWcsStep

from .. import core_utils
from . import assign_wcs_utils


# Set up the fixtures needed for all of the tests, i.e. open up all of the FITS files

# Default names of pipeline output files
prev_file_name = "_subtract_images"
out_file_name = "_assign_wcs"

# fixture to read the config file
@pytest.fixture(scope="module")
def input_hdul(request, config):
    initiate_calwebb_spc2 = "calwebb_spec2_input_file"
    if config.has_option(initiate_calwebb_spc2, "input_file"):
        input_file = config.get(initiate_calwebb_spc2, "input_file")
        data_directory = config.get(initiate_calwebb_spc2, "data_directory")
        step_input_file = os.path.join(data_directory, input_file).replace(".fits", prev_file_name+".fits")
        hdul = core_utils.read_hdrfits(input_file, info=True, show_hdr=True)
        return hdul, step_input_file
    else:
        pytest.skip("Skipping assign_wcs test")


# fixture to read the output file header
@pytest.fixture(scope="module")
def output_hdul(request, config):
    step = "assign_wcs"
    initiate_calwebb_spc2 = "calwebb_spec2_input_file"
    working_directory = config.get(initiate_calwebb_spc2, "working_directory")
    input_file = config.get(initiate_calwebb_spc2, "input_file")
    data_directory = config.get(initiate_calwebb_spc2, "data_directory")
    step_input_file = os.path.join(data_directory, input_file).replace(".fits", prev_file_name+".fits")
    output_file = input_file.replace(".fits", "_assign_wcs.fits")
    output_file = os.path.join(working_directory, output_file)
    stp = AssignWcsStep()
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

def test_wavstart_exists(output_hdul):
    assert assign_wcs_utils.wavstart_exists(output_hdul)

def test_sporder_exists(output_hdul):
    assert assign_wcs_utils.sporder_exists(output_hdul)

def test_s_wcs_exists(output_hdul):
    assert assign_wcs_utils.s_wcs_exists(output_hdul)

