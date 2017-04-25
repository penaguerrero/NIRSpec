from __future__ import print_function, division
"""
py.test module for unit testing the imprint_subtract step.
"""

import pytest
import os
from jwst.imprint.imprint_step import ImprintStep

from .. import core_utils
from . import imprint_subtract_utils


# Set up the fixtures needed for all of the tests, i.e. open up all of the FITS files

# Default names of pipeline output files
prev_file_name = "_assign_wcs_subtract_images"
out_file_name = "_imprint"

# fixture to read the config file
@pytest.fixture(scope="module")
def input_hdul(request, config):
    initiate_calwebb_spc2 = "calwebb_spec2_input_file"
    if config.has_option(initiate_calwebb_spc2, "input_file"):
        working_directory = config.get(initiate_calwebb_spc2, "working_directory")
        initial_input_file = config.get(initiate_calwebb_spc2, "input_file")
        #data_directory = config.get(initiate_calwebb_spc2, "data_directory")
        #initial_input_file_fullpath = os.path.join(data_directory, initial_input_file)
        next_input_file = os.path.join(working_directory, initial_input_file.replace(".fits", prev_file_name+".fits"))
        ###initial_input_file_fullpath = os.path.join(working_directory, initial_input_file)
        ###next_input_file = initial_input_file_fullpath.replace(".fits", prev_file_name+".fits")
        if os.path.isfile(next_input_file):
            hdul = core_utils.read_hdrfits(next_input_file, info=True, show_hdr=True)
            return hdul, next_input_file
        else:
            pytest.skip("skiping imprint_subtract because there are no ASN files")

    else:
        pytest.skip("imprint_subtract needs an input_file from bkg_subtract")


# fixture to read the output file header
@pytest.fixture(scope="module")
def output_hdul(request, config):
    step = "imprint_subtract"
    initiate_calwebb_spc2 = "calwebb_spec2_input_file"
    working_directory = config.get(initiate_calwebb_spc2, "working_directory")
    initial_input_file = config.get(initiate_calwebb_spc2, "input_file")
    initial_input_file_fullpath = os.path.join(working_directory, initial_input_file)
    next_input_file = initial_input_file_fullpath.replace(".fits", prev_file_name+".fits")
    output_file = next_input_file.replace(".fits", out_file_name+".fits")
    stp = ImprintStep()
    run_calwebb_spec2 = config.get("run_calwebb_spec2_in_full", "run_calwebb_spec2")
    # if run_calwebb_spec2 is True calwebb_spec2 will be called, else individual steps will be ran
    if not run_calwebb_spec2:
        print ("Will run "+step+" step ...")
        #result = stp.call(next_input_file)
        #result.save(output_file)
    if config.has_option("steps", step):
        if os.path.isfile(next_input_file):
            hdul = core_utils.read_hdrfits(output_file, info=True, show_hdr=True)
            return hdul
        else:
            pytest.skip("skiping imprint_subtract because there are no ASN files")
    else:
        pytest.skip("needs imprint_subtract output_file")



# Unit tests

def test_s_imprint_exists(output_hdul):
    assert imprint_subtract_utils.s_imprint_exists(output_hdul)
