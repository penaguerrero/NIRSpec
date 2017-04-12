"""
py.test module for unit testing the bkg_subtract step.
"""

import pytest
import os
from jwst.background.background_step import BackgroundStep

from .. import core_utils
from . import bkg_subtract_utils


# Set up the fixtures needed for all of the tests, i.e. open up all of the FITS files

# Default names of pipeline output files
prev_file_name = "_assign_wcs"
out_file_name = "_subtract_images"

# fixture to read the config file
@pytest.fixture(scope="module")
def input_hdul(request, config):
    initiate_calwebb_spc2 = "calwebb_spec2_input_file"
    if config.has_option(initiate_calwebb_spc2, "input_file"):
        working_directory = config.get(initiate_calwebb_spc2, "working_directory")
        initial_input_file = config.get(initiate_calwebb_spc2, "input_file")
        initial_input_file_fullpath = os.path.join(working_directory, initial_input_file)
        next_input_file = initial_input_file_fullpath.replace(".fits", prev_file_name+".fits")
        if os.path.isfile(next_input_file):
            hdul = core_utils.read_hdrfits(next_input_file, info=True, show_hdr=True)
            return hdul, next_input_file
        else:
            pytest.skip("skiping bkg_subtract because there are no ASN files")
    else:
        pytest.skip("bkg_subtract needs an input_file from assign_wcs step")


# fixture to read the output file header
@pytest.fixture(scope="module")
def output_hdul(request, config):
    step = "bkg_subtract"
    initiate_calwebb_spc2 = "calwebb_spec2_input_file"
    working_directory = config.get(initiate_calwebb_spc2, "working_directory")
    initial_input_file = config.get(initiate_calwebb_spc2, "input_file")
    initial_input_file_fullpath = os.path.join(working_directory, initial_input_file)
    next_input_file = initial_input_file_fullpath.replace(".fits", prev_file_name+".fits")
    output_file = next_input_file.replace(".fits", out_file_name+".fits")
    stp = BackgroundStep()
    run_calwebb_spec2 = config.get("run_calwebb_spec2_in_full", "run_calwebb_spec2")
    # if run_calwebb_spec2 is True calwebb_spec2 will be called, else individual steps will be ran
    if not run_calwebb_spec2:
        print ("Will run "+step+" step ...")
        #result = stp.call(next_input_file)
        #result.save(output_file)
    skip_step = True
    if config.has_option("steps", step):
        # since for now this step is not being called, the test will be skipped until
        # I can implement a routine to test if there are ASN files with the input file.
        if not skip_step:
            hdul = core_utils.read_hdrfits(output_file, info=True, show_hdr=True)
            return hdul
        else:
            pytest.skip("skiping bkg_subtract because there are no ASN files")
    else:
        pytest.skip("needs bkg_subtract output_file")



# Unit tests

def test_s_bkdsub_exists(output_hdul):
    assert bkg_subtract_utils.s_bkdsub_exists(output_hdul)
