from __future__ import print_function, division
"""
py.test module for unit testing the bkg_subtract step.
"""

import pytest
import os
from jwst.pipeline import Spec2Pipeline
from jwst.background.background_step import BackgroundStep

from .. import core_utils
from . import bkg_subtract_utils


# Set up the fixtures needed for all of the tests, i.e. open up all of the FITS files

# Default names of pipeline output files
#prev_file_name = ""
out_file_name = "_subtract_images"

# fixture to read the config file
@pytest.fixture(scope="module")
def input_hdul(request, config):
    initiate_calwebb_spc2 = "calwebb_spec2_input_file"
    if config.has_option(initiate_calwebb_spc2, "input_file"):
        working_directory = config.get(initiate_calwebb_spc2, "working_directory")
        initial_input_file = config.get(initiate_calwebb_spc2, "input_file")
        #data_directory = config.get(initiate_calwebb_spc2, "data_directory")
        #initial_input_file_fullpath = os.path.join(data_directory, initial_input_file)
        step_input_file = os.path.join(working_directory, initial_input_file)
        if os.path.isfile(step_input_file):
            hdul = core_utils.read_hdrfits(step_input_file, info=True, show_hdr=True)
            return hdul, step_input_file
        else:
            pytest.skip("skiping bkg_subtract because there are no ASN files")
    else:
        pytest.skip("skiping the bkg_subtract test")


# fixture to read the output file header
@pytest.fixture(scope="module")
def output_hdul(request, config):
    step = "bkg_subtract"
    initiate_calwebb_spc2 = "calwebb_spec2_input_file"
    working_directory = config.get(initiate_calwebb_spc2, "working_directory")
    initial_input_file = config.get(initiate_calwebb_spc2, "input_file")
    step_input_file = os.path.join(working_directory, initial_input_file)
    output_file = step_input_file.replace(".fits", out_file_name+".fits")
    stp = BackgroundStep()
    run_calwebb_spec2 = config.get("run_calwebb_spec2_in_full", "run_calwebb_spec2")
    # if run_calwebb_spec2 is True calwebb_spec2 will be called, else individual steps will be ran
    if run_calwebb_spec2:
        print ("Will run calwebb_spec2... ")
        calwebb_spec2_cfg = config.get("run_calwebb_spec2_in_full", "calwebb_spec2_cfg")
        final_output_name = step_input_file.replace(".fits", "_calwebb_spec2.fits")
        #result_level2B = Spec2Pipeline.call(step_input_file, config_file=calwebb_spec2_cfg)
        #result_level2B.save(final_output_name)
    else:
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

def test_s_bkdsub_exists(output_hdul):
    assert bkg_subtract_utils.s_bkdsub_exists(output_hdul)
