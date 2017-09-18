from __future__ import print_function, division
"""
py.test module for unit testing the assign_wcs step.
"""

import pytest
import os
from jwst.pipeline import Spec2Pipeline
from jwst.assign_wcs.assign_wcs_step import AssignWcsStep

from .. import core_utils
from . import assign_wcs_utils


# Set up the fixtures needed for all of the tests, i.e. open up all of the FITS files


# fixture to read the config file
@pytest.fixture(scope="module")
def input_hdul(request, config):
    initiate_calwebb_spc2 = "calwebb_spec2_input_file"
    if config.has_option(initiate_calwebb_spc2, "input_file"):
        input_file = config.get(initiate_calwebb_spc2, "input_file")
        hdul = core_utils.read_hdrfits(input_file, info=True, show_hdr=True)
        return hdul, input_file
    else:
        pytest.skip("needs an input_file to start running calwebb_spec2")


# fixture to read the output file header
@pytest.fixture(scope="module")
def output_hdul(request, config):
    step = "assign_wcs"
    initiate_calwebb_spc2 = "calwebb_spec2_input_file"
    working_directory = config.get(initiate_calwebb_spc2, "working_directory")
    input_file = config.get(initiate_calwebb_spc2, "input_file")
    data_directory = config.get(initiate_calwebb_spc2, "data_directory")
    input_file = os.path.join(data_directory, input_file)
    output_file = input_file.replace(".fits", "_assign_wcs.fits")
    output_file = os.path.join(working_directory, output_file)
    stp = AssignWcsStep()
    run_calwebb_spec2 = config.get("run_calwebb_spec2_in_full", "run_calwebb_spec2")
    # if run_calwebb_spec2 is True calwebb_spec2 will be called, else individual steps will be ran
    if run_calwebb_spec2:
        print ("Will run calwebb_spec2... ")
        calwebb_spec2_cfg = config.get("run_calwebb_spec2_in_full", "calwebb_spec2_cfg")
        result_level2B = Spec2Pipeline.call(input_file, config_file=calwebb_spec2_cfg)
        #result_level2B.save(newfilename)
    else:
        # run individual step of the pipeline
        print ("Will run assign_wcs step...")
        #result = stp.call(input_file)
        #result.save(output_file)
    if config.has_option("steps", step):
        hdul = core_utils.read_hdrfits(output_file, info=True, show_hdr=True)
        return hdul, output_file
    else:
        pytest.skip("needs assign_wcs output_file")


'''
@pytest.fixture(scope="module")
def input_hdul(request, config):
    step = "B_assign_wcs"
    if  config.has_option(step, "input_file"):
        hdul, fdata = core_utils.read_fits(config.get(step, "input_file"), info=True, show_hdr=True)
        return hdul
    else:
        pytest.skip("needs assign_wcs input_file")

@pytest.fixture(scope="module")
def output_hdul(request, config):
    step = "B_assign_wcs"
    if  config.has_option(step, "output_file"):
        hdul, fdata = core_utils.read_fits(config.get(step, "output_file"), info=True, show_hdr=True)
        return hdul
    else:
        pytest.skip("needs assign_wcs output_file")
'''

# Unit tests

def test_wavstart_exists(output_hdul):
    assert assign_wcs_utils.wavstart_exists(output_hdul)

def test_sporder_exists(output_hdul):
    assert assign_wcs_utils.sporder_exists(output_hdul)

def test_s_wcs_exists(output_hdul):
    assert assign_wcs_utils.s_wcs_exists(output_hdul)

