"""
py.test module for unit testing the assign_wcs step.
"""

import pytest
import os
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
    working_directory = config.get("calwebb_spec2_input_file", "working_directory")
    input_file = config.get("calwebb_spec2_input_file", "input_file")
    input_file = os.path.join(working_directory, input_file)
    output_file = input_file.replace(".fits", "_assign_wcs.fits")
    output_file = os.path.join(working_directory, output_file)
    stp = AssignWcsStep()
    #result = stp.call(input_file)
    #result.save(output_file)
    if config.has_option("steps", step):
        hdul = core_utils.read_hdrfits(output_file, info=True, show_hdr=True)
        return hdul
    else:
        pytest.skip("needs assign_wcs output_file")


'''
@pytest.fixture(scope="module")
def input_hdul(request, config):
    step = "A_assign_wcs"
    if  config.has_option(step, "input_file"):
        hdul = core_utils.read_hdrfits(config.get(step, "input_file"), info=True, show_hdr=True)
        return hdul
    else:
        pytest.skip("needs assign_wcs input_file")

@pytest.fixture(scope="module")
def output_hdul(request, config):
    step = "A_assign_wcs"
    if  config.has_option(step, "output_file"):
        hdul = core_utils.read_hdrfits(config.get(step, "output_file"), info=True, show_hdr=True)
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

