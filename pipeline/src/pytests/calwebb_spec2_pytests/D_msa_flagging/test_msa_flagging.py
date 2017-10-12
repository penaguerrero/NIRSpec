from __future__ import print_function, division
"""
py.test module for unit testing the msa_flagging step.
"""

import pytest
import os
from jwst.msaflagopen.msaflagopen_step import MSAFlagOpenStep

from .. import core_utils
from . import msa_flagging_utils


# Set up the fixtures needed for all of the tests, i.e. open up all of the FITS files

# Default names of pipeline input and output files
@pytest.fixture(scope="module")
def set_inandout_filenames(request, config):
    step = "imprint_subtract"
    step_dict = dict(config.items("steps"))
    initial_input_file = config.get("calwebb_spec2_input_file", "input_file")
    #step_input_filename, step_output_filename = core_utils.get_step_inandout_filename(step, initial_input_file, step_dict)
    suffix_and_filenames = core_utils.get_step_inandout_filename(step, initial_input_file, step_dict)
    in_file_suffix, out_file_suffix, step_input_filename, step_output_filename = suffix_and_filenames
    return initial_input_file, step, step_input_filename, step_output_filename


# fixture to read the output file header
@pytest.fixture(scope="module")
def output_hdul(set_inandout_filenames, config):
    initiate_calwebb_spc2 = "calwebb_spec2_input_file"
    working_directory = config.get(initiate_calwebb_spc2, "working_directory")
    initial_input_file = set_inandout_filenames[0]
    step = set_inandout_filenames[1]
    output_file = set_inandout_filenames[2]
    step_input_file = os.path.join(working_directory, initial_input_file)
    stp = MSAFlagOpenStep()
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
