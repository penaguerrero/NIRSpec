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
from ..auxiliary_code import compare_wcs_mos
from ..auxiliary_code import compare_wcs_ifu


# Set up the fixtures needed for all of the tests, i.e. open up all of the FITS files

# Default names of pipeline input and output files
@pytest.fixture(scope="module")
def set_inandout_filenames(request, config):
    step = "extract_2d"
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
    stp = Extract2dStep()
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



### THESE FUNCTIONS ARE TO VALIDATE BOTH THE WCS AND THE 2D_EXTRACT STEPS

# fixture to validate the WCS and extract 2d steps
@pytest.fixture(scope="module")
def validate_wcs_extract2d(output_hdul):
    # get the input information for the wcs routine
    infile_name = output_hdul.initial_input_file.replace(".fits", "_assign_wcs_extract_2d.fits")
    esa_files_path = output_hdul.esaroot

    # define the threshold difference between the pipeline output and the ESA files for the pytest to pass or fail
    threshold_diff = 1.0e-14

    if extract_2d_utils.check_FS_true(output_hdul):
        # Find what slit the data corresponds to
        ext, slit = extract_2d_utils.find_which_slit(output_hdul)
        if (slit is not None) or (slit != "NULL"):
            median_diff = compare_wcs_fs.compare_wcs(infile_name, esa_files_path=esa_files_path,
                                                     auxiliary_code_path=None, plot_names=None,
                                                     show_figs=False, save_figs=False,
                                                  threshold_diff=threshold_diff)

    elif extract_2d_utils.check_MOS_true(output_hdul):
        msa_conf_root = output_hdul.msa_conf_root
        median_diff = compare_wcs_mos.compare_wcs(infile_name, msa_conf_root=msa_conf_root,
                                                  esa_files_path=esa_files_path, auxiliary_code_path=None,
                                                  plot_names=None, show_figs=False, save_figs=False,
                                                  threshold_diff=threshold_diff)

    elif extract_2d_utils.check_IFU_true(output_hdul):
        median_diff = compare_wcs_ifu.compare_wcs(infile_name, esa_files_path=esa_files_path, auxiliary_code_path=None,
                                                  plot_names=None, show_figs=False, save_figs=False,
                                                  threshold_diff=threshold_diff)
    else:
        pytest.skip("Skipping pytest: The fits file is not FS, MOS, or IFU. Pytest does not yet include the routine to verify this kind of file.")
    return median_diff



### Unit tests

def test_s_ext2d_exists(output_hdul):
    assert extract_2d_utils.s_ext2d_exists(output_hdul)


def test_validate_wcs_extract2d(output_hdul):
    assert validate_wcs_extract2d(output_hdul)
