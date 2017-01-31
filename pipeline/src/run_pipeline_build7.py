#! /usr/bin/env python

from jwst.pipeline import SloperPipeline
from jwst.pipeline import Spec2Pipeline
from configobj import ConfigObj
import argparse
import numpy as np
import os


"""
This script runs the pipeline with the structure of build 7.

Keyword arguments:
    input_file  -- Name of the fits file to pass through the code

Output(s):
    Final product of pipeline

Example usage:
    In a terminal type the following (make sure you call the code from the right directory)
    to run the pipeline up to the calwebb_sloper step:
    > python run_pipeline_build7.py -fs -nrs2 jwtest1008001_01101_00001_NRS2_uncal_mod.fits

    or for the 1b and 2a levels type:
    > python run_pipeline_build7.py -rlev2A -fs -nrs2 jwtest1008001_01101_00001_NRS2_uncal_mod.fits

      Running options:
      -fs = Fixed Slit
      -ifu = Integral Field Unit
      -msa = Multi Shutter Array
      -nrs1 and -nrs2 = use if you want to override the configuration files/steps
      -rlev2A = will stop after CalwebSpec2 (next step from level 1b SloperPipeline)

 * NOTE - This script will run fine as long as you have the following directory structure:
          for the modified reference files --> /nirspec/pipeline/build6/test_data/Reffiles/
          for the configuration files --> /nirspec/pipeline/build7/config_files/

"""

# Header
__author__ = "M. A. Pena-Guerrero"
__version__ = "1.0"

# Get arguments to run script
parser = argparse.ArgumentParser(description='')
parser.add_argument("input_file",
                    action='store',
                    default=None,
                    help='Name of input fits file, i.e. blah.fits')
parser.add_argument("-nrs1",
                    dest="nrs1",
                    action='store_true',
                    default=False,
                    help='Shutter velocity.')
parser.add_argument("-nrs2",
                    dest="nrs2",
                    action='store_true',
                    default=False,
                    help='Shutter velocity.')
parser.add_argument("-fs",
                    dest="fs",
                    action='store_true',
                    default=False,
                    help='Fixed slit.')
parser.add_argument("-fs1600",
                    dest="fs1600",
                    action='store_true',
                    default=False,
                    help='Fixed slit, S1600 subarray')
parser.add_argument("-ifu",
                    dest="ifu",
                    action='store_true',
                    default=False,
                    help='Integral Field Unit.')
parser.add_argument("-msa",
                    dest="msa",
                    action='store_true',
                    default=False,
                    help='Multi Shutter Array.')
parser.add_argument("-rlev2A",
                    dest="rlev2A",
                    action='store_true',
                    default=False,
                    help='Stop the script after the level rlev2A product is obtained, i.e. SloperPipeline is done.')
args = parser.parse_args()

# Set the variables
input_file = args.input_file
nrs1 = args.nrs1
nrs2 = args.nrs2
fs = args.fs
fs1600 = args.fs1600
ifu = args.ifu
msa = args.msa
rlev2A = args.rlev2A

# Paths
pileline_path = os.path.abspath(os.curdir).split('pipeline')[0]
build6_path = os.path.join(pileline_path, 'pipeline/build6')
build7_path = os.path.join(pileline_path, 'pipeline/build7')
config_files_path = os.path.join(build7_path, 'config_files')
test_data_path = os.path.join(build6_path, 'test_data')
FS_data_path = os.path.join(test_data_path, 'FS_data')
calwebb_sloper = os.path.join(config_files_path, 'calwebb_sloper.cfg')
calwebb_spec2 = os.path.join(config_files_path, 'calwebb_spec2.cfg')

# Since there are a couple of steps that we want to override, we want to copy the config file into
# the current directory, modify it, run the pipeline, and then erase it. But only do this if one
# of the shutter velocities is set to true. The velocity will be used to determine which file to
# use in the override.
modify_config_file = False
if nrs1 or nrs2:
    modify_config_file = True

# For the Fixed Slit files we are using subarray=ALLSLITS for all the files I am testing
if fs:
    if nrs1:
        newdark = test_data_path+'/Reffiles/jwst_nirspec_dark_0028_dms.fits'
        newbias = test_data_path+'/Reffiles/jwst_nirspec_superbias_0032_dms.fits'
    if nrs2:
        newdark = test_data_path+'/Reffiles/jwst_nirspec_dark_0036_dms.fits'
        newbias = test_data_path+'/Reffiles/jwst_nirspec_superbias_0040_dms.fits'

if fs1600:
    if nrs1 or nrs2:
        newdark = test_data_path+'/Reffiles/jwst_nirspec_dark_0029_mod.fits'
        newbias = test_data_path+'/Reffiles/jwst_nirspec_superbias_0033_mod.fits'

# For IFU and MSA files these are the new darks and superbias, subarray=FULL for all the files I am testing
if ifu or msa:
    if nrs1:
        newdark = test_data_path+'/Reffiles/jwst_nirspec_dark_0024_mod.fits'
        newbias = test_data_path+'/Reffiles/jwst_nirspec_superbias_0028_mod.fits'
    if nrs2:
        newdark = test_data_path+'/Reffiles/jwst_nirspec_dark_0032_mod.fits'
        newbias = test_data_path+'/Reffiles/jwst_nirspec_superbias_0036_mod.fits'

if modify_config_file:
    # copy the config file into the current directory
    os.system('cp '+calwebb_sloper+' .')
    new_copy_calwebbsloper = os.path.abspath('calwebb_sloper.cfg')
    os.system('cp '+new_copy_calwebbsloper+' '+new_copy_calwebbsloper.replace('.cfg','_original.cfg'))
    # now modify the file
    calwebb_sloper = new_copy_calwebbsloper
    config = ConfigObj(calwebb_sloper)
    # modify paths of config files to use and add relevant lines
    config['save_calibrated_ramp'] = 'True'
    config['steps']['dq_init']['config_file'] = config_files_path+'/dq_init.cfg'
    config['steps']['saturation']['config_file'] = config_files_path+'/saturation.cfg'
    config['steps']['superbias']['override_superbias'] = newbias
    config['steps']['refpix']['config_file'] = config_files_path+'/refpix.cfg'
    config['steps']['rscd']['config_file'] = config_files_path+'/rscd.cfg'
    config['steps']['lastframe']['config_file'] = config_files_path+'/lastframe.cfg'
    config['steps']['linearity']['config_file'] = config_files_path+'/linearity.cfg'
    config['steps']['dark_current']['override_dark'] = newdark
    config['steps']['persistence']['config_file'] = config_files_path+'/persistence.cfg'
    config['steps']['jump']['config_file'] = config_files_path+'/jump.cfg'
    config['steps']['ramp_fit']['config_file'] = config_files_path+'/ramp_fit.cfg'
    # remove unused lines
    config['steps']['superbias'].pop('config_file', None)
    config['steps']['dark_current'].pop('config_file', None)
    config.write()

# Run the pipeline
print ('I am using the followging configuration file: ', calwebb_sloper)
result_level2A = SloperPipeline.call(input_file, config_file=calwebb_sloper)
print ('\n OK, I finished the level 2A, here is the result: ')
print ('Level 2A shape of resulting file:', repr(np.shape(result_level2A)) +'\n')

# Onto level 2B
if rlev2A:
    rlev2A = input_file.split('_uncal')[0]+'_uncal_rateins.fits'
    result_level2B = Spec2Pipeline.call(rlev2A, config_file=calwebb_spec2)
    print ('\n OK, I finished the level 2B, here is the result: ')
    print ('Level 2B shape of resulting file:', repr(np.shape(result_level2B)) +'\n')


print ('\n Script run_pipeline_build7 finished. ')
