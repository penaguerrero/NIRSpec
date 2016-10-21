#! /usr/bin/env python

from jwst_pipeline.pipeline import SloperPipeline
from jwst_pipeline.pipeline import Spec2Pipeline
from astropy.io import fits
from configobj import ConfigObj
from glob import glob
import argparse
import numpy as np
import string
import os


"""
This script runs the pipeline with the structure of build 6.

Keyword arguments:
    input_file  -- Name of the fits file to pass through the code

Output(s):
    final product of pipeline

Example usage:
    In a terminal type the following (make sure you call the code from the right directory):
    > python run_pipeline_build6.py -fs -nrs1 blah.fits

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
parser.add_argument("-ifu",
                    dest="ifu",
                    action='store_true',
                    default=False,
                    help='Integral Field Unit.')
args = parser.parse_args()

# Set the variables
input_file = args.input_file
nrs1 = args.nrs1
nrs2 = args.nrs2
fs = args.fs
ifu = args.ifu

# Paths
pileline_path = string.split(os.getcwd(), sep='pipeline')[0]
build6_path = os.path.join(pileline_path, 'pipeline/build6')
config_files_path = os.path.join(build6_path, 'config_files')
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

# For IFU files these are the new darks and superbias, subarray=FULL for all the files I am testing
if ifu:
    if nrs1:
        newdark = test_data_path+'/Reffiles/jwst_nirspec_dark_0024_mod.fits'
        newbias = test_data_path+'/Reffiles/jwst_nirspec_superbias_0028_mod.fits'
    if nrs2:
        newdark = test_data_path+'/Reffiles/jwst_nirspec_dark_0032_mod.fits'
        newbias = test_data_path+'/Reffiles/jwst_nirspec_superbias_0036_mod.fits'


if modify_config_file:
    # copy the config file into the current directory
    src_dir_path = pileline_path+'pipeline/src'
    os.system('cp '+calwebb_sloper+' '+src_dir_path)
    os.system('cp '+src_dir_path+'/calwebb_sloper.cfg'+' '+src_dir_path+'/calwebb_sloper_original.cfg')
    calwebb_sloper = os.path.join(src_dir_path, 'calwebb_sloper.cfg')
    # now modify the file
    config = ConfigObj(calwebb_sloper)
    # modify paths of config files to use
    config['steps']['dq_init']['config_file'] = config_files_path+'/dq_init.cfg'
    config['steps']['saturation']['config_file'] = config_files_path+'/saturation.cfg'
    config['steps']['refpix']['config_file'] = config_files_path+'/refpix.cfg'
    config['steps']['reset']['config_file'] = config_files_path+'/reset.cfg'
    config['steps']['lastframe']['config_file'] = config_files_path+'/lastframe.cfg'
    config['steps']['linearity']['config_file'] = config_files_path+'/linearity.cfg'
    config['steps']['jump']['config_file'] = config_files_path+'/jump.cfg'
    config['steps']['ramp_fit']['config_file'] = config_files_path+'/ramp_fit.cfg'
    # add relevant lines
    config['steps']['superbias']['override_superbias'] = newbias
    config['steps']['dark_current']['override_dark'] = newdark
    # remove unused lines
    config['steps']['superbias'].pop('config_file', None)
    config['steps']['dark_current'].pop('config_file', None)
    config.write()

# Run the pipeline
result_level2A = SloperPipeline.call(input_file, config_file=calwebb_sloper)
print ('\n OK, I finished the level 2A, here is the result: ')
print ('Level 2A shape of resulting file:', repr(np.shape(result_level2A)) +'\n')

# Onto level 2B
#rlev2A = string.split(input_file, sep='_uncal')[0]+'_uncal_rateints.fits'
#result_level2B = Spec2Pipeline.call(rlev2A, config_file=calwebb_spec2)
#print ('\n OK, I finished the level 2B, here is the result: ')
#print ('Level 2B shape of resulting file:', repr(np.shape(result_level2B)) +'\n')


print ('\n Script run_pipeline_build6 finished. ')
