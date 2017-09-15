from __future__ import print_function, division
import json
import os
import matplotlib.pyplot as plt
from pandeia.engine.calc_utils import build_default_calc
from pandeia.engine.io_utils import read_json, write_json

# utility functions for building calculations and sources.
from pandeia.engine.calc_utils import build_default_calc, build_default_source

# this is the core function that takes the calculation input, peforms the ETC calculation, and returns the results.
from pandeia.engine.perform_calculation import perform_calculation


"""
This script creates a json file with a NIRSpec TA calculation for a given filter and 
normalization magnitude.
It assumes a flat source with the normalization on abmag and on a bandpass.
"""


#######################################################################################

"""
    
Allowed values and strings for the setup variables

filter = "f110w", "f140x", "clear"
readmode = "nrsrapid", "nrs", "nrsrapidd1", nrsrapidd2"
subarray = "full", "sub32"
norm_mag = float in abmag

* NOTE: the json file to change the aperture size lives at:
/nirspec/python_ETC_engine/pandeia/v1.1.1/pandeia_data-1.1.1/jwst/nirspec/config.json

"""

# Script setup variables
json_file = "nirspec_ta_test.json"
filter = "f110w"
readmode = "nrsrapid"
subarray = "full"
norm_mag = 22.0 # in abmag


#######################################################################################


# other configurable params
aperture_size = 0.53   # in arcsec, 1 pixel = 0.106 arcsec
background = 'medium'  # choices of 'low', 'medium', 'high', 'none'
ngroups = 3
nint = 1
nexp = 1


# Get the paths right, this routine will work as long as:
#    1- The jason files are in directory batch_TA, which is expected to be in the
#       following path structure: /path_to_this_dir/python_ETC_engine/general_examples
#    2- This script is in a src directory at the same level as /python_ETC_engine/
working_dir_path_for_this_script = os.getcwd()
t = working_dir_path_for_this_script.split("/")
json_files_path = working_dir_path_for_this_script.replace(t[-1], "general_examples")

# Create the json file
# make a default NIRSpec TA calculation
c = build_default_calc(telescope="jwst", instrument="nirspec", mode="target_acq")

# use NIRCam SW imaging normalization and use a flat source
# set the source
c['scene'][0]['spectrum']['sed']['sed_type'] = 'flat'
c['scene'][0]['spectrum']['sed']['unit'] = 'fnu'
# set the configuration, background, and strategy params
c['configuration']['instrument']['filter'] = filter
c['configuration']['detector']['readmode'] = readmode
c['configuration']['detector']['subarray'] = subarray
c['configuration']['detector']['ngroup'] = ngroups
c['configuration']['detector']['nint'] = nint
c['configuration']['detector']['nexp'] = nexp
c['background'] = background
c['strategy']['aperture_size'] = aperture_size
# set normalization params
c['scene'][0]['spectrum']['normalization']['type'] = 'jwst'
c['scene'][0]['spectrum']['normalization']['bandpass'] = "nircam,sw_imaging,f150w2"
c['scene'][0]['spectrum']['normalization']['norm_flux'] = norm_mag
c['scene'][0]['spectrum']['normalization']['norm_fluxunit'] = 'abmag'

# write the json file
write_json(c, json_file)
print(json_file, " written")

# Load scene file
jf = os.path.join(json_files_path, json_file)
with open(jf) as f:
    imgr_data = json.load(f)   # this is a python dictionary

# Perform calculation
results = perform_calculation(imgr_data)   # results is a dictionary

# Print the results on-screen
print ("filter =", results["scalar"]["filter"])
print ("subarray: ", results["input"]["configuration"]["detector"]["subarray"])
print ("readout mode: ", results["input"]["configuration"]["detector"]["readmode"])
print ("warnings: ", results["warnings"])
print ("number of groups = ", results["information"]["exposure_specification"]["ngroup"])
print ("number of integrations = ", results["information"]["exposure_specification"]["nint"])
print ("number of exposures = ", results["information"]["exposure_specification"]["nexp"])
print ("aperture size [arcsec]: ", results["scalar"]["aperture_size"])
print ("total exposure time [s]: ", results["scalar"]["total_exposure_time"])
print ("extracted flux [-e/s] =", results["scalar"]["extracted_flux"])
print ("S/N =", results["scalar"]["sn"])

# to investigate where other parameter might live, print the results dictionary.

"""
# Plot the results
f, ax = plt.subplots()
ax.imshow(results["2d"]["saturation"])
f.savefig("some_name.jpg")
plt.show()
"""

print("\n Script   run_single_calc_create_json.py   finished.")
