from pandeia.engine.perform_calculation import perform_calculation
import json
import os
import matplotlib.pyplot as plt


#######################################################################################

# Script setup variables
json_file = "nirspec_ta_1431_2.json"

#######################################################################################

# Get the paths right, this routine will work as long as:
#    1- The jason files are in directory batch_TA, which is expected to be in the
#       following path structure: /path_to_this_dir/python_ETC_engine/general_examples
#    2- This script is in a src directory at the same level as /python_ETC_engine/
working_dir_path_for_this_script = os.getcwd()
t = working_dir_path_for_this_script.split("/")
json_files_path = working_dir_path_for_this_script.replace(t[-1], "general_examples")

"""
The norm_configuration.json contains the format for normalizing to different filters.
"""

# Load scene file
jf = os.path.join(json_files_path, json_file)
with open(jf) as f:
    imgr_data = json.load(f)   # this is a python dictionary

# Change the flux value of the source
#imgr_data['scene'][0]['spectrum']['normalization']['norm_flux'] = 100

# Perform calculation
results = perform_calculation(imgr_data)   # results is a dictionary

# Plot the results
f, ax = plt.subplots()
ax.imshow(results["2d"]["saturation"])
f.savefig("nirspec_1431_2d_saturation.jpg")
plt.show()

print("Script run_nirspec_single_calculation.py  finished.")