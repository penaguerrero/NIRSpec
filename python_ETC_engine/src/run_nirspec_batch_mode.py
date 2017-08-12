from pandeia.engine.perform_calculation import perform_calculation
import json
import os
import time
import numpy as np
import matplotlib.pyplot as plt


#######################################################################################

# Script setup variables
json_file = "nirspec_ta.json"

#######################################################################################

# Get the paths right, this routine will work as long as:
#    1- The jason files are in directory batch_TA, which is expected to be in the
#       following path structure: /path_to_this_dir/python_ETC_engine/batch_TA
#    2- This script is in a src directory at the same level as /python_ETC_engine/
working_dir_path_for_this_script = os.getcwd()
t = working_dir_path_for_this_script.split("/")
json_files_path = working_dir_path_for_this_script.replace(t[-1], "batch_TA")

"""
The norm_configuration.json contains the format for normalizing to different filters.
"""

# Start the timer to compute the whole running time
start_time = time.time()

# Load scene file
jf = os.path.join(json_files_path, json_file)
with open(jf) as f:
    data = json.load(f)   # this is a python dictionary

# Define some variables for saturation exploration
tot_none, tot_soft, tot_hard = [], [], []
flux_range = np.linspace(0.1, 1e8, 10)

# Loop over some flux range to discern saturation levels
for flux in flux_range:
    # Set the flux value
    data['scene'][0]['spectrum']['normalization']['norm_flux'] = flux
    
    # Perform a calculation
    results = perform_calculation(data)
    
    # Store the saturation so it is easily accessible
    saturation = results['2d']['saturation']
    tot_none.append(len(saturation[saturation==0]))
    tot_soft.append(len(saturation[saturation==1]))
    tot_hard.append(len(saturation[saturation==2]))

# Plot the values
f, ax = plt.subplots()
ax.plot(flux_range, tot_none, label="No saturation")
ax.plot(flux_range, tot_soft, label="Soft saturation")
ax.plot(flux_range, tot_hard, label="Hard saturation")
ax.set_ylabel("Pixels [counts]")
ax.set_xlabel("Flux [mJy]")
ax.legend(loc=0)
f.savefig("saturation_comparison.jpg")
#plt.show()

time2run = 'ETC engine finished! Took  %s  seconds tu run calculations.' % ( (time.time() - start_time) )
print(time2run)

print("Script run_nirspec_batch_mode.py  finished.")

