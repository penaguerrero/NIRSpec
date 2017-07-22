from pandeia.engine.perform_calculation import perform_calculation
import json
import matplotlib.pyplot as plt


#######################################################################################

# Script setup variables
json_files_path = "../general_examples"

#######################################################################################


# Load scene file
with open(json_files_path) as f:
    imgr_data = json.load(f)   # this is a python dictionary

# Change the flux value of the source
imgr_data['scene'][0]['spectrum']['normalization']['norm_flux'] = 100

# Perform calculation
results = perform_calculation(imgr_data)   # results is a dictionary

f, ax = plt.subplots()
ax.imshow(results["2d"]["saturation"])