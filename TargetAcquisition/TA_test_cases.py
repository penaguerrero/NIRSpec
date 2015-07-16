from __future__ import print_function, division
from astropy.io import fits
from glob import glob
import numpy as np
import matplotlib.pyplot as plt
import os

# Tommy's code
import tautils as tu
import jwst_targloc as jtl

# Header
__author__ = "Maria A. Pena-Guerrero"
__version__ = "1.0"


"""

DESCRIPTION:
    This script runs a series of test cases for the Target Locate code jwst_targloc.py 
    written by T. Le Blanc, and modified by M. A. Pena-Guerrero. Test cases are:
        - Performing centroiding for widths of 3, 5, and 7 pixels
        - Fractional backgrounds from 0.0 to 1.0 in increments of 0.1, and a fixed 
        or None background=0 case
        - Moments are included in the code though not currently being used
        
    The tests are performed on ESA data in the following directories on central sotrage:
        - /grp/jwst/wit4/nirspec/PFforMaria/Scene_1_AB23
            (which contains 200 closer to "real" stars of magnitude 23: cosmic rays, hot 
            pixels, noise, and some shutters closed; and an ideal case)
        - /grp/jwst/wit4/nirspec/PFforMaria/Scene_2_AB1823
            (which contains 200 stars of different magnitudes and also has the "real" data
            as well as the ideal case)
        * There are 2 sub-folders in each scene: rapid and slow. This is the shutter speed. Both
        sub-cases will be tested.
        ** The simulated data is described in detail in the NIRSpec Technical Note NTN-2015-013, which
        is in /grp/jwst/wit4/nirspec/PFforMaria/Documentation.
        
        
"""

# Paths to Scenes 1 and 2 local directories: /Users/pena/Documents/AptanaStudio3/NIRSpec/TargetAcquisition/
path_scene1_slow = "PFforMaria/Scene_1_AB23/NIRSpec_TA_Sim_AB23 first NRS/postage"
path_scene1_slow_nonoise = os.path.abspath("PFforMaria/Scene_1_AB23/NIRSpec_TA_Sim_AB23 first NRS no_noise/postage")
path_scene1_rapid = os.path.abspath("PFforMaria/Scene_1_AB23/NIRSpec_TA_Sim_AB23 first NRSRAPID/postage")
path_scene1_rapid_nonoise = os.path.abspath("PFforMaria/Scene_1_AB23/NIRSpec_TA_Sim_AB23 first NRS no_noise/postage")
path_scene2_slow = os.path.abspath("PFforMaria/Scene_2_AB1823/NIRSpec_TA_Sim_AB1823 first NRS/postage")
path_scene2_slow_nonoise = os.path.abspath("PFforMaria/Scene_2_AB1823/NIRSpec_TA_Sim_AB1823 first NRS no_noise/postage")
path_scene2_rapid = os.path.abspath("PFforMaria/Scene_2_AB1823/NIRSpec_TA_Sim_AB1823 first NRSRAPID/postage")
path_scene2_rapid_nonoise = os.path.abspath("PFforMaria/Scene_2_AB1823/NIRSpec_TA_Sim_AB1823 first NRSRAPID no_noise/postage")
#                      0                   1                     2                      3            
paths_list = [path_scene1_slow, path_scene1_slow_nonoise, path_scene1_rapid, path_scene1_rapid_nonoise, 
             path_scene2_slow, path_scene2_slow_nonoise, path_scene2_rapid, path_scene2_rapid_nonoise]
#                      4                   5                     6                      7            

##########################################################################################

# Set test parameters
path_number = 0   # select 0 through 7 from paths_list above
checkbox_size = 3
background_method = None#"fix"   # Select either 'fractional', 'fixed', or None   
xwidth_list = [3, 5, 7]   # Number of rows of the centroid region
ywidth_list = [3, 5, 7]   # Number of columns of the centroid region
max_iter = 50
threshold = 1e-5
detector = 491
save_text_file = False
determine_moments = False   # Want to determine 2nd and 3rd moments?
display_master_img = False   # Want to see the combined ramped images for every star?
debug = False   # see all debug messages (i.e. values of all calculations)

##########################################################################################

##### FUNCTIONS
def run_recursive_centroids(master_img_bgcorr, background, xwidth_list, ywidth_list, checkbox_size, max_iter, 
                            threshold, debug, save_text_file, determine_moments):    
    # Obtain and display the combined FITS image that combines all frames into one image
    psf = tu.readimage(master_img_bgcorr)
    if display_master_img: 
        tu.display_ns_psf(psf)
    # Test checkbox piece
    cb_centroid_list = []
    for xwidth, ywidth in zip(xwidth_list, ywidth_list):
        print ("Testing centroid width: ", checkbox_size)
        print ("     xwidth = ", xwidth, "  ywidth = ", ywidth)
        cb_cen, cb_hw = jtl.checkbox_2D(psf, checkbox_size, xwidth, ywidth, debug=debug)
        print ('Got coarse location for checkbox_size {} \n'.format(checkbox_size))
        # Checkbox center, in base 1
        print('Checkbox Output:')
        print('Checkbox center: [{}, {}]'.format(cb_cen[0], cb_cen[1]))
        print('Checkbox halfwidths: xhw: {}, yhw: {}'.format(cb_hw[0], cb_hw[1]))
        print()
        # Calculate the centroid based on the checkbox region calculated above
        cb_centroid, cb_sum = jtl.centroid_2D(psf, cb_cen, cb_hw, max_iter=max_iter, threshold=threshold, debug=debug)
        cb_centroid_list.append(cb_centroid)
        print('Final sum: ', cb_sum)
        print('cb_centroid: ', cb_centroid)
        print()
    # Write output into text file
    bg = background
    line1 = "{:<12} {:<14} {:<18} {:<14} {:<18} {:<14} {:<18}\n".format(bg, cb_centroid_list[0][0], cb_centroid_list[0][1],
                                                              cb_centroid_list[1][0], cb_centroid_list[1][1],
                                                              cb_centroid_list[2][0], cb_centroid_list[2][1])
    if save_text_file:
        f = open(output_file, "a")
        f.write(line1+"\n")
        f.close()
    # Find the 2nd and 3rd moments
    if determine_moments:
        x_mom, y_mom = jtl.find2D_higher_moments(psf, cb_centroid, cb_hw, cb_sum)
        print('Higher moments(2nd, 3rd):')
        print('x_moments: ', x_mom)
        print('y moments: ', y_mom)
        print('---------------------------------------------------------------')
        print()
    return line1

##### CODE

# Background cases
bg_frac, bg_value = None, None   # for the None case
bg_choice = "_bgNone"
background = 0.0
if background_method is not None and "frac" in background_method:
    fractional_background_list = [x*0.1 for x in range(11)]
    bg_choice = "_bgFrac"
elif background_method is not None and "fix" in background_method:
    bg_value = 0.0
    bg_choice = "_bgFixed"
    background = bg_value

# Set path of test directory
dir2test = paths_list[path_number]

# Set the name for the output file
if path_number == 0:
    case = "scene1_slow_real"
elif path_number == 1:
    case = "scene1_slow_nonoise"
elif path_number == 2:
    case = "scene1_rapid_real"
elif path_number == 3:
    case = "scene1_rapid_nonoise"
elif path_number == 4:
    case = "scene2_slow_real"
elif path_number == 5:
    case = "scene2_slow_nonoise"
elif path_number == 6:
    case = "scene2_rapid_real"
elif path_number == 7:
    case = "scene2_rapid_nonoise"
output_file = "TA_testcases_for_"+case+bg_choice+".txt"
line0a = "{:<15} {:<16} {:>30} {:>32}".format("Background", "Centroid width: 3", "5", "7")
line0b = "{:>20} {:>14} {:>18} {:>14} {:>18} {:>14}".format("x", "y", "x", "y", "x", "y")
if save_text_file:
    f = open(output_file, "w+")
    f.write(line0a+"\n")
    f.write(line0b+"\n")
    f.close()

# Stars of detector 492
stars_492 = range(1, 101)
# Stars of detector 491
stars_491 = [x+100 for x in range(1, 101)]
# Select appropriate set of stars to test according to chosen detector
if detector == 491:
    stars_detector = stars_491
elif detector == 492:
    stars_detector = stars_492

# Check if the directory path exists
dir_exist = os.path.isdir(dir2test)
if dir_exist == False:
    print ("The directory: ", dir2test, "\n    does NOT exist. Exiting the script.")
    exit()

# Read the star parameters file to compare results
star_param_txt = os.path.join(dir2test,"star parameters.txt")
benchmark_data = np.loadtxt(star_param_txt, skiprows=2, unpack=True)
bench_star, quadrant, star_in_quad, x_491, y_491, x_492, y_492, V2, V3 = benchmark_data
# The x and y positions in the 491 and 492 detectors are in coordinates of the whole detector
# transform to 32 by 32 coordinates:
bench_x = np.floor(x_491) - 16.0
bench_y = np.floor(y_491) - 16.0
if detector == 492:
    bench_x = np.floor(x_492) - 16.0
    bench_y = np.floor(y_492) - 16.0
# Corrections for offsets in positions (see section 2.5 of Technical Notes in Documentation directory)
offset_491 = (-0.086, -0.077)
offset_942 = (0.086, 0.077)

# Start the loop in the given directory
dir_stars = glob(os.path.join(dir2test,"postageout_star_*.fits"))   # get all star fits files in that directory
for star in dir_stars:
    # Test stars of detector of choice
    for st in stars_detector:
        if str(st)+" quad_       " in star:
            print ("Will test star in directory: \n     ", dir2test)
            print ("Star: ", os.path.basename(star))
            # Make sure the file actually exists
            star_exists = os.path.isfile(star)
            if star_exists == False:
                print ("The file: ", star, "\n    does NOT exist. Exiting the script.")
                exit() 
            
            # Read FITS image 
            #hdr = fits.getheader(star, 0)
            #print("** HEADER:", hdr)
            master_img = fits.getdata(star, 0)
            print ('Master image shape: ', np.shape(master_img))
            # Do background correction on each of 3 ramp images
            if background_method is not None and "frac" in background_method:
                # If fractional method is selected, loop over backgrounds from 0.0 to 1.0 in increments of 0.1
                for bg_frac in fractional_background_list:
                    master_img_bgcorr = jtl.bg_correction(master_img, bg_method=background_method, 
                                                          bg_value=bg_value, bg_frac=bg_frac)
                    line1 = run_recursive_centroids(master_img_bgcorr, bg_frac, xwidth_list, ywidth_list, checkbox_size,
                                                    max_iter, threshold, debug, save_text_file, determine_moments)
            else:
                master_img_bgcorr = jtl.bg_correction(master_img, bg_method=background_method, 
                                                      bg_value=bg_value, bg_frac=bg_frac)
                line1 = run_recursive_centroids(master_img_bgcorr, background, xwidth_list, ywidth_list, checkbox_size, 
                                                max_iter, threshold, debug, save_text_file, determine_moments)
            print(line0a)
            print(line0b)
            print(line1) 
            #exit()
            
print ("Recursive test script finished. \n")