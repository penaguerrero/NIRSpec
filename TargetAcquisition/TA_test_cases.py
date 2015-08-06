from __future__ import print_function, division
from astropy.io import fits
from glob import glob
import numpy as np
import matplotlib.pyplot as plt
import time
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
path_scene1_slow_shifted = "PFforMaria/Scene_1_AB23/NIRSpec_TA_Sim_AB23 shifted NRS/postage"
path_scene1_slow_shifted_nonoise = os.path.abspath("PFforMaria/Scene_1_AB23/NIRSpec_TA_Sim_AB23 shifted NRS no_noise/postage")
path_scene1_rapid_shifted = os.path.abspath("PFforMaria/Scene_1_AB23/NIRSpec_TA_Sim_AB23 shifted NRSRAPID/postage")
path_scene1_rapid_shifted_nonoise = os.path.abspath("PFforMaria/Scene_1_AB23/NIRSpec_TA_Sim_AB23 shifted NRS no_noise/postage")
path_scene2_slow = os.path.abspath("PFforMaria/Scene_2_AB1823/NIRSpec_TA_Sim_AB1823 first NRS/postage")
path_scene2_slow_nonoise = os.path.abspath("PFforMaria/Scene_2_AB1823/NIRSpec_TA_Sim_AB1823 first NRS no_noise/postage")
path_scene2_rapid = os.path.abspath("PFforMaria/Scene_2_AB1823/NIRSpec_TA_Sim_AB1823 first NRSRAPID/postage")
path_scene2_rapid_nonoise = os.path.abspath("PFforMaria/Scene_2_AB1823/NIRSpec_TA_Sim_AB1823 first NRSRAPID no_noise/postage")
path_scene2_slow_shifted = os.path.abspath("PFforMaria/Scene_2_AB1823/NIRSpec_TA_Sim_AB1823 shifted NRS/postage")
path_scene2_slow_shifted_nonoise = os.path.abspath("PFforMaria/Scene_2_AB1823/NIRSpec_TA_Sim_AB1823 shifted NRS no_noise/postage")
path_scene2_rapid_shifted = os.path.abspath("PFforMaria/Scene_2_AB1823/NIRSpec_TA_Sim_AB1823 shifted NRSRAPID/postage")
path_scene2_rapid_shifted_nonoise = os.path.abspath("PFforMaria/Scene_2_AB1823/NIRSpec_TA_Sim_AB1823 shifted NRSRAPID no_noise/postage")
#                      0                   1                     2                      3            
paths_list = [path_scene1_slow, path_scene1_slow_nonoise, path_scene1_rapid, path_scene1_rapid_nonoise,
#                      4                        5                                   6                           7            
              path_scene1_slow_shifted, path_scene1_slow_shifted_nonoise, path_scene1_rapid_shifted, path_scene1_rapid_shifted_nonoise, 
#                      8                   9                     10                      11            
              path_scene2_slow, path_scene2_slow_nonoise, path_scene2_rapid, path_scene2_rapid_nonoise,
#                      12                       13                                  14                          15            
              path_scene2_slow_shifted, path_scene2_slow_shifted_nonoise, path_scene2_rapid_shifted, path_scene2_rapid_shifted_nonoise]


###########################################################################################################


# Set test parameters
path_number = 15             # Select 0 through 7 from paths_list above OR  
single_star = False          # If only want to test one star set to True and type the path for a single star 
single_star_path = paths_list[10]+'/postageout_star_     103 quad_       3 quad_star        3.fits'
display_master_img = False  # Want to see the combined ramped images for every star?
debug = False               # see all debug messages (i.e. values of all calculations)
centroid_in_full_detector = False   # Want results in full detector coordinates or 32x32? 
checkbox_size = 3
background_method = 'frac'  # Select either 'fractional', 'fixed', or None   
xwidth_list = [3, 5, 7]     # Number of rows of the centroid region
ywidth_list = [3, 5, 7]     # Number of columns of the centroid region
max_iter = 50
threshold = 1e-5
detector = 491
determine_moments = False   # Want to determine 2nd and 3rd moments?
save_text_file = False
just_read_text_file = True  # skip the for loop to the plotting part
# make plot of magnitude (in x) versus radial offset distance (in y) for Scene2 
show_plot = True
save_plot = False
plot_type = '.jpg'


###########################################################################################################

if single_star:
    #print('got here')
    save_text_file = False
    display_master_img = True   # Want to see the combined ramped images for every star?
    debug = True   # see all debug messages (i.e. values of all calculations)
    vlim = ()#(0.1, 5)
else:
    vlim = ()

##### FUNCTIONS
def run_recursive_centroids(psf, background, xwidth_list, ywidth_list, checkbox_size, max_iter, 
                            threshold, determine_moments, debug, vlim=()):   
    """
    Determine the centroid location given the that the background is already subtracted. 
    """ 
    # Display the combined FITS image that combines all frames into one image
    if display_master_img: 
        tu.display_ns_psf(psf, vlim=vlim)
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
        #raw_input()
    # Find the 2nd and 3rd moments
    if determine_moments:
        x_mom, y_mom = jtl.find2D_higher_moments(psf, cb_centroid, cb_hw, cb_sum)
        print('Higher moments(2nd, 3rd):')
        print('x_moments: ', x_mom)
        print('y moments: ', y_mom)
        print('---------------------------------------------------------------')
        print()
    return cb_centroid_list


def transform2fulldetector(detector, centroid_in_full_detector, cb_centroid_list, ESA_center, true_center):
    """
    Transform centroid coordinates into full detector coordinates.
    
    Keyword arguments:
    detector                   -- Either 491 or 492
    centroid_in_full_detector  -- Resuling coordinates in terms of full detector, True or False
    cb_centroid_list           -- Centroid determined by target acquisition (TA) algorithm in terms of 32 by 32 pixels
    ESA_center                 -- Centroid determined with the ESA version of TA algorithm in terms of full detector
    true_center                -- Actual (true) position of star in terms of full detector  
    
    Output(s):
    cb_centroid_list_fulldetector  -- List of centroid locations determined with the TA algorithm in 
                                      terms of full detector. List is for positions determined with
                                      3, 5, and 7 checkbox sizes. 
    
    """
    # Corrections for offsets in positions (see section 2.5 of Technical Notes in Documentation directory)
    offset_491 = (-0.086, -0.077)
    offset_492 = (0.086, 0.077)
    if detector == 491:
        corrected_x = true_center[0] + offset_491[0]
        corrected_y = true_center[1] + offset_491[1]
    elif detector == 492:
        corrected_x = true_center[0] + offset_492[0]
        corrected_y = true_center[1] + offset_492[1]
        
    # Get the lower left corner coordinates in terms of full detector. We subtract 16.0 because indexing
    # from centroid function starts with 1
    loleft_x = np.floor(corrected_x) - 16.0
    loleft_y = np.floor(corrected_y) - 16.0

    if centroid_in_full_detector:    
        # Add lower left corner to centroid location to get it in terms of full detector
        cb_centroid_list_fulldetector = []
        for centroid_location in cb_centroid_list:
            centroid_fulldetector_x = centroid_location[0] + loleft_x
            centroid_fulldetector_y = centroid_location[1] + loleft_y
            centroid_fulldetector = [centroid_fulldetector_x, centroid_fulldetector_y]
            cb_centroid_list_fulldetector.append(centroid_fulldetector)
        corr_cb_centroid_list = cb_centroid_list_fulldetector
        corr_true_center_centroid = true_center
    else:
        corr_cb_centroid_list = cb_centroid_list
        # Add lower left corner to centroid location to get it in terms of full detector
        corr_true_center_x = corrected_x - loleft_x
        corr_true_center_y = corrected_y - loleft_y
        true_center_centroid_32x32 = [corr_true_center_x, corr_true_center_y]
        corr_true_center_centroid = true_center_centroid_32x32

    # Determine difference between center locations
    differences_true_TA = []
    d3_x = corr_true_center_centroid[0] - corr_cb_centroid_list[0][0]
    d5_x = corr_true_center_centroid[0] - corr_cb_centroid_list[1][0]
    d7_x = corr_true_center_centroid[0] - corr_cb_centroid_list[2][0]
    d3_y = corr_true_center_centroid[1] - corr_cb_centroid_list[0][1]
    d5_y = corr_true_center_centroid[1] - corr_cb_centroid_list[1][1]
    d7_y = corr_true_center_centroid[1] - corr_cb_centroid_list[2][1]
    d3 = [d3_x, d3_y]
    d5 = [d5_x, d5_y]
    d7 = [d7_x, d7_y]
    diffs = [d3, d5, d7]
    differences_true_TA.append(diffs)
    return corr_true_center_centroid, corr_cb_centroid_list, differences_true_TA


def write2file(save_text_file, output_file, st, bg, corr_cb_centroid_list, corr_true_center_centroid, differences_true_TA):
    line1 = "{:<5} {:<10} {:<14} {:<16} {:<14} {:<16} {:<14} {:<16} {:<12} {:<14} {:<20} {:<22} {:<20} {:<22} {:<20} {:<22}\n".format(st, bg, 
                                                    corr_cb_centroid_list[0][0], corr_cb_centroid_list[0][1],
                                                    corr_cb_centroid_list[1][0], corr_cb_centroid_list[1][1],
                                                    corr_cb_centroid_list[2][0], corr_cb_centroid_list[2][1],
                                                    corr_true_center_centroid[0], corr_true_center_centroid[1],
                                                    differences_true_TA[0][0][0], differences_true_TA[0][0][1],
                                                    differences_true_TA[0][1][0], differences_true_TA[0][1][1],
                                                    differences_true_TA[0][2][0], differences_true_TA[0][2][1])
    if save_text_file:
        f = open(output_file, "a")
        f.write(line1)
        f.close()
    print(line0a)
    print(line0b)
    print(line1) 


def read_listfile(list_file_name):    
    """ This function reads the fits table that contains the flux and converts to magnitude for the 
    simulated stars. """
    listfiledata = fits.getdata(list_file_name)
    star_number, xpos, ypos, orient, factor = np.array([]), np.array([]), np.array([]), np.array([]), np.array([])  
    for row in listfiledata:
        star_number = np.append(star_number, row[0]) 
        xpos = np.append(xpos, row[1]) 
        ypos = np.append(ypos, row[2])
        orient = np.append(orient, row[3])
        factor = np.append(factor, row[4])
    # convert the flux into magnitude (factor=1.0 is equivalent to magnitude=18.0)
    mag = 2.5*np.log10(factor) + 18.0
    return star_number, xpos, ypos, mag
    

def get_fracdata(offsets):
    """ This function gets arrays for each fractional background for the same star. """
    frac003x, frac005x, frac007x = [], [], []
    frac003y, frac005y, frac007y = [], [], []
    frac013x, frac015x, frac017x = [], [], []
    frac013y, frac015y, frac017y = [], [], []
    frac023x, frac025x, frac027x = [], [], []
    frac023y, frac025y, frac027y = [], [], []
    frac033x, frac035x, frac037x = [], [], []
    frac033y, frac035y, frac037y = [], [], []
    frac043x, frac045x, frac047x = [], [], []
    frac043y, frac045y, frac047y = [], [], []
    frac053x, frac055x, frac057x = [], [], []
    frac053y, frac055y, frac057y = [], [], []
    frac063x, frac065x, frac067x = [], [], []
    frac063y, frac065y, frac067y = [], [], []
    frac073x, frac075x, frac077x = [], [], []
    frac073y, frac075y, frac077y = [], [], []
    frac083x, frac085x, frac087x = [], [], []
    frac083y, frac085y, frac087y = [], [], []
    frac093x, frac095x, frac097x = [], [], []
    frac093y, frac095y, frac097y = [], [], []
    frac103x, frac105x, frac107x = [], [], []
    frac103y, frac105y, frac107y = [], [], []
    i0,i1,i2,i3,i4,i5,i6,i7,i8,i9,i10 = 0,1,2,3,4,5,6,7,8,9,10
    for i, _ in enumerate(offsets[0]):
        #row = [offsets[0][i], offsets[1][i], offsets[2][i], offsets[3][i], offsets[4][i], offsets[5][i]]
        #row = np.array(row).reshape(1,6)
        if i == i0:
            frac003x.append(offsets[0][i])
            frac003y.append(offsets[1][i])
            frac005x.append(offsets[2][i])
            frac005y.append(offsets[3][i])
            frac007x.append(offsets[4][i])
            frac007y.append(offsets[5][i])
            i0 = i0+11
        if i == i1:
            frac013x.append(offsets[0][i])
            frac013y.append(offsets[1][i])
            frac015x.append(offsets[2][i])
            frac015y.append(offsets[3][i])
            frac017x.append(offsets[4][i])
            frac017y.append(offsets[5][i])
            i1 = i1+11
        if i == i2:
            frac023x.append(offsets[0][i])
            frac023y.append(offsets[1][i])
            frac025x.append(offsets[2][i])
            frac025y.append(offsets[3][i])
            frac027x.append(offsets[4][i])
            frac027y.append(offsets[5][i])
            i2 = i2+11
        if i == i3:
            frac033x.append(offsets[0][i])
            frac033y.append(offsets[1][i])
            frac035x.append(offsets[2][i])
            frac035y.append(offsets[3][i])
            frac037x.append(offsets[4][i])
            frac037y.append(offsets[5][i])
            i3 = i3+11
        if i == i4:
            frac043x.append(offsets[0][i])
            frac043y.append(offsets[1][i])
            frac045x.append(offsets[2][i])
            frac045y.append(offsets[3][i])
            frac047x.append(offsets[4][i])
            frac047y.append(offsets[5][i])
            i4 = i4+11
        if i == i5:
            frac053x.append(offsets[0][i])
            frac053y.append(offsets[1][i])
            frac055x.append(offsets[2][i])
            frac055y.append(offsets[3][i])
            frac057x.append(offsets[4][i])
            frac057y.append(offsets[5][i])
            i5 = i5+11
        if i == i6:
            frac063x.append(offsets[0][i])
            frac063y.append(offsets[1][i])
            frac065x.append(offsets[2][i])
            frac065y.append(offsets[3][i])
            frac067x.append(offsets[4][i])
            frac067y.append(offsets[5][i])
            i6 = i6+11
        if i == i7:
            frac073x.append(offsets[0][i])
            frac073y.append(offsets[1][i])
            frac075x.append(offsets[2][i])
            frac075y.append(offsets[3][i])
            frac077x.append(offsets[4][i])
            frac077y.append(offsets[5][i])
            i7 = i7+11
        if i == i8:
            frac083x.append(offsets[0][i])
            frac083y.append(offsets[1][i])
            frac085x.append(offsets[2][i])
            frac085y.append(offsets[3][i])
            frac087x.append(offsets[4][i])
            frac087y.append(offsets[5][i])
            i8 = i8+11
        if i == i9:
            frac093x.append(offsets[0][i])
            frac093y.append(offsets[1][i])
            frac095x.append(offsets[2][i])
            frac095y.append(offsets[3][i])
            frac097x.append(offsets[4][i])
            frac097y.append(offsets[5][i])
            i9 = i9+11
        if i == i10:
            frac103x.append(offsets[0][i])
            frac103y.append(offsets[1][i])
            frac105x.append(offsets[2][i])
            frac105y.append(offsets[3][i])
            frac107x.append(offsets[4][i])
            frac107y.append(offsets[5][i])
            i10 = i10+11
    frac00 = np.array([frac003x, frac003y, frac005x, frac005y, frac007x, frac007y])
    frac01 = np.array([frac013x, frac013y, frac015x, frac015y, frac017x, frac017y])
    frac02 = np.array([frac023x, frac023y, frac025x, frac025y, frac027x, frac027y])
    frac03 = np.array([frac033x, frac033y, frac035x, frac035y, frac037x, frac037y])
    frac04 = np.array([frac043x, frac043y, frac045x, frac045y, frac047x, frac047y])
    frac05 = np.array([frac053x, frac053y, frac055x, frac055y, frac057x, frac057y])
    frac06 = np.array([frac063x, frac063y, frac065x, frac065y, frac067x, frac067y])
    frac07 = np.array([frac073x, frac073y, frac075x, frac075y, frac077x, frac077y])
    frac08 = np.array([frac083x, frac083y, frac085x, frac085y, frac087x, frac087y])
    frac09 = np.array([frac093x, frac093y, frac095x, frac095y, frac097x, frac097y])
    frac10 = np.array([frac103x, frac103y, frac105x, frac105y, frac107x, frac107y])
    return frac00, frac01, frac02, frac03, frac04, frac05, frac06, frac07, frac08, frac09, frac10


##### CODE

# start the timer to compute the whole running time
start_time = time.time()

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
    case = "scene1_slow_real_shifted"
elif path_number == 5:
    case = "scene1_slow_nonoise_shifted"
elif path_number == 6:
    case = "scene1_rapid_real_shifted"
elif path_number == 7:
    case = "scene1_rapid_nonoise_shifted"
elif path_number == 8:
    case = "scene2_slow_real"
elif path_number == 9:
    case = "scene2_slow_nonoise"
elif path_number == 10:
    case = "scene2_rapid_real"
elif path_number == 11:
    case = "scene2_rapid_nonoise"
elif path_number == 12:
    case = "scene2_slow_real_shifted"
elif path_number == 13:
    case = "scene2_slow_nonoise_shifted"
elif path_number == 14:
    case = "scene2_rapid_real_shifted"
elif path_number == 15:
    case = "scene2_rapid_nonoise_shifted"
outpul_file_path = "PFforMaria/Resulting_centroid_txt_files/"
output_file = os.path.join(outpul_file_path, "TA_testcases_for_"+case+bg_choice+".txt")
line0 = "Centroid indexing starting at 1 !"
line0a = "{:<5} {:<15} {:<16} {:>23} {:>32} {:>33} {:>35} {:>38} {:>43}".format("Star", "Background", 
                                                                  "Centroid width: 3", "5", "7", "TruePositions", 
                                                                  "Difference with: checkbox 3", "checkbox 5", "checkbox 7")
line0b = "{:>25} {:>12} {:>16} {:>14} {:>16} {:>14} {:>16} {:>14} {:>16} {:>14} {:>26} {:>14} {:>28} {:>14}".format(
                                                                       "x", "y", "x", "y", "x", "y", "TrueX", "TrueY", 
                                                                       "x", "y", "x", "y", "x", "y")
if save_text_file:
    f = open(output_file, "w+")
    f.write(line0+"\n")
    f.write(line0a+"\n")
    f.write(line0b+"\n")
    f.close()

# Stars of detector 492
stars_492 = range(1, 101)
# Stars of detector 491
stars_491 = [x+100 for x in range(1, 101)]

# Check if the directory path exists
dir_exist = os.path.isdir(dir2test)
if dir_exist == False:
    print ("The directory: ", dir2test, "\n    does NOT exist. Exiting the script.")
    exit()

# Read the star parameters file to compare results
star_param_txt = os.path.join(dir2test,"star parameters.txt")
benchmark_data = np.loadtxt(star_param_txt, skiprows=2, unpack=True)
bench_star, quadrant, star_in_quad, x_491, y_491, x_492, y_492, V2, V3 = benchmark_data
# Select appropriate set of stars to test according to chosen detector
if detector == 491:
    stars_detector = stars_491
    true_x = x_491
    true_y = y_491
elif detector == 492:
    stars_detector = stars_492
    true_x = x_492
    true_y = y_492

if just_read_text_file != True:
    # Start the loop in the given directory
    dir_stars = glob(os.path.join(dir2test,"postageout_star_*.fits"))   # get all star fits files in that directory
    for star in dir_stars:
        if single_star:
            star = single_star_path
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
                
                # Obtain real star position
                idx_star = stars_detector.index(st)
                true_center = [true_x[idx_star], true_y[idx_star]]
                
                # Read FITS image 
                #hdr = fits.getheader(star, 0)
                #print("** HEADER:", hdr)
                master_img = fits.getdata(star, 0)
                print ('Master image shape: ', np.shape(master_img))
                # Do background correction on each of 3 ramp images
                if background_method is not None and "frac" in background_method:
                    # If fractional method is selected, loop over backgrounds from 0.0 to 1.0 in increments of 0.1
                    for bg_frac in fractional_background_list:
                        print ("* Using fractional background value of: ", bg_frac)
                        master_img_bgcorr = jtl.bg_correction(master_img, bg_method=background_method, 
                                                              bg_value=bg_value, bg_frac=bg_frac, debug=debug)
                        # Obtain the combined FITS image that combines all frames into one image AND
                        # check if all image is zeros, take the image that still has a max value
                        psf = tu.readimage(master_img_bgcorr, debug=debug)
                        master_img_bgcorr_max = psf.max()
                        while master_img_bgcorr_max == 0.0:
                            print('  IMPORTANT WARNING!!! Combined ramped images have a max of 0.0 with bg_frac=', bg_frac)
                            bg_frac = bg_frac - 0.02
                            if bg_frac < 0.0:   # prevent an infinite loop
                                print ('   ERROR - Cannot subtract from  bg_frac < 0.0 ...')
                                bg_frac = 0.0
                                break
                            print('       *** Setting  NEW  bg_frac = ', bg_frac)
                            master_img_bgcorr = jtl.bg_correction(master_img, bg_method=background_method, 
                                                                  bg_value=bg_value, bg_frac=bg_frac, debug=debug)
                            psf = tu.readimage(master_img_bgcorr, debug=debug)
                            master_img_bgcorr_max = psf.max()
                        cb_centroid_list = run_recursive_centroids(psf, bg_frac, xwidth_list, ywidth_list, 
                                                                   checkbox_size, max_iter, threshold, 
                                                                   determine_moments, debug, vlim=vlim)
                        # Transform to full detector coordinates in order to compare with real centers
                        ESA_center = [0,0]
                        corr_true_center_centroid, corr_cb_centroid_list, differences_true_TA = transform2fulldetector(detector, 
                                                                                                      centroid_in_full_detector,
                                                                                                      cb_centroid_list, ESA_center, 
                                                                                                      true_center)
                        # Write output into text file
                        bg = bg_frac
                        write2file(save_text_file, output_file, st, bg, corr_cb_centroid_list, corr_true_center_centroid, differences_true_TA)
                        #raw_input()
                else:
                    master_img_bgcorr = jtl.bg_correction(master_img, bg_method=background_method, 
                                                          bg_value=bg_value, bg_frac=bg_frac)
                    # Obtain the combined FITS image that combines all frames into one image AND
                    # check if all image is zeros, take the image that still has a max value
                    psf = tu.readimage(master_img_bgcorr, debug=debug)                
                    cb_centroid_list = run_recursive_centroids(psf, bg_frac, xwidth_list, ywidth_list, 
                                                               checkbox_size, max_iter, threshold, 
                                                               determine_moments, debug, vlim=vlim)
                    # Transform to full detector coordinates in order to compare with real centers
                    ESA_center = [0,0]
                    corr_true_center_centroid, corr_cb_centroid_list, differences_true_TA = transform2fulldetector(detector, 
                                                                                                  centroid_in_full_detector,
                                                                                                  cb_centroid_list, ESA_center, 
                                                                                                  true_center)
                    # Write output into text file
                    bg = background
                    write2file(save_text_file, output_file, st, bg, corr_cb_centroid_list, corr_true_center_centroid, differences_true_TA) 
                
                if single_star:
                    print ("Recursive test script finished. \n")
                    exit()

# Make the plot of magnitude (in x) versus radial offset distance (in y) for Scene2
if "scene2" in case:
    # Read the text file just written to get the offsets from the "real" positions of the fake stars
    offsets = np.loadtxt(output_file, skiprows=3, usecols=(10,11,12,13,14,15), unpack=True)
    path2listfile = "PFforMaria/Scene_2_AB1823"
    list_file = "simuTA20150528-F140X-S50-K-AB18to23.list"
    lf = os.path.join(path2listfile,list_file)
    star_number, xpos, ypos, mag = read_listfile(lf)
    if detector == 491:   # slice from star 101 to 200
        star_number, xpos, ypos, mag = star_number[100:], xpos[100:], ypos[100:], mag[100:]
    elif detector == 492:   # slice from star 1 to 100
        star_number, xpos, ypos, mag = star_number[:100], xpos[:100], ypos[:100:], mag[:100]
    bg_method = background_method
    if background_method is None:   # convert the None value to string
        bg_method = 'None'
    if 'frac' not in bg_method:
        fig1 = plt.figure(1, figsize=(12, 10))
        plt.title(case+'_BG'+bg_method)
        plt.xlabel('Magnitude')
        plt.ylabel('Radial offset in Y')
        plt.hlines(0.0, 18.0, 23.0, colors='k', linestyles='dashed')
        plt.plot(mag, offsets[1], 'bo', ms=8, alpha=0.7, label='Checkbox=3')
        plt.plot(mag, offsets[3], 'go', ms=8, alpha=0.7, label='Checkbox=5')
        plt.plot(mag, offsets[5], 'ro', ms=8, alpha=0.7, label='Checkbox=7')
        plt.legend(loc='upper left')
        if save_plot:
            if background_method is None:
                bg = 'None_'
            else:
                bg = 'fix_'
            destination = os.path.abspath("PFforMaria/plots/MagVsYoffset_"+bg+case+plot_type)
            plt.savefig(destination)
            print ("\n Plot saved: ", destination)
        if show_plot:
            plt.show()
    else:
        print ('frac!')
        print ("Reading text file: ",  output_file)
        frac00, frac01, frac02, frac03, frac04, frac05, frac06, frac07, frac08, frac09, frac10 = get_fracdata(offsets)
        fig1 = plt.figure(1, figsize=(12, 10))
        fig1.subplots_adjust(hspace=0.10)
        ax1 = fig1.add_subplot(311)
        ax1.set_title(case+"_BGfrac")
        ax1.set_xlabel('Magnitude')
        ax1.set_ylabel('Radial offset in Y: Checkbox=3')
        plt.hlines(0.0, 18.0, 23.0, colors='k', linestyles='dashed')
        ax1.plot(mag, frac00[1], 'bo', ms=8, alpha=0.7, label='bg_frac=0.0')
        ax1.plot(mag, frac01[1], 'ro', ms=8, alpha=0.7, label='bg_frac=0.1')
        ax1.plot(mag, frac02[1], 'mo', ms=8, alpha=0.7, label='bg_frac=0.2')
        ax1.plot(mag, frac03[1], 'go', ms=5, alpha=0.7, label='bg_frac=0.3')
        ax1.plot(mag, frac04[1], 'ko', ms=8, alpha=0.7, label='bg_frac=0.4')
        ax1.plot(mag, frac05[1], 'yo', ms=8, alpha=0.7, label='bg_frac=0.5')
        ax1.plot(mag, frac06[1], 'co', ms=8, alpha=0.7, label='bg_frac=0.6')
        ax1.plot(mag, frac07[1], 'b+', ms=10, alpha=0.7, label='bg_frac=0.7')
        ax1.plot(mag, frac08[1], 'r+', ms=8, alpha=0.7, label='bg_frac=0.8')
        ax1.plot(mag, frac09[1], 'm+', ms=5, alpha=0.7, label='bg_frac=0.9')
        ax1.plot(mag, frac10[1], 'k+', ms=5, alpha=0.7, label='bg_frac=1.0')
        # Shrink current axis by 10%
        box = ax1.get_position()
        ax1.set_position([box.x0, box.y0, box.width * 0.9, box.height])
        ax2 = fig1.add_subplot(312)
        ax2.set_xlabel('Magnitude')
        ax2.set_ylabel('Radial offset in Y: Checkbox=5')
        plt.hlines(0.0, 18.0, 23.0, colors='k', linestyles='dashed')
        ax2.plot(mag, frac00[3], 'bo', ms=8, alpha=0.7, label='bg_frac=0.0')
        ax2.plot(mag, frac01[3], 'ro', ms=8, alpha=0.7, label='bg_frac=0.1')
        ax2.plot(mag, frac02[3], 'mo', ms=8, alpha=0.7, label='bg_frac=0.2')
        ax2.plot(mag, frac03[3], 'go', ms=5, alpha=0.7, label='bg_frac=0.3')
        ax2.plot(mag, frac04[3], 'ko', ms=8, alpha=0.7, label='bg_frac=0.4')
        ax2.plot(mag, frac05[3], 'yo', ms=8, alpha=0.7, label='bg_frac=0.5')
        ax2.plot(mag, frac06[3], 'co', ms=8, alpha=0.7, label='bg_frac=0.6')
        ax2.plot(mag, frac07[3], 'b+', ms=10, alpha=0.7, label='bg_frac=0.7')
        ax2.plot(mag, frac08[3], 'r+', ms=8, alpha=0.7, label='bg_frac=0.8')
        ax2.plot(mag, frac09[3], 'm+', ms=5, alpha=0.7, label='bg_frac=0.9')
        ax2.plot(mag, frac10[3], 'k+', ms=5, alpha=0.7, label='bg_frac=1.0')
        # Shrink current axis by 10%
        box = ax2.get_position()
        ax2.set_position([box.x0, box.y0, box.width * 0.9, box.height])
        # put legend out of the plot box
        ax2.legend(loc='center left', bbox_to_anchor=(1, 0.5))            
        ax3 = fig1.add_subplot(313)
        ax3.set_xlabel('Magnitude')
        ax3.set_ylabel('Radial offset in Y: Checkbox=7')
        plt.hlines(0.0, 18.0, 23.0, colors='k', linestyles='dashed')
        ax3.plot(mag, frac00[5], 'bo', ms=8, alpha=0.7, label='bg_frac=0.0')
        ax3.plot(mag, frac01[5], 'ro', ms=8, alpha=0.7, label='bg_frac=0.1')
        ax3.plot(mag, frac02[5], 'mo', ms=8, alpha=0.7, label='bg_frac=0.2')
        ax3.plot(mag, frac03[5], 'go', ms=5, alpha=0.7, label='bg_frac=0.3')
        ax3.plot(mag, frac04[5], 'ko', ms=8, alpha=0.7, label='bg_frac=0.4')
        ax3.plot(mag, frac05[5], 'yo', ms=8, alpha=0.7, label='bg_frac=0.5')
        ax3.plot(mag, frac06[5], 'co', ms=8, alpha=0.7, label='bg_frac=0.6')
        ax3.plot(mag, frac07[5], 'b+', ms=10, alpha=0.7, label='bg_frac=0.7')
        ax3.plot(mag, frac08[5], 'r+', ms=8, alpha=0.7, label='bg_frac=0.8')
        ax3.plot(mag, frac09[5], 'm+', ms=5, alpha=0.7, label='bg_frac=0.9')
        ax3.plot(mag, frac10[5], 'k+', ms=5, alpha=0.7, label='bg_frac=1.0')
        # Shrink current axis by 10%
        box = ax3.get_position()
        ax3.set_position([box.x0, box.y0, box.width * 0.9, box.height])
        if save_plot:
            destination = os.path.abspath("PFforMaria/plots/MagVsYoffset_frac_"+case+plot_type)
            plt.savefig(destination)
            print ("\n Plot saved: ", destination)
        if show_plot:
            plt.show()
        

if not single_star:
    print ("\n Centroids and differences were written into: \n  {}".format(output_file))
print ("\n Recursive test script finished. Took  %s  seconds to finish. \n" % ((time.time() - start_time)) )
            