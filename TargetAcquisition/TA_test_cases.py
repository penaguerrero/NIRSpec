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

# other code
import testing_functions as tf

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
path_number = 15             # Select 0 through 15 from paths_list above OR  
save_text_file = True
just_read_text_file = False  # skip the for loop to the plotting part
single_star = False          # If only want to test one star set to True and type the path for a single star
save_centroid_disp = True   # To modify go to lines 306 and 307
# NOTE: for the names of stars, single numbers before the .fits require 8 spaces after quad_star, 
#       while 2 numbers require 7 spaces.  
single_star_path = paths_list[path_number]+'/postageout_star_     105 quad_       3 quad_star        5.fits'
display_master_img = False  # Want to see the combined ramped images for every star?
vlim = (0.001,10)             # sensitivity limits of image, i.e. (0.001, 0.1) 
debug = False               # see all debug messages (i.e. values of all calculations)
centroid_in_full_detector = False   # Want results in full detector coordinates or 32x32? 
checkbox_size = 3
background_method = None#'frac'  # Select either 'fractional', 'fixed', or None   
xwidth_list = [3, 5, 7]     # Number of rows of the centroid region
ywidth_list = [3, 5, 7]     # Number of columns of the centroid region
max_iter = 50
threshold = 1e-5
detector = 491
determine_moments = False   # Want to determine 2nd and 3rd moments?
# make plot of magnitude (in x) versus radial offset distance (in y) for Scene2 
show_plot = False
save_plot = True   # legend can be moved in line 789
plot_type = '.jpg'


###########################################################################################################

show_disp = False
if single_star:
    save_text_file = False
    show_disp = True
    #display_master_img = True   # Want to see the combined ramped images for every star?
    #debug = True   # see all debug messages (i.e. values of all calculations)
    #vlim = (1,10)#(0.001, .01)
    #       146

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
line0a = "{:<5} {:<15} {:<16} {:>23} {:>32} {:>33} {:>26} {:>15} {:>35} {:>38} {:>43}".format("Star", "Background", 
                                                                  "Centroid width: 3", "5", "7", 
                                                                  "TruePositions", "LoLeftCoords",
                                                                  "Factor",
                                                                  "Difference with: checkbox 3", "checkbox 5", "checkbox 7")
line0b = "{:>25} {:>12} {:>16} {:>14} {:>16} {:>14} {:>16} {:>14} {:>12} {:>10} {:>26} {:>16} {:>26} {:>16} {:>28} {:>16}".format(
                                                                       "x", "y", "x", "y", "x", "y", 
                                                                       "TrueX", "TrueY", "LoLeftX", "LoLeftY",
                                                                       "x", "y", "x", "y", "x", "y")
lines4screenandfile = [line0, line0a, line0b]

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

# Read fits table with benchmark data
if "scene1" in case:
    path2listfile = "PFforMaria/Scene_1_AB23"
    list_file = "simuTA20150528-F140X-S50-K-AB23.list"
    if 'shifted' in case: 
        list_file = "simuTA20150528-F140X-S50-K-AB23-shifted.list"
if "scene2" in case:
    # Read the text file just written to get the offsets from the "real" positions of the fake stars
    path2listfile = "PFforMaria/Scene_2_AB1823"
    list_file = "simuTA20150528-F140X-S50-K-AB18to23.list"
    if 'shifted' in case: 
        list_file = "simuTA20150528-F140X-S50-K-AB18to23-shifted.list"
lf = os.path.join(path2listfile,list_file)
star_number, xpos, ypos, factor, mag, bg_method = tf.read_listfile(lf, detector, background_method)

# Start the loop in the given directory
if just_read_text_file != True:
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
                factor_i = factor[idx_star]
                
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
                        cb_centroid_list = tf.run_recursive_centroids(psf, bg_frac, xwidth_list, ywidth_list, 
                                                               checkbox_size, max_iter, threshold, 
                                                               determine_moments, debug, display_master_img, vlim=vlim)
                        # Transform to full detector coordinates in order to compare with real centers
                        ESA_center = [0,0]
                        corr_true_center_centroid, corr_cb_centroid_list, loleftcoords, differences_true_TA = tf.transform2fulldetector(detector, 
                                                                                                      centroid_in_full_detector,
                                                                                                      cb_centroid_list, ESA_center, 
                                                                                                      true_center, perform_avgcorr=False)
                        # Write output into text file
                        bg = bg_frac
                        data2write = [save_text_file, output_file, st, bg, corr_cb_centroid_list, corr_true_center_centroid, loleftcoords, factor_i, differences_true_TA]
                        tf.write2file(data2write, lines4screenandfile)
                        #raw_input()
                else:
                    master_img_bgcorr = jtl.bg_correction(master_img, bg_method=background_method, 
                                                          bg_value=bg_value, bg_frac=bg_frac)
                    # Obtain the combined FITS image that combines all frames into one image AND
                    # check if all image is zeros, take the image that still has a max value
                    psf = tu.readimage(master_img_bgcorr, debug=debug)                
                    cb_centroid_list = tf.run_recursive_centroids(psf, bg_frac, xwidth_list, ywidth_list, 
                                                               checkbox_size, max_iter, threshold, 
                                                               determine_moments, debug, display_master_img, vlim=vlim)
                    # Transform to full detector coordinates in order to compare with real centers
                    ESA_center = [0,0]
                    corr_true_center_centroid, corr_cb_centroid_list, loleftcoords, differences_true_TA = tf.transform2fulldetector(detector, 
                                                                                                  centroid_in_full_detector,
                                                                                                  cb_centroid_list, ESA_center, 
                                                                                                  true_center, perform_avgcorr=False)
                    # Write output into text file
                    bg = background
                    data2write = [save_text_file, output_file, st, bg, corr_cb_centroid_list, corr_true_center_centroid, loleftcoords, factor_i, differences_true_TA]
                    tf.write2file(data2write, lines4screenandfile) 
                
                tf.display_centroids(st, case, psf, corr_true_center_centroid, corr_cb_centroid_list, show_disp, 
                                     vlim, savefile=save_centroid_disp)  
                if single_star:
                    tf.display_centroids(st, case, psf, corr_true_center_centroid, corr_cb_centroid_list, show_disp, 
                     vlim, savefile=save_centroid_disp)  
                    print ("Recursive test script finished. \n")
                    exit()

### Obtain standard deviation from true star positions
# Read the text file just written to get the offsets from the "real" positions of the fake stars
offsets = np.loadtxt(output_file, skiprows=3, usecols=(13,14,15,16,17,18), unpack=True)
sig3, mean3 = tf.find_std(offsets[1])
sig5, mean5 = tf.find_std(offsets[3])
sig7, mean7 = tf.find_std(offsets[5])
if 'frac' not in bg_method:
    fig1 = plt.figure(1, figsize=(12, 10))
    ax1 = fig1.add_subplot(111)
    plt.title(case+'_BG'+bg_method)
    plt.xlabel('Radial offset in X')
    plt.ylabel('Radial offset in Y')
    plt.plot(offsets[0], offsets[1], 'bo', ms=8, alpha=0.7, label='Checkbox=3')
    plt.plot(offsets[2], offsets[1], 'go', ms=8, alpha=0.7, label='Checkbox=5')
    plt.plot(offsets[4], offsets[1], 'ro', ms=8, alpha=0.7, label='Checkbox=7')
    xmin, xmax = ax1.get_xlim()
    plt.hlines(0.0, xmin, xmax, colors='k', linestyles='dashed')
    ymin, ymax = ax1.get_ylim()
    plt.vlines(0.0, ymin, ymax, colors='k', linestyles='dashed')
    plt.legend(loc='lower right')
    textinfig3 = r'$\sigma3$ = %0.2f    $\mu3$ = %0.2f' % (sig3, mean3)
    textinfig5 = r'$\sigma5$ = %0.2f    $\mu5$ = %0.2f' % (sig5, mean5)
    textinfig7 = r'$\sigma7$ = %0.2f    $\mu7$ = %0.2f' % (sig7, mean7)
    ax1.annotate(textinfig3, xy=(0.15, 0.055), xycoords='axes fraction' )
    ax1.annotate(textinfig5, xy=(0.15, 0.03), xycoords='axes fraction' )
    ax1.annotate(textinfig7, xy=(0.15, 0.005), xycoords='axes fraction' )
    if save_plot:
        if background_method is None:
            bg = 'None_'
        else:
            bg = 'fix_'
        destination = os.path.abspath("PFforMaria/plots/XoffsetVsYoffset_"+bg+case+plot_type)
        fig1.savefig(destination)
        print ("\n Plot saved: ", destination)
    if show_plot:
        plt.show()
    else:
        plt.close('all')
else:
    frac00, frac01, frac02, frac03, frac04, frac05, frac06, frac07, frac08, frac09, frac10 = tf.get_fracdata(offsets)
    frac_data  = [frac00, frac01, frac02, frac03, frac04, frac05, frac06, frac07, frac08, frac09, frac10]
    sig3, mean3, sig5, mean5, sig7, mean7 = tf.get_frac_stdevs(frac_data)
    frac_bgs = ['0.0', '0.1', '0.2', '0.3', '0.4', '0.5' ,'0.6' ,'0.7', '0.8', '0.9', '1.0']
    print ('\n{:<4} {:<20} {:>10}'.format('FrBG', 'Standard_deviation', 'Mean_y-offset'))
    print ('Checkbox sizes:')
    print ('{:<4} {:<6} {:<6} {:<6} {:<6} {:<6} {:<6}'.format('', '3', '5', '7', '3', '5', '7'))
    for fbg, s3, s5, s7, m3, m5, m7 in zip(frac_bgs, sig3, sig5, sig7, mean3, mean5, mean7):
        print('{:<4} {:<6.2f} {:<6.2f} {:<6.2f} {:<6.2f} {:<6.2f} {:<6.2f}'.format(fbg, s3, s5, s7, m3, m5, m7))
    fig2 = plt.figure(1, figsize=(12, 10))
    fig2.subplots_adjust(hspace=0.30)
    ax1 = fig2.add_subplot(311)
    ax1.set_title(case+"_BGfrac")
    ax1.set_xlabel('Radial offset in X: Checkbox=3')
    ax1.set_ylabel('Radial offset in Y: Checkbox=3')
    ax1.plot(frac00[0], frac00[1], 'bo', ms=8, alpha=0.7, label='bg_frac=0.0')
    ax1.plot(frac01[0], frac01[1], 'ro', ms=8, alpha=0.7, label='bg_frac=0.1')
    ax1.plot(frac02[0], frac02[1], 'mo', ms=8, alpha=0.7, label='bg_frac=0.2')
    ax1.plot(frac03[0], frac03[1], 'go', ms=5, alpha=0.7, label='bg_frac=0.3')
    ax1.plot(frac04[0], frac04[1], 'ko', ms=8, alpha=0.7, label='bg_frac=0.4')
    ax1.plot(frac05[0], frac05[1], 'yo', ms=8, alpha=0.7, label='bg_frac=0.5')
    ax1.plot(frac06[0], frac06[1], 'co', ms=8, alpha=0.7, label='bg_frac=0.6')
    ax1.plot(frac07[0], frac07[1], 'b+', ms=10, alpha=0.7, label='bg_frac=0.7')
    ax1.plot(frac08[0], frac08[1], 'r+', ms=8, alpha=0.7, label='bg_frac=0.8')
    ax1.plot(frac09[0], frac09[1], 'm+', ms=5, alpha=0.7, label='bg_frac=0.9')
    ax1.plot(frac10[0], frac10[1], 'k+', ms=5, alpha=0.7, label='bg_frac=1.0')
    xmin, xmax = ax1.get_xlim()
    plt.hlines(0.0, xmin, xmax, colors='k', linestyles='dashed')
    ymin, ymax = ax1.get_ylim()
    plt.vlines(0.0, ymin, ymax, colors='k', linestyles='dashed')
    #textinfig = r'$\sigma$ = %0.2f    $\mu$ = %0.2f' % (sig3, mean3)
    #ax1.annotate(textinfig, xy=(0.75, 0.05), xycoords='axes fraction' )
    # Shrink current axis by 10%
    box = ax1.get_position()
    ax1.set_position([box.x0, box.y0, box.width * 0.9, box.height])
    ax1.legend(loc='center left', bbox_to_anchor=(1, 0.5))   # put legend out of the plot box   
    ax2 = fig2.add_subplot(312)
    ax2.set_xlabel('Radial offset in X: Checkbox=5')
    ax2.set_ylabel('Radial offset in Y: Checkbox=5')
    ax2.plot(frac00[2], frac00[3], 'bo', ms=8, alpha=0.7, label='bg_frac=0.0')
    ax2.plot(frac01[2], frac01[3], 'ro', ms=8, alpha=0.7, label='bg_frac=0.1')
    ax2.plot(frac02[2], frac02[3], 'mo', ms=8, alpha=0.7, label='bg_frac=0.2')
    ax2.plot(frac03[2], frac03[3], 'go', ms=5, alpha=0.7, label='bg_frac=0.3')
    ax2.plot(frac04[2], frac04[3], 'ko', ms=8, alpha=0.7, label='bg_frac=0.4')
    ax2.plot(frac05[2], frac05[3], 'yo', ms=8, alpha=0.7, label='bg_frac=0.5')
    ax2.plot(frac06[2], frac06[3], 'co', ms=8, alpha=0.7, label='bg_frac=0.6')
    ax2.plot(frac07[2], frac07[3], 'b+', ms=10, alpha=0.7, label='bg_frac=0.7')
    ax2.plot(frac08[2], frac08[3], 'r+', ms=8, alpha=0.7, label='bg_frac=0.8')
    ax2.plot(frac09[2], frac09[3], 'm+', ms=5, alpha=0.7, label='bg_frac=0.9')
    ax2.plot(frac10[2], frac10[3], 'k+', ms=5, alpha=0.7, label='bg_frac=1.0')
    xmin, xmax = ax2.get_xlim()
    plt.hlines(0.0, xmin, xmax, colors='k', linestyles='dashed')
    ymin, ymax = ax2.get_ylim()
    plt.vlines(0.0, ymin, ymax, colors='k', linestyles='dashed')
    #textinfig = r'$\sigma$ = %0.2f    $\mu$ = %0.2f' % (sig5, mean5)
    #ax2.annotate(textinfig, xy=(0.75, 0.05), xycoords='axes fraction' )
    textinfig = r'BG      $\sigma$3     $\sigma$5     $\sigma$7'
    ax2.annotate(textinfig, xy=(1.02, 0.90), xycoords='axes fraction' )
    sigx = 1.02
    sigy = 0.9
    for fbg, s3, s5, s7 in zip(frac_bgs, sig3, sig5, sig7):
        line = ('{:<7} {:<6.2f} {:<6.2f} {:<6.2f}'.format(fbg, s3, s5, s7))
        sigy = sigy - 0.08
        ax2.annotate(line, xy=(sigx, sigy), xycoords='axes fraction' )
    # Shrink current axis by 10%
    box = ax2.get_position()
    ax2.set_position([box.x0, box.y0, box.width * 0.9, box.height])
    # put legend out of the plot box
    #ax2.legend(loc='center left', bbox_to_anchor=(1, 0.5))            
    ax3 = fig2.add_subplot(313)
    ax3.set_xlabel('Radial offset in X: Checkbox=7')
    ax3.set_ylabel('Radial offset in Y: Checkbox=7')
    ax3.plot(frac00[4], frac00[5], 'bo', ms=8, alpha=0.7, label='bg_frac=0.0')
    ax3.plot(frac01[4], frac01[5], 'ro', ms=8, alpha=0.7, label='bg_frac=0.1')
    ax3.plot(frac02[4], frac02[5], 'mo', ms=8, alpha=0.7, label='bg_frac=0.2')
    ax3.plot(frac03[4], frac03[5], 'go', ms=5, alpha=0.7, label='bg_frac=0.3')
    ax3.plot(frac04[4], frac04[5], 'ko', ms=8, alpha=0.7, label='bg_frac=0.4')
    ax3.plot(frac05[4], frac05[5], 'yo', ms=8, alpha=0.7, label='bg_frac=0.5')
    ax3.plot(frac06[4], frac06[5], 'co', ms=8, alpha=0.7, label='bg_frac=0.6')
    ax3.plot(frac07[4], frac07[5], 'b+', ms=10, alpha=0.7, label='bg_frac=0.7')
    ax3.plot(frac08[4], frac08[5], 'r+', ms=8, alpha=0.7, label='bg_frac=0.8')
    ax3.plot(frac09[4], frac09[5], 'm+', ms=5, alpha=0.7, label='bg_frac=0.9')
    ax3.plot(frac10[4], frac10[5], 'k+', ms=5, alpha=0.7, label='bg_frac=1.0')
    xmin, xmax = ax3.get_xlim()
    plt.hlines(0.0, xmin, xmax, colors='k', linestyles='dashed')
    ymin, ymax = ax3.get_ylim()
    plt.vlines(0.0, ymin, ymax, colors='k', linestyles='dashed')
    #textinfig = r'$\sigma$ = %0.2f    $\mu$ = %0.2f' % (sig7, mean7)
    #ax3.annotate(textinfig, xy=(0.75, 0.05), xycoords='axes fraction' )
    # Shrink current axis by 10%
    textinfig = r'BG      $\mu$3     $\mu$5     $\mu$7'
    ax3.annotate(textinfig, xy=(1.02, 0.90), xycoords='axes fraction' )
    sigx = 1.02
    sigy = 0.9
    for fbg, m3, m5, m7 in zip(frac_bgs, mean3, mean5, mean7):
        line = ('{:<7} {:<6.2f} {:<6.2f} {:<6.2f}'.format(fbg, m3, m5, m7))
        sigy = sigy - 0.08
        ax3.annotate(line, xy=(sigx, sigy), xycoords='axes fraction' )
    box = ax3.get_position()
    ax3.set_position([box.x0, box.y0, box.width * 0.9, box.height])
    if save_plot:
        destination = os.path.abspath("PFforMaria/plots/XoffsetVsYoffset_frac_"+case+plot_type)
        fig2.savefig(destination)
        print ("\n Plot saved: ", destination)
    if show_plot:
        plt.show()        
    plt.close('all')

# Make the plot of magnitude (in x) versus radial offset distance (in y) for Scene2
if "scene2" in case:
    # Read the text file just written to get the offsets from the "real" positions of the fake stars
    path2listfile = "PFforMaria/Scene_2_AB1823"
    list_file = "simuTA20150528-F140X-S50-K-AB18to23.list"
    if 'shifted' in case: 
        list_file = "simuTA20150528-F140X-S50-K-AB18to23-shifted.list"
    lf = os.path.join(path2listfile,list_file)
    star_number, xpos, ypos, factor, mag, bg_method = tf.read_listfile(lf, detector, background_method)
    print ('For Checkbox=3: ')
    sig3, mean3 = tf.find_std(offsets[1])
    print ('For Checkbox=5: ')
    sig5, mean5 = tf.find_std(offsets[3])
    print ('For Checkbox=7: ')
    sig7, mean7 = tf.find_std(offsets[5])
    if 'frac' not in bg_method:
        fig3 = plt.figure(1, figsize=(12, 10))
        ax1 = fig3.add_subplot(111)
        plt.title(case+'_BG'+bg_method)
        plt.xlabel('Magnitude')
        plt.ylabel('Radial offset in Y')
        plt.plot(mag, offsets[1], 'bo', ms=8, alpha=0.7, label='Checkbox=3')
        plt.plot(mag, offsets[3], 'go', ms=8, alpha=0.7, label='Checkbox=5')
        plt.plot(mag, offsets[5], 'ro', ms=8, alpha=0.7, label='Checkbox=7')
        #plt.legend(loc='lower left')
        plt.legend(loc='upper right')
        textinfig3 = r'$\sigma3$ = %0.2f    $\mu3$ = %0.2f' % (sig3, mean3)
        textinfig5 = r'$\sigma5$ = %0.2f    $\mu5$ = %0.2f' % (sig5, mean5)
        textinfig7 = r'$\sigma7$ = %0.2f    $\mu7$ = %0.2f' % (sig7, mean7)
        ax1.annotate(textinfig3, xy=(0.75, 0.055), xycoords='axes fraction' )
        ax1.annotate(textinfig5, xy=(0.75, 0.03), xycoords='axes fraction' )
        ax1.annotate(textinfig7, xy=(0.75, 0.005), xycoords='axes fraction' )
        xmin, xmax = ax1.get_xlim()
        plt.hlines(0.0, xmin, xmax, colors='k', linestyles='dashed')
        if save_plot:
            if background_method is None:
                bg = 'None_'
            else:
                bg = 'fix_'
            destination = os.path.abspath("PFforMaria/plots/MagVsYoffset_"+bg+case+plot_type)
            fig3.savefig(destination)
            print ("\n Plot saved: ", destination)
        if show_plot:
            plt.show()
        else:
            plt.close('all')
    else:
        #print ("Reading text file: ",  output_file)
        frac00, frac01, frac02, frac03, frac04, frac05, frac06, frac07, frac08, frac09, frac10 = tf.get_fracdata(offsets)
        frac_data  = [frac00, frac01, frac02, frac03, frac04, frac05, frac06, frac07, frac08, frac09, frac10]
        sig3, mean3, sig5, mean5, sig7, mean7 = tf.get_frac_stdevs(frac_data)
        frac_bgs = ['0.0', '0.1', '0.2', '0.3', '0.4', '0.5' ,'0.6' ,'0.7', '0.8', '0.9', '1.0']
        print ('\n{:<4} {:<20} {:>10}'.format('FrBG', 'Standard_deviation', 'Mean_y-offset'))
        print ('Checkbox sizes:')
        print ('{:<4} {:<6} {:<6} {:<6} {:<6} {:<6} {:<6}'.format('', '3', '5', '7', '3', '5', '7'))
        for fbg, s3, s5, s7, m3, m5, m7 in zip(frac_bgs, sig3, sig5, sig7, mean3, mean5, mean7):
            print('{:<4} {:<6.2f} {:<6.2f} {:<6.2f} {:<6.2f} {:<6.2f} {:<6.2f}'.format(fbg, s3, s5, s7, m3, m5, m7))
        fig4 = plt.figure(1, figsize=(12, 10))
        fig4.subplots_adjust(hspace=0.10)
        ax1 = fig4.add_subplot(311)
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
        #textinfig = r'$\sigma$ = %0.2f    $\mu$ = %0.2f' % (sig3, mean3)
        #ax1.annotate(textinfig, xy=(0.75, 0.05), xycoords='axes fraction' )
        # Shrink current axis by 10%
        box = ax1.get_position()
        ax1.set_position([box.x0, box.y0, box.width * 0.9, box.height])
        ax1.legend(loc='center left', bbox_to_anchor=(1, 0.5))            
        ax2 = fig4.add_subplot(312)
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
        textinfig = r'BG      $\sigma$3     $\sigma$5     $\sigma$7'
        ax2.annotate(textinfig, xy=(1.02, 0.88), xycoords='axes fraction' )
        sigx = 1.02
        sigy = 0.9
        for fbg, s3, s5, s7 in zip(frac_bgs, sig3, sig5, sig7):
            line = ('{:<7} {:<6.2f} {:<6.2f} {:<6.2f}'.format(fbg, s3, s5, s7))
            sigy = sigy - 0.08
            ax2.annotate(line, xy=(sigx, sigy), xycoords='axes fraction' )
        # Shrink current axis by 10%
        box = ax2.get_position()
        ax2.set_position([box.x0, box.y0, box.width * 0.9, box.height])
        # put legend out of the plot box
        #ax2.legend(loc='center left', bbox_to_anchor=(1, 0.5))            
        ax3 = fig4.add_subplot(313)
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
        textinfig = r'BG      $\mu$3     $\mu$5     $\mu$7'
        ax3.annotate(textinfig, xy=(1.02, 0.90), xycoords='axes fraction' )
        sigx = 1.02
        sigy = 0.9
        for fbg, m3, m5, m7 in zip(frac_bgs, mean3, mean5, mean7):
            line = ('{:<7} {:<6.2f} {:<6.2f} {:<6.2f}'.format(fbg, m3, m5, m7))
            sigy = sigy - 0.08
            ax3.annotate(line, xy=(sigx, sigy), xycoords='axes fraction' )
        # Shrink current axis by 10%
        box = ax3.get_position()
        ax3.set_position([box.x0, box.y0, box.width * 0.9, box.height])
        if save_plot:
            destination = os.path.abspath("PFforMaria/plots/MagVsYoffset_frac_"+case+plot_type)
            fig4.savefig(destination)
            print ("\n Plot saved: ", destination)
        if show_plot:
            plt.show()
        else:
            plt.close('all')
        
if not single_star:
    print ("\n Centroids and differences were written into: \n  {}".format(output_file))
print ("\n Recursive test script finished. Took  %s  seconds to finish. \n" % ((time.time() - start_time)) )
            