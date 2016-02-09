from __future__ import print_function, division
from glob import glob
from astropy.io import fits
import numpy as np
import os
import time
import random

# other code
import TA_functions as TAf 

print("Modules correctly imported! \n")

"""
DESCRIPTION:
    This script runs the centroid algorithm fully from beginning to end for the given initial
    conditions (which are set in the first part of the script). Once the centroids have been
    located, the script will choose X random stars from either detector and run the selected
    transformations for the specific test and case for the set of X stars:  
    
     TEST1 - Average positions P1 and P2, transform to V2-V3 space, and compare to average 
             reference positions (V2-V3 space)
     TEST2 - Transform individual positions P1 and P2 to V2-V3 space, average V2-V3 space 
             positions, and compare to average reference positions.
     TEST3 - Transform P1 and P2 individually to V2-V3 space and compare star by star and 
             position by position.

NOTES:
    * Scenes are defined by 100 stars in each detector:
        Scene or scenario 1 = All stars with magnitude 23.0
        Scene or scenario 2 = Magnitude range from 18.0 to 23.0
    
    * case depends on scene, noise, background value, and shutter velocity; results in 36 files per scene.
    
    * The tests are performed on ESA data in the following directories on central sotrage:
        - /grp/jwst/wit4/nirspec/PFforMaria/Scene_1_AB23
            (which contains 200 closer to "real" stars of magnitude 23: cosmic rays, hot 
            pixels, noise, and some shutters closed; and an ideal case)
        - /grp/jwst/wit4/nirspec/PFforMaria/Scene_2_AB1823
            (which contains 200 stars of different magnitudes and also has the "real" data
            as well as the ideal case)
        ** There are 2 sub-folders in each scene: rapid and slow. This is the shutter speed. Both
        sub-cases will be tested.
        ** The simulated data is described in detail in the NIRSpec Technical Note NTN-2015-013, which
        is in /grp/jwst/wit4/nirspec/PFforMaria/Documentation.



OUTPUT:
    - display images with true and calculated centroid
    - text file for the test ran with standard deviations and means for centroid window sizes 3, 5, and 7,
      sigma-clipped standard deviations and means, iterative least squares standard deviations and means,
      and the list of stars, background value used, the differences (in arcsecs or degrees) with respect
      to true or benchmark sky positions, and the centroid window size that has the minimum difference with
      respect to the true value.
"""


#######################################################################################################################


# INITIAL CONDITIONS

output_full_detector = True        # Give resulting coordinates in terms of full detector: True or False
save_text_file = False             # Want to save the text file of comparison? True or False
save_centroid_disp = False         # Save the display with measured and true positions?
keep_bad_stars = True              # Keep the bad stars in the sample? True or False
stars_in_sample = 20               # Number of stars in sample
scene = 2                          # Integer or string, scene=1 is constant Mag 23, scene=2 is stars with Mag 18-23
background_method = "frac"         # Select either 'fractional', 'fixed', or None   
background2use = 0.3               # Background to use for analysis: None or float
shutters = "rapid"                 # Shutter velocity, string: "rapid" or "slow"
noise = "real"                     # Noise level, string: "nonoise" or "real"
filter_input = "F140X"             # Filter, string: for now only test case is "F140X"
test2perform = "T3"                # Test to perform, string: "T1", "T2", "T3" for test 1, 2, and 3, respectively
Nsigma = 3                         # N-sigma rejection of bad stars: integer or float
max_iters_Nsig = 10                # Max number of iterations for N-sigma function: integer


# SECONDARY PARAMETERS THAT CAN BE ADJUSTED

checkbox_size = 3                  # Real checkbox size
xwidth_list = [3, 5, 7]            # Number of rows of the centroid region
ywidth_list = [3, 5, 7]            # Number of columns of the centroid region
vlim = (1, 100)                    # Sensitivity limits of image, i.e. (0.001, 0.1) 
threshold = 0.3                    # Convergence threshold of accepted difference between checkbox centroid and coarse location
max_iter = 10                      # Maximum number of iterations for finding coarse location
debug = False                      # See all debug messages (i.e. values of all calculations)
diffs_in_arcsecs = True            # Print the differences in arcsecs? True or False (=degrees) 
determine_moments = False          # Want to determine 2nd and 3rd moments?
display_master_img = False         # Want to see the combined ramped images for every star?
show_centroids = False             # Print measured centroid on screen: True or False
show_disp = False                  # Show display of resulting positions? (will show 2 figs, same but different contrast)
Pier_corr = True                   # Include Pier's corrections to measured positions
tilt = False                       # Tilt angle: True or False
backgnd_subtraction_method = 1     # 1    = Do background subtraction on final image (after subtracting 3-2 and 2-1), 
#                                           before converting negative values into zeros
#                                    2    = Do background subtraction on 3-2 and 2-1 individually
#                                    None = Do not subtract background

random_sample = True               # choose a random sample of stars from either detector: True or False
# control samples to be used when random is set to False
stars_sample = [7, 24, 51, 56, 66, 68, 71, 72, 74, 91, 106, 109, 120, 125, 127, 128, 138, 154, 187, 188]
# OLNY detector 491
#stars_sample = [101, 105, 108, 109, 111, 113, 114, 133, 136, 147, 150, 157, 158, 161, 181, 184, 185, 186, 194, 199]
#stars_sample = [101, 104, 105, 112, 117, 118, 133, 135, 136, 140, 145, 151, 152, 157, 159, 161, 174, 178, 184, 200]   
# ONLY detector 492
#stars_sample = [8, 11, 19, 24, 30, 37, 39, 41, 48, 51, 55, 65, 73, 85, 87, 88, 90, 91, 93, 98]
#stars_sample = [2, 4, 8, 10, 11, 22, 25, 28, 33, 37, 54, 64, 68, 76, 80, 89, 96, 97, 99, 100]
# all stars of one detector or both
#stars_sample = [s+101 for s in range(100)]
# Known bad stars in X and Y: 103, 105, 106, 112, 134, 152, 156, 170, 188
#6, 23, 50, 55, 65, 67, 70, 71, 73, 90, 105, 108, 119, 124, 126, 127, 137, 153, 186, 187

#######################################################################################################################


detectors = [491, 492]

# Stars of detector 491 and 492
stars_detectors = range(1, 201)

if random_sample:
    # select stars_in_sample stars from 1 to 200
    stars_sample = []
    for i in range(stars_in_sample):
        random_star = random.choice(stars_detectors)
        stars_sample.append(random_star)
    # make sure that there are no repetitions
    stars_sample = list(set(stars_sample))
    while len(stars_sample) != stars_in_sample:
        random_star = random.choice(stars_detectors)
        stars_sample.append(random_star)
        stars_sample = list(set(stars_sample))
        # remove the bad stars
        if not keep_bad_stars:
            TAf.remove_bad_stars(stars_sample)
# remove the bad stars
if not keep_bad_stars:
    TAf.remove_bad_stars(stars_sample)
    print (" * Sample has %i stars left." % len(stars_sample))

# order the star list
stars_sample.sort(key=lambda xx: xx)
print ("stars_sample =", stars_sample)
#raw_input("\nPress enter to continue...")

# start the timer to compute the whole running time
start_time = time.time()

# Background cases variable setting
bg_frac, bg_value = None, None   # for the None case
bg_choice = "_bgNone"
if background_method is not None:
    if "frac" in background_method:
        bg_frac = background2use
        bg_choice = "_bgFrac"
    elif "fix" in background_method:
        bg_value = background2use
        bg_choice = "_bgFixed"
else:
    background2use = 0.0

# Paths to Scenes 1 and 2 local directories: /Users/pena/Documents/AptanaStudio3/NIRSpec/TargetAcquisition/
path4starfiles = "../PFforMaria/"

# Define the paths for results
path4results = "../resultsXrandomstars/"

# Set the case to study according to the selected scene
scene2study = "Scene"+str(scene)+"_"
case = "Scene"+str(scene)+"_"+str(shutters)+"_"+noise

# get the benchmark data according to Scene selected
benchmark_data, magnitudes = TAf.read_star_param_files(scene2study)
bench_P1, bench_P2 = benchmark_data
allbench_starP1, allbench_xP1, allbench_yP1, allbench_V2P1, allbench_V3P1, allbench_xLP1, allbench_yLP1 = bench_P1
allbench_starP2, allbench_xP2, allbench_yP2, allbench_V2P2, allbench_V3P2, allbench_xLP2, allbench_yLP2 = bench_P2
allbench_stars = allbench_starP1.tolist()

# get the index for the sample stars
star_idx_list = []
for st in stars_sample:
    st_idx = allbench_stars.index(st)
    star_idx_list.append(st_idx)

# get the benchmark for star sample
bench_starP1, bench_xP1, bench_yP1, bench_V2P1, bench_V3P1, bench_xLP1, bench_yLP1 = np.array([]),np.array([]),np.array([]),np.array([]),np.array([]),np.array([]),np.array([])
bench_starP2, bench_xP2, bench_yP2, bench_V2P2, bench_V3P2, bench_xLP2, bench_yLP2 = np.array([]),np.array([]),np.array([]),np.array([]),np.array([]),np.array([]),np.array([])
for i in star_idx_list:
    bench_starP1 = np.append(bench_starP1, allbench_starP1[i])
    bench_xP1 = np.append(bench_xP1, allbench_xP1[i])
    bench_yP1 = np.append(bench_yP1, allbench_yP1[i])
    bench_V2P1 = np.append(bench_V2P1, allbench_V2P1[i])
    bench_V3P1 = np.append(bench_V3P1, allbench_V3P1[i])
    bench_xLP1 = np.append(bench_xLP1, allbench_xLP1[i])
    bench_yLP1 = np.append(bench_yLP1, allbench_yLP1[i])
    bench_starP2 = np.append(bench_starP2, allbench_starP2[i])
    bench_xP2 = np.append(bench_xP2, allbench_xP2[i])
    bench_yP2 = np.append(bench_yP2, allbench_yP2[i])
    bench_V2P2 = np.append(bench_V2P2, allbench_V2P2[i])
    bench_V3P2 = np.append(bench_V3P2, allbench_V3P2[i])
    bench_xLP2 = np.append(bench_xLP2, allbench_xLP2[i])
    bench_yLP2 = np.append(bench_yLP2, allbench_yLP2[i])
trueVsP1 = [bench_V2P1, bench_V3P1]
trueVsP2 = [bench_V2P2, bench_V3P2]
LoLeftCornersP1 = [bench_xLP1, bench_yLP1]
LoLeftCornersP2 = [bench_xLP2, bench_yLP2]


### Perform centroid algorithm for stars sample

# start the text file with the measured centroids
output_file_path = "../resultsXrandomstars/centroid_txt_files/"
line0 = "Centroid indexing starting at 1 !"
line0a = "{:<5} {:<15} {:<16} {:>23} {:>32} {:>40} {:>25} {:>15} {:>9}".format("Star", "Background",
                                                                  "Centroid width: 3", "5", "7",
                                                                  "TruePositions", "LoLeftCoords",
                                                                  "Magnitude",
                                                                  "MinDiff")
line0b = "{:>25} {:>12} {:>16} {:>14} {:>16} {:>20} {:>16} {:>17} {:>11} {:>11} {:>16} {:>2}".format(
                                                                       "x", "y", "x", "y", "x", "y",
                                                                       "TrueX", "TrueY", "LoLeftX", "LoLeftY",
                                                                       "x", "y")
lines4screenandfile = [line0, line0a, line0b]
# write the file
positions = ["_Position1", "_Position2"]
if save_text_file:
    for pos in positions:
        output_file = os.path.join(output_file_path, "centroids_Scene"+repr(scene)+bg_choice+pos+".txt")
        f = open(output_file, "w+")
        f.write(line0+"\n")
        f.close()

# get the star files to run the TA algorithm on
dir2test_list = TAf.get_raw_star_directory(path4starfiles, scene, shutters, noise)

# run centroid algorithm on each position and save them into a text file
x13, x15, x17 = np.array([]), np.array([]), np.array([])
y13, y15, y17 = np.array([]), np.array([]), np.array([])
x23, x25, x27 = np.array([]), np.array([]), np.array([])
y23, y25, y27 = np.array([]), np.array([]), np.array([])
min_diff_pixposX, min_diff_pixposY, loleftcoords_list, mag_list = [], [], [], []
true_centers=[]
for pos, dir2test in zip(positions, dir2test_list):
    dir_stars = glob(os.path.join(dir2test,"postageout_star_*.fits"))   # get all star fits files in that directory
    #print("does dir2test exist?", os.path.isdir(dir2test))
    for star in dir_stars:
        dir_star_number = int(os.path.basename(star).split()[1])
        # Test stars of detector of choice
        for st in stars_sample:
            if st == dir_star_number: #if str(st)+" quad_       " in star:
                print ("Will test stars in directory: \n     ", dir2test)
                print ("Star: ", os.path.basename(star))
                # Make sure the file actually exists
                star_exists = os.path.isfile(star)
                if not star_exists:
                    print ("The file: ", star, "\n    does NOT exist. Exiting the script.")
                    exit()

                # Obtain real star position and corresponding detector
                if st < 100:
                    detector = detectors[1]
                else:
                    detector = detectors[0]
                idx_star = stars_sample.index(st)
                mag_i = magnitudes[idx_star]
                true_center = [bench_xP1[idx_star], bench_yP1[idx_star]]
                if pos == "_Position2":
                    true_center = [bench_xP2[idx_star], bench_yP2[idx_star]]

                # Read FITS image
                print ("Running centroid algorithm... ")
                #hdr = fits.getheader(star, 0)
                #print("** HEADER:", hdr)
                master_img = fits.getdata(star, 0)
                print ('Master image shape: ', np.shape(master_img))
                # Obtain the combined FITS image that combines all frames into one image
                # background subtraction is done here
                psf = TAf.readimage(master_img, backgnd_subtraction_method, bg_method=background_method,
                                    bg_value=bg_value, bg_frac=bg_frac, debug=debug)
                cb_centroid_list_in32x32pix = TAf.run_recursive_centroids(psf, bg_frac, xwidth_list, ywidth_list,
                                                           checkbox_size, max_iter, threshold,
                                                           determine_moments, debug)
                cb_centroid_list, loleftcoords, true_center32x32, differences_true_TA = TAf.centroid2fulldetector(cb_centroid_list_in32x32pix,
                                                                                                    true_center)
                if not output_full_detector:
                    cb_centroid_list = cb_centroid_list_in32x32pix
                    true_center = true_center32x32
                # Correct true centers for average value given by Pier
                if Pier_corr:
                    corr_cb_centroid_list = TAf.do_Piers_correction(detector, cb_centroid_list)
                else:
                    corr_cb_centroid_list = cb_centroid_list
                if show_centroids:
                    print ('***** Measured centroids for centroid window sizes 3, 5, and 7, respectively:')
                    if output_full_detector:
                        print ('      cb_centroid_list = ', corr_cb_centroid_list)
                    else:
                        print ('      cb_centroid_list = ', cb_centroid_list_in32x32pix)
                    print ('           True center = ', true_center)
                # Show the display with the measured and true positions
                fig_name = os.path.join("../resultsXrandomstars", "centroid_displays/Star"+repr(st)+"_Scene"+repr(scene)+bg_choice+pos+".jpg")
                # Display the combined FITS image that combines all frames into one image
                m_img = display_master_img
                if display_master_img:
                    m_img = TAf.readimage(master_img, backgnd_subtraction_method=None, bg_method=None,
                                      bg_value=None, bg_frac=None, debug=False)
                TAf.display_centroids(detector, st, case, psf, true_center32x32, cb_centroid_list_in32x32pix,
                                     show_disp, vlim, savefile=save_centroid_disp, fig_name=fig_name, display_master_img=m_img)
                # Find the best centroid window size = minimum difference with true values
                min_diff, _ = TAf.get_mindiff(differences_true_TA[0][0], differences_true_TA[0][1], differences_true_TA[0][2])
                # Save output
                true_centers.append(true_center)
                loleftcoords_list.append(loleftcoords)
                mag_list.append(mag_i)
                min_diff_pixposX.append(min_diff[0])
                min_diff_pixposY.append(min_diff[1])
                if pos == "_Position1":
                    x13 = np.append(x13, corr_cb_centroid_list[0][0])
                    x15 = np.append(x15, corr_cb_centroid_list[1][0])
                    x17 = np.append(x17, corr_cb_centroid_list[2][0])
                    y13 = np.append(y13, corr_cb_centroid_list[0][1])
                    y15 = np.append(y15, corr_cb_centroid_list[1][1])
                    y17 = np.append(y17, corr_cb_centroid_list[2][1])
                if pos == "_Position2":
                    x23 = np.append(x23, corr_cb_centroid_list[0][0])
                    x25 = np.append(x25, corr_cb_centroid_list[1][0])
                    x27 = np.append(x27, corr_cb_centroid_list[2][0])
                    y23 = np.append(y23, corr_cb_centroid_list[0][1])
                    y25 = np.append(y25, corr_cb_centroid_list[1][1])
                    y27 = np.append(y27, corr_cb_centroid_list[2][1])
    # Write output into text file
    position = "_Position1"
    x_pixpos = [x13, x15, x17]
    y_pixpos = [y13, y15, y17]
    if pos == "_Position2":
        x2_pixpos = [x23, x25, x27]
        y2_pixpos = [y23, y25, y27]
        position = "_Position2"
    output_file = os.path.join(output_file_path, "centroids_Scene"+repr(scene)+bg_choice+position+".txt")
    data2write = [x_pixpos, y_pixpos, true_centers, loleftcoords_list, mag_list, min_diff_pixposX, min_diff_pixposY]
    TAf.writePixPos(save_text_file, show_centroids, output_file, lines4screenandfile, stars_sample, background2use, data2write)

# transform into sky coordinates
case2study = [scene, shutters, noise, bg_choice]
case = "Scene"+str(scene)+"_"+shutters+"_"+noise+bg_choice

if debug:
    print ("Check that read BENCHMARK values correspond to expected for case: ", case)
    print ("Star, xP1, yP1, V2P1, V3P1, xLP1, yLP1")
    print (bench_starP1[0], bench_xP1[0], bench_yP1[0], bench_V2P1[0], bench_V3P1[0], bench_xLP1[0], bench_yLP1[0])
    print ("Star, xP2, yP2, V2P2, V3P2, xLP2, yLP2")
    print (bench_starP2[0], bench_xP2[0], bench_yP2[0], bench_V2P2[0], bench_V3P2[0], bench_xLP2[0], bench_yLP2[0])
    print ("Check that read MEASURED values correspond to expected for the same case: ", case)
    print ("   -> reading measured info from: ", case)
    print ("Star, BG, x13, y13, x15, y15, x17, y17, LoLeftP1 (x, y), TrueP1 (x, y)")
    print (stars_sample[0], bg_choice, x13[0], y13[0], x15[0], y15[0], x17[0], y17[0], bench_xLP1[0], bench_yLP1[0], bench_xP1[0], bench_yP1[0])
    print ("Star, BG, x23, y23, x25, y25, x27, y27, LoLeftP2 (x, y), TrueP2 (x, y)")
    print (stars_sample[0], bg_choice, x23[0], y23[0], x25[0], y25[0], x27[0], y27[0], bench_xLP2[0], bench_yLP2[0], bench_xP2[0], bench_yP2[0])
    raw_input(" * press enter to continue... \n")

# show positions on screen
line0 = "\n Centroid indexing starting at 1 !"
line0a = "{:<5} {:<15} {:<16} {:>23} {:>30} {:>44} {:>17} {:>15}".format("Star", "Background",
                                                                  "Centroid windows: 3", "5", "7",
                                                                  "TruePositions", "LoLeftCoords",
                                                                  "Mag")
line0b = "{:>25} {:>12} {:>16} {:>14} {:>16} {:>14} {:>16} {:>18} {:>12} {:>10}".format(
                                                                       "x", "y", "x", "y", "x", "y",
                                                                       "TrueX", "TrueY", "LoLeftX", "LoLeftY")
print ("Analyzing case: ", case)
print (line0)
print (line0a)
print (line0b)
for i, st in enumerate(stars_sample):
    line1 = "{:<5} {:<10} {:<14} {:<16} {:<14} {:<16} {:<14} {:<16} {:<14} {:<16} {:<8} {:<12} {:<10.2f}".format(
                                                                int(st), background2use,
                                                                x13[i], y13[i], x15[i], y15[i], x17[i], y17[i],
                                                                bench_xP1[i]-bench_xLP1[i], bench_yP1[i]-bench_yLP1[i],
                                                                bench_xLP1[i], bench_yLP1[i],
                                                                magnitudes[i])
    print (line1)

# compact results for functions
P1P2data = [x13,y13, x23,y23, x15,y15, x25,y25, x17,y17, x27,y27]

# Now run the tests
transf_direction = "forward"
# TEST 1: (a) Avg P1 and P2, (b) transform to V2-V3, (c) compare to avg reference positions (V2-V3 space)
if test2perform == "T1":
    resultsTEST1 = TAf.runTEST(test2perform, detectors, transf_direction, case, stars_sample, P1P2data, bench_starP1,
                               trueVsP1, trueVsP2, filter_input, tilt, diffs_in_arcsecs, debug)
    T1P1P2data, T1_transformations, T1_diffs, T1_benchVs_list = resultsTEST1
    T1_V2_3, T1_V3_3, T1_V2_5, T1_V3_5, T1_V2_7, T1_V3_7 = T1_transformations
    T1_diffV2_3, T1_diffV3_3, T1_diffV2_5, T1_diffV3_5, T1_diffV2_7, T1_diffV3_7 = T1_diffs
    T1bench_V2_list, T1bench_V3_list = T1_benchVs_list
    # Get the statistics
    results_stats = TAf.get_stats(case, T1_transformations, T1_diffs, T1_benchVs_list, Nsigma, max_iters_Nsig)
    # unfold results
    T1_st_devsAndMeans, T1_diff_counter, T1_bench_values, T1_sigmas_deltas, T1_sigma_reject, rejected_elementsLS, rejected_eleNsig = results_stats
    T1stdev_V2_3, T1mean_V2_3, T1stdev_V2_5, T1mean_V2_5, T1stdev_V2_7, T1mean_V2_7, T1stdev_V3_3, T1mean_V3_3, T1stdev_V3_5, T1mean_V3_5, T1stdev_V3_7, T1mean_V3_7 = T1_st_devsAndMeans
    T1_min_diff, T1_counter = T1_diff_counter
    T1bench_V2, T1bench_V3 = T1_bench_values
    T1LSdeltas_3, T1LSsigmas_3, T1LSlines2print_3, T1LSdeltas_5, T1LSsigmas_5, T1LSlines2print_5, T1LSdeltas_7, T1LSsigmas_7, T1LSlines2print_7 = T1_sigmas_deltas
    T1sigmaV2_3, T1meanV2_3, T1sigmaV3_3, T1meanV3_3, T1newV2_3, T1newV3_3, T1niter_3, T1lines2print_3, T1sigmaV2_5, T1meanV2_5, T1sigmaV3_5, T1meanV3_5, T1newV2_5, T1newV3_5, T1niter_5, T1lines2print_5, T1sigmaV2_7, T1meanV2_7, T1sigmaV3_7, T1meanV3_7, T1newV2_7, T1newV3_7, T1niter_7, T1lines2print_7 = T1_sigma_reject

# TEST 2: (a) Transform individual P1 and P2 to V2-V3, (b) avg V2-V3 space positions, (c) compare to avg reference positions
if test2perform == "T2":
    resultsTEST2 = TAf.runTEST(test2perform, detectors, transf_direction, case, stars_sample, P1P2data, bench_starP1,
                               trueVsP1, trueVsP2, filter_input, tilt, diffs_in_arcsecs, debug)
    T2P1P2data, T2_transformations, T2_diffs, T2_benchVs_list = resultsTEST2
    T2_V2_3, T2_V3_3, T2_V2_5, T2_V3_5, T2_V2_7, T2_V3_7 = T2_transformations
    T2_diffV2_3, T2_diffV3_3, T2_diffV2_5, T2_diffV3_5, T2_diffV2_7, T2_diffV3_7 = T2_diffs
    T2bench_V2_list, T2bench_V3_list = T2_benchVs_list
    # Get the statistics
    results_stats = TAf.get_stats(case, T2_transformations, T2_diffs, T2_benchVs_list, Nsigma, max_iters_Nsig)
    # unfold results
    T2_st_devsAndMeans, T2_diff_counter, T2_bench_values, T2_sigmas_deltas, T2_sigma_reject, rejected_elementsLS, rejected_eleNsig = results_stats
    T2stdev_V2_3, T2mean_V2_3, T2stdev_V2_5, T2mean_V2_5, T2stdev_V2_7, T2mean_V2_7, T2stdev_V3_3, T2mean_V3_3, T2stdev_V3_5, T2mean_V3_5, T2stdev_V3_7, T2mean_V3_7 = T2_st_devsAndMeans
    T2_min_diff, T2_counter = T2_diff_counter
    T2bench_V2, T2bench_V3 = T2_bench_values
    T2LSdeltas_3, T2LSsigmas_3, T2LSlines2print_3, T2LSdeltas_5, T2LSsigmas_5, T2LSlines2print_5, T2LSdeltas_7, T2LSsigmas_7, T2LSlines2print_7 = T2_sigmas_deltas
    T2sigmaV2_3, T2meanV2_3, T2sigmaV3_3, T2meanV3_3, T2newV2_3, T2newV3_3, T2niter_3, T2lines2print_3, T2sigmaV2_5, T2meanV2_5, T2sigmaV3_5, T2meanV3_5, T2newV2_5, T2newV3_5, T2niter_5, T2lines2print_5, T2sigmaV2_7, T2meanV2_7, T2sigmaV3_7, T2meanV3_7, T2newV2_7, T2newV3_7, T2niter_7, T2lines2print_7 = T2_sigma_reject

# TEST 3: (a) Transform P1 and P2 individually to V2-V3 (b) compare star by star and position by position
if test2perform == "T3":
    resultsTEST3 = TAf.runTEST(test2perform, detectors, transf_direction, case, stars_sample, P1P2data, bench_starP1,
                               trueVsP1, trueVsP2, filter_input, tilt, diffs_in_arcsecs, debug)
    T3P1P2data, T3_transformations, T3_diffs, T3_benchVs_list = resultsTEST3
    x13,y13, x23,y23, x15,y15, x25,y25, x17,y17, x27,y27 = T3P1P2data
    T_V2_3, T_V3_3, T_V2_5, T_V3_5, T_V2_7, T_V3_7 = T3_transformations
    T3_V2_13, T3_V2_23 = T_V2_3
    T3_V3_13, T3_V3_23 = T_V3_3
    T3_V2_15, T3_V2_25 = T_V2_5
    T3_V3_15, T3_V3_25 = T_V3_5
    T3_V2_17, T3_V2_27 = T_V2_7
    T3_V3_17, T3_V3_27 = T_V3_7
    T_diffV2_3, T_diffV3_3, T_diffV2_5, T_diffV3_5, T_diffV2_7, T_diffV3_7 = T3_diffs
    T3_diffV2_13, T3_diffV2_23 = T_diffV2_3
    T3_diffV3_13, T3_diffV3_23 = T_diffV3_3
    T3_diffV2_15, T3_diffV2_25 = T_diffV2_5
    T3_diffV3_15, T3_diffV3_25 = T_diffV3_5
    T3_diffV2_17, T3_diffV2_27 = T_diffV2_7
    T3_diffV3_17, T3_diffV3_27 = T_diffV3_7
    T3bench_V2_list, T3bench_V3_list = T3_benchVs_list
    T3bench_V2_listP1, T3bench_V2_listP2 = T3bench_V2_list
    T3bench_V3_listP1, T3bench_V3_listP2 = T3bench_V3_list
    # combine the arrays (positions 1 and 2)
    T3_V2_3, T3_V2_5, T3_V2_7 = np.array([]), np.array([]), np.array([])
    T3_V2_3 = TAf.combine2arrays(T3_V2_13, T3_V2_23, T3_V2_3)
    T3_V2_5 = TAf.combine2arrays(T3_V2_15, T3_V2_25, T3_V2_5)
    T3_V2_7 = TAf.combine2arrays(T3_V2_17, T3_V2_27, T3_V2_7)
    T3_V3_3, T3_V3_5, T3_V3_7 = np.array([]), np.array([]), np.array([])
    T3_V3_3 = TAf.combine2arrays(T3_V3_13, T3_V3_23, T3_V3_3)
    T3_V3_5 = TAf.combine2arrays(T3_V3_15, T3_V3_25, T3_V3_5)
    T3_V3_7 = TAf.combine2arrays(T3_V3_17, T3_V3_27, T3_V3_7)
    T3_diffV2_3, T3_diffV2_5, T3_diffV2_7 = np.array([]), np.array([]), np.array([])
    T3_diffV2_3 = TAf.combine2arrays(T3_diffV2_13, T3_diffV2_23, T3_diffV2_3)
    T3_diffV2_5 = TAf.combine2arrays(T3_diffV2_15, T3_diffV2_25, T3_diffV2_5)
    T3_diffV2_7 = TAf.combine2arrays(T3_diffV2_17, T3_diffV2_27, T3_diffV2_7)
    T3_diffV3_3, T3_diffV3_5, T3_diffV3_7 = np.array([]), np.array([]), np.array([])
    T3_diffV3_3 = TAf.combine2arrays(T3_diffV3_13, T3_diffV3_23, T3_diffV3_3)
    T3_diffV3_5 = TAf.combine2arrays(T3_diffV3_15, T3_diffV3_25, T3_diffV3_5)
    T3_diffV3_7 = TAf.combine2arrays(T3_diffV3_17, T3_diffV3_27, T3_diffV3_7)
    T3bench_V2_list, T3bench_V3_list = np.array([]), np.array([])
    T3bench_V2_list = TAf.combine2arrays(np.array(T3bench_V2_listP1), np.array(T3bench_V2_listP2), T3bench_V2_list)
    T3bench_V3_list = TAf.combine2arrays(np.array(T3bench_V3_listP1), np.array(T3bench_V3_listP2), T3bench_V3_list)
    T3bench_V2_list.tolist()
    T3bench_V3_list.tolist()
    # Get the statistics
    T3_transformations = [T3_V2_3, T3_V3_3, T3_V2_5, T3_V3_5, T3_V2_7, T3_V3_7]
    T3_diffs = [T3_diffV2_3, T3_diffV3_3, T3_diffV2_5, T3_diffV3_5, T3_diffV2_7, T3_diffV3_7]
    T3_benchVs_list = [T3bench_V2_list, T3bench_V3_list]
    results_stats = TAf.get_stats(case, T3_transformations, T3_diffs, T3_benchVs_list, Nsigma, max_iters_Nsig)
    # unfold results
    T3_st_devsAndMeans, T3_diff_counter, T3_bench_values, T3_sigmas_deltas, T3_sigma_reject, rejected_elementsLS, rejected_eleNsig = results_stats
    T3stdev_V2_3, T3mean_V2_3, T3stdev_V2_5, T3mean_V2_5, T3stdev_V2_7, T3mean_V2_7, T3stdev_V3_3, T3mean_V3_3, T3stdev_V3_5, T3mean_V3_5, T3stdev_V3_7, T3mean_V3_7 = T3_st_devsAndMeans
    T3_min_diff, T3_counter = T3_diff_counter
    T3bench_V2, T3bench_V3 = T3_bench_values
    T3LSdeltas_3, T3LSsigmas_3, T3LSlines2print_3, T3LSdeltas_5, T3LSsigmas_5, T3LSlines2print_5, T3LSdeltas_7, T3LSsigmas_7, T3LSlines2print_7 = T3_sigmas_deltas
    T3sigmaV2_3, T3meanV2_3, T3sigmaV3_3, T3meanV3_3, T3newV2_3, T3newV3_3, T3niter_3, T3lines2print_3, T3sigmaV2_5, T3meanV2_5, T3sigmaV3_5, T3meanV3_5, T3newV2_5, T3newV3_5, T3niter_5, T3lines2print_5, T3sigmaV2_7, T3meanV2_7, T3sigmaV3_7, T3meanV3_7, T3newV2_7, T3newV3_7, T3niter_7, T3lines2print_7 = T3_sigma_reject

# Print results to screen and save into a text file if told so
if test2perform == "T1":
    Tstdev_Vs = [T1stdev_V2_3, T1stdev_V3_3, T1stdev_V2_5, T1stdev_V3_5, T1stdev_V2_7, T1stdev_V3_7]
    Tmean_Vs = [T1mean_V2_3, T1mean_V3_3, T1mean_V2_5, T1mean_V3_5, T1mean_V2_7, T1mean_V3_7]
    T_diff_counter = [T1_min_diff, T1_counter]
    TLSlines2print = [T1LSlines2print_3, T1LSlines2print_5, T1LSlines2print_7]
    Tlines2print = [T1lines2print_3, T1lines2print_5, T1lines2print_7]
    Tbench_Vs_list = [T1bench_V2_list, T1bench_V3_list]
    T_Vs = [T1_V2_3, T1_V3_3, T1_V2_5, T1_V3_5, T1_V2_7, T1_V3_7]
    T_diffVs = [T1_diffV2_3, T1_diffV3_3, T1_diffV2_5, T1_diffV3_5, T1_diffV2_7, T1_diffV3_7]

if test2perform == "T2":
    Tstdev_Vs = [T2stdev_V2_3, T2stdev_V3_3, T2stdev_V2_5, T2stdev_V3_5, T2stdev_V2_7, T2stdev_V3_7]
    Tmean_Vs = [T2mean_V2_3, T2mean_V3_3, T2mean_V2_5, T2mean_V3_5, T2mean_V2_7, T2mean_V3_7]
    T_diff_counter = [T2_min_diff, T2_counter]
    TLSlines2print = [T2LSlines2print_3, T2LSlines2print_5, T2LSlines2print_7]
    Tlines2print = [T2lines2print_3, T2lines2print_5, T2lines2print_7]
    Tbench_Vs_list = [T2bench_V2_list, T2bench_V3_list]
    T_Vs = [T2_V2_3, T2_V3_3, T2_V2_5, T2_V3_5, T2_V2_7, T2_V3_7]
    T_diffVs = [T2_diffV2_3, T2_diffV3_3, T2_diffV2_5, T2_diffV3_5, T2_diffV2_7, T2_diffV3_7]

if test2perform == "T3":
    Tstdev_Vs = [T3stdev_V2_3, T3stdev_V3_3, T3stdev_V2_5, T3stdev_V3_5, T3stdev_V2_7, T3stdev_V3_7]
    Tmean_Vs = [T3mean_V2_3, T3mean_V3_3, T3mean_V2_5, T3mean_V3_5, T3mean_V2_7, T3mean_V3_7]
    T_diff_counter = [T3_min_diff, T3_counter]
    TLSlines2print = [T3LSlines2print_3, T3LSlines2print_5, T3LSlines2print_7]
    Tlines2print = [T3lines2print_3, T3lines2print_5, T3lines2print_7]
    Tbench_Vs_list = [T3bench_V2_list, T3bench_V3_list]
    T_Vs = [T3_V2_3, T3_V3_3, T3_V2_5, T3_V3_5, T3_V2_7, T3_V3_7]
    T_diffVs = [T3_diffV2_3, T3_diffV3_3, T3_diffV2_5, T3_diffV3_5, T3_diffV2_7, T3_diffV3_7]


TAf.printTESTresults(stars_sample, case, test2perform, diffs_in_arcsecs, Tstdev_Vs, Tmean_Vs, T_diff_counter,
              save_text_file, TLSlines2print, Tlines2print, Tbench_Vs_list, T_Vs, T_diffVs,
              rejected_elementsLS, rejected_eleNsig, background_method, background2use, path4results)


print ("\n Script 'testXrandom_stars.py' finished! Took  %s  seconds to finish. \n" % (time.time() - start_time))
