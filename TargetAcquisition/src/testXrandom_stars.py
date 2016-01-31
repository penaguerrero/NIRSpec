from __future__ import print_function, division
from glob import glob
from astropy.io import fits
import numpy as np
import os
import time
#import collections
import random
#import PIL.Image as Image
# other code
import coords_transform as ct
import TA_functions as tf 
import least_squares_iterate as lsi
# Tommy's code
import tautils as tu
import jwst_targloc as jtl

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
    - text file for the test ran with standard deviations and means for checkboxes 3, 5, and 7,
      sigma-clipped standard deviations and means, iterative least squares standard deviations 
      and means, and the list of stars, background value used, the differences (in arcsecs or
      degrees) with respect to true or benchmark sky positions, and the checkbox size that has
      the minimum difference with respect to the true value.  
"""


#######################################################################################################################


# INITIAL CONDITIONS

output_full_detector = True        # Give resulting coordinates in terms of full detector: True or False
save_text_file = False             # Want to save the text file of comparison? True or False
save_centroid_disp = False         # Save the display with measured and true positions?
keep_bad_stars = True              # Keep the bad stars in the sample? True or False
stars_in_sample = 20               # Number of stars in sample
scene = 2                          # Integer or string, scene=1 is constant Mag 23, scene=2 is stars with Mag 18-23
background_method = 'frac'         # Select either 'fractional', 'fixed', or None   
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
vlim = (1.0,30)                    # Sensitivity limits of image, i.e. (0.001, 0.1) 
threshold = 1e-5                   # Convergence threshold of accepted difference between checkbox centroid and coarse location
max_iter = 50                      # Maximum number of iterations for finding coarse location
debug = False                      # See all debug messages (i.e. values of all calculations)
diffs_in_arcsecs = True            # Print the differences in arcsecs? True or False (=degrees) 
determine_moments = False          # Want to determine 2nd and 3rd moments?
display_master_img = False         # Want to see the combined ramped images for every star?
show_centroids = False             # Print measured centroid on screen: True or False
show_disp = False                  # Show display of resulting positions? (will show 2 figs, same but different contrast)
Pier_corr = True                   # Include Pier's corrections to measured positions
tilt = False                       # Tilt angle: True or False

random_sample = False               # choose a random sample of stars from either detector: True or False
# control samples to be used when random is set to False
#stars_sample = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20]
#stars_sample = [1, 15, 60, 65, 67, 72, 81, 124, 132, 133, 139, 156, 166, 167, 182, 183, 187, 189, 198, 200]
stars_sample = [7, 24, 51, 56, 66, 68, 71, 72, 74, 91, 106, 109, 120, 125, 127, 128, 138, 154, 187, 188]
# Known bad stars in X and Y: 103, 105, 106, 112, 134, 152, 156, 170, 188
#6, 23, 50, 55, 65, 67, 70, 71, 73, 90, 105, 108, 119, 124, 126, 127, 137, 153, 186, 187

#######################################################################################################################


#  --> FUNCTIONS       
    
def writedatafile(show_centroids, data2write, lines4screenandfile):
    line0, line0a, line0b = lines4screenandfile
    save_text_file, output_file, st, bg, corr_cb_centroid_list, corr_true_center_centroid, loleftcoords, factor, mindiffs = data2write
    line1 = "{:<5} {:<10} {:<14} {:<16} {:<14} {:<16} {:<14} {:<18} {:<16} {:<16} {:<10} {:<10} {:<10.2f} {:<10}\n".format(
                                                    st, bg, 
                                                    corr_cb_centroid_list[0][0], corr_cb_centroid_list[0][1],
                                                    corr_cb_centroid_list[1][0], corr_cb_centroid_list[1][1],
                                                    corr_cb_centroid_list[2][0], corr_cb_centroid_list[2][1],
                                                    corr_true_center_centroid[0], corr_true_center_centroid[1],
                                                    loleftcoords[0], loleftcoords[1],
                                                    factor,
                                                    mindiffs)
    if save_text_file:
        f = open(output_file, "a")
        f.write(line1)
        f.close()
    if show_centroids:
        print(line0)
        print(line0a)
        print(line0b)
        print(line1) 
    

def runTEST(test2run, detectors, transf_direction, case, stars, P1P2data, bench_starP1, trueVsP1, trueVsP2):
    """ This function runs the test for both detectors and returns the results for the 20 star sample. 
    Pier_corr is set to False in case it was corrected before."""
    x13,y13, x23,y23, x15,y15, x25,y25, x17,y17, x27,y27 = P1P2data
    bench_V2P1, bench_V3P1 = trueVsP1
    bench_V2P2, bench_V3P2 = trueVsP2
    V2_3, V3_3, V2_5, V3_5, V2_7, V3_7 = [], [], [], [], [], [] 
    diffV2_3, diffV3_3, diffV2_5, diffV3_5, diffV2_7, diffV3_7 = [], [], [], [], [], []
    bench_V2_list, bench_V3_list = [], []
    if test2run == "T3":
        V2_13, V3_13, V2_15, V3_15, V2_17, V3_17 = [], [], [], [], [], [] 
        V2_23, V3_23, V2_25, V3_25, V2_27, V3_27 = [], [], [], [], [], [] 
        V2_3, V3_3 = [V2_13, V2_23], [V3_13, V3_23]
        V2_5, V3_5 = [V2_15, V2_25], [V3_15, V3_25]
        V2_7, V3_7 = [V2_17, V2_27], [V3_17, V3_27]
        diffV2_13, diffV3_13, diffV2_15, diffV3_15, diffV2_17, diffV3_17 = [], [], [], [], [], []
        diffV2_23, diffV3_23, diffV2_25, diffV3_25, diffV2_27, diffV3_27 = [], [], [], [], [], []
        diffV2_3, diffV3_3 = [diffV2_13, diffV2_23], [diffV3_13, diffV3_23]
        diffV2_5, diffV3_5 = [diffV2_15, diffV2_25], [diffV3_15, diffV3_25]
        diffV2_7, diffV3_7 = [diffV2_17, diffV2_27], [diffV3_17, diffV3_27]
        bench_V2_listP1, bench_V2_listP2 = [], []
        bench_V3_listP1, bench_V3_listP2 = [], []
        bench_V2_list, bench_V3_list = [bench_V2_listP1, bench_V2_listP2], [bench_V3_listP1, bench_V3_listP2]
    Vs = [V2_3, V3_3, V2_5, V3_5, V2_7, V3_7]
    diffs = [diffV2_3, diffV3_3, diffV2_5, diffV3_5, diffV2_7, diffV3_7]
    benchVs = [bench_V2_list, bench_V3_list]
    # Find the index at which to change detector
    change_detector_idx = len(stars)   # just in case all stars are from the same detector
    for st in stars:
        if st >= 100:
            if type(stars) is not list:
                change_detector_idx = stars.tolist().index(st)
            else:
                change_detector_idx = stars.index(st)
            break
    # slice arrays according to detector and run test
    # detector 492
    detector = detectors[1]   
    d2x13, d2y13 = x13[:change_detector_idx], y13[:change_detector_idx]
    d2x23, d2y23 = x23[:change_detector_idx], y23[:change_detector_idx]
    d2x15, d2y15 = x15[:change_detector_idx], y15[:change_detector_idx]
    d2x25, d2y25 = x25[:change_detector_idx], y25[:change_detector_idx]
    d2x17, d2y17 = x17[:change_detector_idx], y17[:change_detector_idx]
    d2x27, d2y27 = x27[:change_detector_idx], y27[:change_detector_idx]
    P1P2data = [d2x13, d2y13, d2x23, d2y23, d2x15, d2y15, d2x25, d2y25, d2x17, d2y17, d2x27, d2y27]
    d2bench_starP1 = bench_starP1[:change_detector_idx]
    d2bench_V2P1, d2bench_V3P1  = bench_V2P1[:change_detector_idx], bench_V3P1[:change_detector_idx]
    d2bench_V2P2, d2bench_V3P2  = bench_V2P2[:change_detector_idx], bench_V3P2[:change_detector_idx]
    d2benchV23 = [d2bench_V2P1, d2bench_V3P1, d2bench_V2P2, d2bench_V3P2]
    d2stars = stars[:change_detector_idx]
    data4test = [detector, transf_direction, case, d2stars, P1P2data, d2bench_starP1, d2benchV23]
    P1P2data, Vs, diffs, benchVs = runTest_and_append_results(test2run, data4test, Vs, diffs, benchVs)
    # detector 491
    detector = detectors[0]  
    if change_detector_idx != len(stars):   # in case all stars are from the same detector skip this part
        d1x13, d1y13 = x13[change_detector_idx:], y13[change_detector_idx:]
        d1x23, d1y23 = x23[change_detector_idx:], y23[change_detector_idx:]
        d1x15, d1y15 = x15[change_detector_idx:], y15[change_detector_idx:]
        d1x25, d1y25 = x25[change_detector_idx:], y25[change_detector_idx:]
        d1x17, d1y17 = x17[change_detector_idx:], y17[change_detector_idx:]
        d1x27, d1y27 = x27[change_detector_idx:], y27[change_detector_idx:]
        P1P2data = [d1x13, d1y13, d1x23, d1y23, d1x15, d1y15, d1x25, d1y25, d1x17, d1y17, d1x27, d1y27]
        d1bench_starP1 = bench_starP1[change_detector_idx:]
        d1bench_V2P1, d1bench_V3P1  = bench_V2P1[change_detector_idx:], bench_V3P1[change_detector_idx:]
        d1bench_V2P2, d1bench_V3P2  = bench_V2P2[change_detector_idx:], bench_V3P2[change_detector_idx:]
        d1benchV23 = [d1bench_V2P1, d1bench_V3P1, d1bench_V2P2, d1bench_V3P2]
        d1stars = stars[change_detector_idx:]
        data4test = [detector, transf_direction, case, d1stars, P1P2data, d1bench_starP1, d1benchV23]
        P1P2data, Vs, diffs, benchVs = runTest_and_append_results(test2run, data4test, Vs, diffs, benchVs)
    resultsTEST = [P1P2data, Vs, diffs, benchVs]
    return resultsTEST


def runTest_and_append_results(test2run, data4test, Vs, diffs, benchVs):
    """ This function runs the test for the specified detector and sliced arrays, and appends it to the results. """
    detector, transf_direction, case, stars, P1P2data, bench_starP1, benchV23 = data4test
    bench_V2P1, bench_V3P1, bench_V2P2, bench_V3P2 = benchV23
    avg_benchV2 = (bench_V2P1 + bench_V2P2)/2.0
    avg_benchV3 = (bench_V3P1 + bench_V3P2)/2.0
    avg_benchV23 = [avg_benchV2, avg_benchV3]
    T_V2_3, T_V3_3, T_V2_5, T_V3_5, T_V2_7, T_V3_7 = Vs
    T_diffV2_3, T_diffV3_3, T_diffV2_5, T_diffV3_5, T_diffV2_7, T_diffV3_7 = diffs
    Tbench_V2_list, Tbench_V3_list = benchVs
    if test2run == "T1":
        transformations, diffs, benchVs_list = TEST1(detector, transf_direction, stars, case, bench_starP1, avg_benchV23, P1P2data)
    if test2run == "T2":
        transformations, diffs, benchVs_list = TEST2(detector, transf_direction, stars, case, bench_starP1, avg_benchV23, P1P2data)
    if test2run == "T3":
        transformations, diffs, benchVs_list = TEST3(detector, transf_direction, stars, case, bench_starP1, benchV23, P1P2data)
    # append appropriate number of arrays
    if test2run == "T1" or test2run == "T2":
        V2_3, V3_3, V2_5, V3_5, V2_7, V3_7 = transformations
        diffV2_3, diffV3_3, diffV2_5, diffV3_5, diffV2_7, diffV3_7 = diffs
        bench_V2_list, bench_V3_list = benchVs_list
        for v23, v33, v25, v35, v27, v37 in zip(V2_3, V3_3, V2_5, V3_5, V2_7, V3_7):
            T_V2_3.append(v23)
            T_V3_3.append(v33)
            T_V2_5.append(v25)
            T_V3_5.append(v35)
            T_V2_7.append(v27)
            T_V3_7.append(v37)
        for dv23, dv33, dv25, dv35, dv27, dv37 in zip(diffV2_3, diffV3_3, diffV2_5, diffV3_5, diffV2_7, diffV3_7):
            T_diffV2_3.append(dv23)
            T_diffV3_3.append(dv33)
            T_diffV2_5.append(dv25)
            T_diffV3_5.append(dv35)
            T_diffV2_7.append(dv27)
            T_diffV3_7.append(dv37)
        for bv2, bv3 in zip(bench_V2_list, bench_V3_list):
            Tbench_V2_list.append(bv2)
            Tbench_V3_list.append(bv3)
    if test2run == "T3":
        # unfold results from test 3
        transformationsP1, transformationsP2 = transformations
        V2_13, V3_13, V2_15, V3_15, V2_17, V3_17 = transformationsP1
        V2_23, V3_23, V2_25, V3_25, V2_27, V3_27 = transformationsP2
        diffsP1, diffsP2 = diffs
        diffV2_13, diffV3_13, diffV2_15, diffV3_15, diffV2_17, diffV3_17 = diffsP1
        diffV2_23, diffV3_23, diffV2_25, diffV3_25, diffV2_27, diffV3_27 = diffsP2
        bench_V2_listP1, bench_V3_listP1, bench_V2_listP2, bench_V3_listP2 = benchVs_list
        # unfold empty lists to append to 
        T_V2_13, T_V2_23 = T_V2_3
        T_V3_13, T_V3_23 = T_V3_3
        T_V2_15, T_V2_25 = T_V2_5
        T_V3_15, T_V3_25 = T_V3_5
        T_V2_17, T_V2_27 = T_V2_7
        T_V3_17, T_V3_27 = T_V3_7
        T_diffV2_13, T_diffV2_23 = T_diffV2_3
        T_diffV3_13, T_diffV3_23 = T_diffV3_3
        T_diffV2_15, T_diffV2_25 = T_diffV2_5
        T_diffV3_15, T_diffV3_25 = T_diffV3_5
        T_diffV2_17, T_diffV2_27 = T_diffV2_7
        T_diffV3_17, T_diffV3_27 = T_diffV3_7
        Tbench_V2_listP1, Tbench_V2_listP2 = Tbench_V2_list
        Tbench_V3_listP1, Tbench_V3_listP2 = Tbench_V3_list
        # append to individual position lists       
        for v23, v33, v25, v35, v27, v37 in zip(V2_13, V3_13, V2_15, V3_15, V2_17, V3_17):
            T_V2_13.append(v23)
            T_V3_13.append(v33)
            T_V2_15.append(v25)
            T_V3_15.append(v35)
            T_V2_17.append(v27)
            T_V3_17.append(v37)
        for v23, v33, v25, v35, v27, v37 in zip(V2_23, V3_23, V2_25, V3_25, V2_27, V3_27):
            T_V2_23.append(v23)
            T_V3_23.append(v33)
            T_V2_25.append(v25)
            T_V3_25.append(v35)
            T_V2_27.append(v27)
            T_V3_27.append(v37)
        T_V2_3, T_V3_3 = [T_V2_13, T_V2_23], [T_V3_13, T_V3_23]
        T_V2_5, T_V3_5 = [T_V2_15, T_V2_25], [T_V3_15, T_V3_25]
        T_V2_7, T_V3_7 = [T_V2_17, T_V2_27], [T_V3_17, T_V3_27]
        for dv23, dv33, dv25, dv35, dv27, dv37 in zip(diffV2_13, diffV3_13, diffV2_15, diffV3_15, diffV2_17, diffV3_17):
            T_diffV2_13.append(dv23)
            T_diffV3_13.append(dv33)
            T_diffV2_15.append(dv25)
            T_diffV3_15.append(dv35)
            T_diffV2_17.append(dv27)
            T_diffV3_17.append(dv37)
        for dv23, dv33, dv25, dv35, dv27, dv37 in zip(diffV2_23, diffV3_23, diffV2_25, diffV3_25, diffV2_27, diffV3_27):
            T_diffV2_23.append(dv23)
            T_diffV3_23.append(dv33)
            T_diffV2_25.append(dv25)
            T_diffV3_25.append(dv35)
            T_diffV2_27.append(dv27)
            T_diffV3_27.append(dv37)
        T_diffV2_3, T_diffV3_3 = [T_diffV2_13, T_diffV2_23], [T_diffV3_13, T_diffV3_23]
        T_diffV2_5, T_diffV3_5 = [T_diffV2_15, T_diffV2_25], [T_diffV3_15, T_diffV3_25]
        T_diffV2_7, T_diffV3_7 = [T_diffV2_17, T_diffV2_27], [T_diffV3_17, T_diffV3_27]
        for bv2, bv3 in zip(bench_V2_listP1, bench_V3_listP1):
            Tbench_V2_listP1.append(bv2)
            Tbench_V3_listP1.append(bv3)
        for bv2, bv3 in zip(bench_V2_listP2, bench_V3_listP2):
            Tbench_V2_listP2.append(bv2)
            Tbench_V3_listP2.append(bv3)
        Tbench_V2_list, Tbench_V3_list = [Tbench_V2_listP1, Tbench_V2_listP2], [Tbench_V3_listP1, Tbench_V3_listP2]
    Vs = [T_V2_3, T_V3_3, T_V2_5, T_V3_5, T_V2_7, T_V3_7]
    diffs = [T_diffV2_3, T_diffV3_3, T_diffV2_5, T_diffV3_5, T_diffV2_7, T_diffV3_7]
    benchVs = [Tbench_V2_list, Tbench_V3_list]
    return P1P2data, Vs, diffs, benchVs


def TEST1(detector, transf_direction, stars, case, bench_starP1, avg_benchV23, P1P2data):
    # TEST 1: (a) Avg P1 and P2, (b) transform to V2-V3, (c) compare to avg reference positions (V2-V3 space)
    avg_benchV2, avg_benchV3 = avg_benchV23
    x13,y13, x23,y23, x15,y15, x25,y25, x17,y17, x27,y27 = P1P2data
    # Step (a) - averages
    avgx3 = (x13+x23)/2.0
    avgy3 = (y13+y23)/2.0
    avgx5 = (x15+x25)/2.0
    avgy5 = (y15+y25)/2.0
    avgx7 = (x17+x27)/2.0
    avgy7 = (y17+y27)/2.0
    # Step (b) - transformations to degrees
    T1_V2_3, T1_V3_3 = ct.coords_transf(transf_direction, detector, filter_input, avgx3, avgy3, tilt, debug)
    T1_V2_5, T1_V3_5 = ct.coords_transf(transf_direction, detector, filter_input, avgx5, avgy5, tilt, debug)
    T1_V2_7, T1_V3_7 = ct.coords_transf(transf_direction, detector, filter_input, avgx7, avgy7, tilt, debug)
    # TEST 1: (a) Avg P1 and P2, (b) transform to V2-V3, (c) compare to avg reference positions (V2-V3 space)
    # Step (c) - comparison
    T1_diffV2_3, T1_diffV3_3, T1bench_V2_list, T1bench_V3_list = tf.compare2ref(case, bench_starP1, avg_benchV2, avg_benchV3, stars, T1_V2_3, T1_V3_3, arcsecs=diffs_in_arcsecs)
    T1_diffV2_5, T1_diffV3_5, _, _ = tf.compare2ref(case, bench_starP1, avg_benchV2, avg_benchV3, stars, T1_V2_5, T1_V3_5, arcsecs=diffs_in_arcsecs)
    T1_diffV2_7, T1_diffV3_7, _, _ = tf.compare2ref(case, bench_starP1, avg_benchV2, avg_benchV3, stars, T1_V2_7, T1_V3_7, arcsecs=diffs_in_arcsecs)
    if debug:
        print ("TEST 1: ")
        print ("transformations: detector (avgx, avgy),  sky (V2, V3),  true (avgV2, avgV3)")
        print ("            ChBx3: ", avgx3[0], avgy3[0], T1_V2_3[0], T1_V3_3[0], avg_benchV2[0], avg_benchV3[0])
        print ("            ChBx5: ", avgx5[0], avgy5[0], T1_V2_5[0], T1_V3_5[0], avg_benchV2[0], avg_benchV3[0])
        print ("            ChBx7: ", avgx7[0], avgy7[0], T1_V2_7[0], T1_V3_7[0], avg_benchV2[0], avg_benchV3[0])
        raw_input(" * press enter to continue... \n")
    # Organize results
    T1_transformations = [T1_V2_3, T1_V3_3, T1_V2_5, T1_V3_5, T1_V2_7, T1_V3_7]
    T1_diffs = [T1_diffV2_3, T1_diffV3_3, T1_diffV2_5, T1_diffV3_5, T1_diffV2_7, T1_diffV3_7]
    T1_benchVs_list = [T1bench_V2_list, T1bench_V3_list]
    return T1_transformations, T1_diffs, T1_benchVs_list
    
def TEST2(detector, transf_direction, stars, case, bench_starP1, avg_benchV23, P1P2data):
    # TEST 2: (a) Transform individual P1 and P2 to V2-V3, (b) avg V2-V3 space positions, (c) compare to avg reference positions
    x13, y13, x23, y23, x15, y15, x25, y25, x17, y17, x27, y27 = P1P2data
    avg_benchV2, avg_benchV3 = avg_benchV23
    # Step (a) - transformations
    T2_V2_13, T2_V3_13 = ct.coords_transf(transf_direction, detector, filter_input, x13, y13, tilt, debug)
    T2_V2_15, T2_V3_15 = ct.coords_transf(transf_direction, detector, filter_input, x15, y15, tilt, debug)
    T2_V2_17, T2_V3_17 = ct.coords_transf(transf_direction, detector, filter_input, x17, y17, tilt, debug)
    T2_V2_23, T2_V3_23 = ct.coords_transf(transf_direction, detector, filter_input, x23, y23, tilt, debug)
    T2_V2_25, T2_V3_25 = ct.coords_transf(transf_direction, detector, filter_input, x25, y25, tilt, debug)
    T2_V2_27, T2_V3_27 = ct.coords_transf(transf_direction, detector, filter_input, x27, y27, tilt, debug)
    # Step (b) - averages
    T2_V2_3 = (T2_V2_13 + T2_V2_23)/2.0
    T2_V3_3 = (T2_V3_13 + T2_V3_23)/2.0
    T2_V2_5 = (T2_V2_15 + T2_V2_25)/2.0
    T2_V3_5 = (T2_V3_15 + T2_V3_25)/2.0
    T2_V2_7 = (T2_V2_17 + T2_V2_27)/2.0
    T2_V3_7 = (T2_V3_17 + T2_V3_27)/2.0
    # Step (c) - comparison
    T2_diffV2_3, T2_diffV3_3, T2bench_V2_list, T2bench_V3_list = tf.compare2ref(case, bench_starP1, avg_benchV2, avg_benchV3, stars, T2_V2_3, T2_V3_3, arcsecs=diffs_in_arcsecs)
    T2_diffV2_5, T2_diffV3_5, _, _ = tf.compare2ref(case, bench_starP1, avg_benchV2, avg_benchV3, stars, T2_V2_5, T2_V3_5, arcsecs=diffs_in_arcsecs)
    T2_diffV2_7, T2_diffV3_7, _, _ = tf.compare2ref(case, bench_starP1, avg_benchV2, avg_benchV3, stars, T2_V2_7, T2_V3_7, arcsecs=diffs_in_arcsecs)
    if debug:
        print ("TEST 2: ")
        print ("transformations: detector P1 and P2 (x, y),  sky (avgV2, avgV3),  true (avgV2, avgV3)")
        print ("            ChBx3: ", x13[0], y13[0], x23[0], y23[0], T2_V2_3[0], T2_V3_3[0], avg_benchV2[0], avg_benchV3[0])
        print ("            ChBx5: ", x15[0], y15[0], x25[0], y25[0], T2_V2_5[0], T2_V3_5[0], avg_benchV2[0], avg_benchV3[0])
        print ("            ChBx7: ", x17[0], y17[0], x27[0], y27[0], T2_V2_7[0], T2_V3_7[0], avg_benchV2[0], avg_benchV3[0])
        raw_input(" * press enter to continue... \n")
    # Organize results
    T2_transformations = [T2_V2_3, T2_V3_3, T2_V2_5, T2_V3_5, T2_V2_7, T2_V3_7]
    T2_diffs = [T2_diffV2_3, T2_diffV3_3, T2_diffV2_5, T2_diffV3_5, T2_diffV2_7, T2_diffV3_7]
    T2_benchVs_list = [T2bench_V2_list, T2bench_V3_list]
    return T2_transformations, T2_diffs, T2_benchVs_list
    

def TEST3(detector, transf_direction, stars, case, bench_starP1, bench_Vs, P1P2data):
    # TEST 3: (a) Transform P1 and P2 individually to V2-V3 (b) compare star by star and position by position
    x13, y13, x23, y23, x15, y15, x25, y25, x17, y17, x27, y27 = P1P2data
    bench_V2P1, bench_V3P1, bench_V2P2, bench_V3P2 = bench_Vs
    # Step (a) - transformations
    T3_V2_13, T3_V3_13 = ct.coords_transf(transf_direction, detector, filter_input, x13, y13, tilt, debug)
    T3_V2_15, T3_V3_15 = ct.coords_transf(transf_direction, detector, filter_input, x15, y15, tilt, debug)
    T3_V2_17, T3_V3_17 = ct.coords_transf(transf_direction, detector, filter_input, x17, y17, tilt, debug)
    T3_V2_23, T3_V3_23 = ct.coords_transf(transf_direction, detector, filter_input, x23, y23, tilt, debug)
    T3_V2_25, T3_V3_25 = ct.coords_transf(transf_direction, detector, filter_input, x25, y25, tilt, debug)
    T3_V2_27, T3_V3_27 = ct.coords_transf(transf_direction, detector, filter_input, x27, y27, tilt, debug)
    # Step (b) - comparison
    T3_diffV2_13, T3_diffV3_13, T3bench_V2_listP1, T3bench_V3_listP1 = tf.compare2ref(case, bench_starP1, bench_V2P1, bench_V3P1, stars, T3_V2_13, T3_V3_13, arcsecs=diffs_in_arcsecs)
    T3_diffV2_23, T3_diffV3_23, T3bench_V2_listP2, T3bench_V3_listP2 = tf.compare2ref(case, bench_starP1, bench_V2P2, bench_V3P2, stars, T3_V2_23, T3_V3_23, arcsecs=diffs_in_arcsecs)
    T3_diffV2_15, T3_diffV3_15, _, _ = tf.compare2ref(case, bench_starP1, bench_V2P1, bench_V3P1, stars, T3_V2_15, T3_V3_15, arcsecs=diffs_in_arcsecs)
    T3_diffV2_25, T3_diffV3_25, _, _ = tf.compare2ref(case, bench_starP1, bench_V2P2, bench_V3P2, stars, T3_V2_25, T3_V3_25, arcsecs=diffs_in_arcsecs)
    T3_diffV2_17, T3_diffV3_17, _, _ = tf.compare2ref(case, bench_starP1, bench_V2P1, bench_V3P1, stars, T3_V2_17, T3_V3_17, arcsecs=diffs_in_arcsecs)
    T3_diffV2_27, T3_diffV3_27, _, _ = tf.compare2ref(case, bench_starP1, bench_V2P2, bench_V3P2, stars, T3_V2_27, T3_V3_27, arcsecs=diffs_in_arcsecs)
    if debug:
        print ("TEST 3: ")
        print ("transformations: detector P1 and P2 (x, y),  sky P1 and P2 (V2, V3),  true P1 and P2 (V2, V3)")
        print ("            ChBx3: ", x13[0], y13[0], x23[0], y23[0], T3_V2_13[0], T3_V3_13[0], T3_V2_23[0], T3_V3_23[0], bench_V2P1[0], bench_V3P1[0], bench_V2P2[0], bench_V3P2[0])
        print ("            ChBx5: ", x15[0], y15[0], x25[0], y25[0], T3_V2_13[0], T3_V3_13[0], T3_V2_23[0], T3_V3_23[0], bench_V2P1[0], bench_V3P1[0], bench_V2P2[0], bench_V3P2[0])
        print ("            ChBx7: ", x17[0], y17[0], x27[0], y27[0], T3_V2_13[0], T3_V3_13[0], T3_V2_23[0], T3_V3_23[0], bench_V2P1[0], bench_V3P1[0], bench_V2P2[0], bench_V3P2[0])
        raw_input(" * press enter to continue... \n")
    # Organize results
    T3_transformationsP1 = [T3_V2_13, T3_V3_13, T3_V2_15, T3_V3_15, T3_V2_17, T3_V3_17]
    T3_transformationsP2 = [T3_V2_23, T3_V3_23, T3_V2_25, T3_V3_25, T3_V2_27, T3_V3_27]
    T3_transformations = [T3_transformationsP1, T3_transformationsP2]
    T3_diffsP1 = [T3_diffV2_13, T3_diffV3_13, T3_diffV2_15, T3_diffV3_15, T3_diffV2_17, T3_diffV3_17]
    T3_diffsP2 = [T3_diffV2_23, T3_diffV3_23, T3_diffV2_25, T3_diffV3_25, T3_diffV2_27, T3_diffV3_27]
    T3_diffs = [T3_diffsP1, T3_diffsP2]
    T3_benchVs_list = [T3bench_V2_listP1, T3bench_V3_listP1, T3bench_V2_listP2, T3bench_V3_listP2]
    return T3_transformations, T3_diffs, T3_benchVs_list


def convert2MSAcenter(xin, yin, xtin, ytin):
    """ This function converts the measured coordinates of each star into the frame relative
    to the center of the MSA. """
    # Center coordinates of NIRSpec V2, V3 
    x0_XAN = 376.769       # V2 in arcsec
    y0_YAN = -428.453      # V3 in arcsec
    # measured V2 V3 in arcsec
    x0 = x0_XAN/3600.      # conversion of V2 to XAN in degrees
    y0_YANd = y0_YAN/3600. # intermediate conversion: V3 arcsec to V3 degrees=-0.119015
    y0 = -y0_YANd -0.13    # convert V3 degrees to YAN=+0.249015 
    # convert inputs to MSA center
    x = xin - x0
    y = yin - y0
    xt = xtin - x0
    yt = ytin - y0
    return x, y, xt, yt


def get_stats(case, T_transformations, T_diffs, T_benchVs_list, Nsigma, max_iterations):
    """ This function obtains the standard deviations through regular statistics as well as through
    a sigma clipping algorithm and an iterative least square algorithm. It also obtains the minimum
    differences from checkbox sizes 3, 5, and 7, and returns the counter for each."""
    T_V2_3, T_V3_3, T_V2_5, T_V3_5, T_V2_7, T_V3_7 = T_transformations
    T_V2_3, T_V3_3 = np.array(T_V2_3), np.array(T_V3_3)
    T_V2_5, T_V3_5 = np.array(T_V2_5), np.array(T_V3_5)
    T_V2_7, T_V3_7 = np.array(T_V2_7), np.array(T_V3_7)
    T_diffV2_3, T_diffV3_3, T_diffV2_5, T_diffV3_5, T_diffV2_7, T_diffV3_7 = T_diffs
    T_diffV2_3, T_diffV3_3 = np.array(T_diffV2_3), np.array(T_diffV3_3)
    T_diffV2_5, T_diffV3_5 = np.array(T_diffV2_5), np.array(T_diffV3_5)
    T_diffV2_7, T_diffV3_7 = np.array(T_diffV2_7), np.array(T_diffV3_7)
    Tbench_V2_list, Tbench_V3_list = T_benchVs_list
    # calculate standard deviations and means
    Tstdev_V2_3, Tmean_V2_3 = tf.find_std(T_diffV2_3)
    Tstdev_V2_5, Tmean_V2_5 = tf.find_std(T_diffV2_5)
    Tstdev_V2_7, Tmean_V2_7 = tf.find_std(T_diffV2_7)
    Tstdev_V3_3, Tmean_V3_3 = tf.find_std(T_diffV3_3)
    Tstdev_V3_5, Tmean_V3_5 = tf.find_std(T_diffV3_5)
    Tstdev_V3_7, Tmean_V3_7 = tf.find_std(T_diffV3_7)
    Tbench_V2, Tbench_V3 = np.array(Tbench_V2_list), np.array(Tbench_V3_list)
    # get the minimum of the differences
    T_min_diff, T_counter = tf.get_mindiff(T_diffV2_3, T_diffV2_5, T_diffV2_7)
    # calculate least squares but first convert to MSA center
    T_V2_3, T_V3_3, Tbench_V2, Tbench_V3 = convert2MSAcenter(T_V2_3, T_V3_3, Tbench_V2, Tbench_V3)
    T_V2_5, T_V3_5, _, _ = convert2MSAcenter(T_V2_5, T_V3_5, Tbench_V2, Tbench_V3)
    T_V2_7, T_V3_7, _, _ = convert2MSAcenter(T_V2_7, T_V3_7, Tbench_V2, Tbench_V3)
    # to express in arcsecs multiply by 3600.0
    TLSdeltas_3, TLSsigmas_3, TLSlines2print_3, rejected_elements_3 = lsi.ls_fit_iter(max_iterations, T_V2_3*3600.0, T_V3_3*3600.0, Tbench_V2*3600.0, Tbench_V3*3600.0)
    TLSdeltas_5, TLSsigmas_5, TLSlines2print_5, rejected_elements_5 = lsi.ls_fit_iter(max_iterations, T_V2_5*3600.0, T_V3_5*3600.0, Tbench_V2*3600.0, Tbench_V3*3600.0)
    TLSdeltas_7, TLSsigmas_7, TLSlines2print_7, rejected_elements_7 = lsi.ls_fit_iter(max_iterations, T_V2_7*3600.0, T_V3_7*3600.0, Tbench_V2*3600.0, Tbench_V3*3600.0)
    # Do N-sigma rejection
    TsigmaV2_3, TmeanV2_3, TsigmaV3_3, TmeanV3_3, TnewV2_3, TnewV3_3, Tniter_3, Tlines2print_3, rej_elements_3 = tf.Nsigma_rejection(Nsigma, T_diffV2_3, T_diffV3_3, max_iterations)
    TsigmaV2_5, TmeanV2_5, TsigmaV3_5, TmeanV3_5, TnewV2_5, TnewV3_5, Tniter_5, Tlines2print_5, rej_elements_5 = tf.Nsigma_rejection(Nsigma, T_diffV2_5, T_diffV3_5, max_iterations)
    TsigmaV2_7, TmeanV2_7, TsigmaV3_7, TmeanV3_7, TnewV2_7, TnewV3_7, Tniter_7, Tlines2print_7, rej_elements_7 = tf.Nsigma_rejection(Nsigma, T_diffV2_7, T_diffV3_7, max_iterations)
    # organize the results
    st_devsAndMeans = [Tstdev_V2_3, Tmean_V2_3, Tstdev_V2_5, Tmean_V2_5, Tstdev_V2_7, Tmean_V2_7,
                       Tstdev_V3_3, Tmean_V3_3, Tstdev_V3_5, Tmean_V3_5, Tstdev_V3_7, Tmean_V3_7]
    diff_counter = [T_min_diff, T_counter]
    bench_values = [Tbench_V2, Tbench_V3]
    sigmas_deltas = [TLSdeltas_3, TLSsigmas_3, TLSlines2print_3, 
                     TLSdeltas_5, TLSsigmas_5, TLSlines2print_5,
                     TLSdeltas_7, TLSsigmas_7, TLSlines2print_7]
    sigma_reject = [TsigmaV2_3, TmeanV2_3, TsigmaV3_3, TmeanV3_3, TnewV2_3, TnewV3_3, Tniter_3, Tlines2print_3,
                    TsigmaV2_5, TmeanV2_5, TsigmaV3_5, TmeanV3_5, TnewV2_5, TnewV3_5, Tniter_5, Tlines2print_5,
                    TsigmaV2_7, TmeanV2_7, TsigmaV3_7, TmeanV3_7, TnewV2_7, TnewV3_7, Tniter_7, Tlines2print_7]
    rejected_elementsLS = [rejected_elements_3, rejected_elements_5, rejected_elements_7]
    rejected_elementsNsig = [rej_elements_3, rej_elements_5, rej_elements_7]
    results_stats = [st_devsAndMeans, diff_counter, bench_values, sigmas_deltas, sigma_reject, rejected_elementsLS, rejected_elementsNsig]
    return results_stats


def combine2arrays(arr1, arr2, combined_arr):
    for item in arr1:
        combined_arr = np.append(combined_arr, item)
    for item in arr2:
        combined_arr = np.append(combined_arr, item)
    return combined_arr

def get_rejected_stars(stars_samplex2, rejected_elements_idx):
    rejected_elements = []
    if len(rejected_elements_idx) != len(stars_samplex2):
        for i in rejected_elements_idx:
            rejected_elements.append(stars_samplex2[i])
    return rejected_elements


def print_results(stars_sample, case, test2perform, diffs_in_arcsecs, Tstdev_Vs, Tmean_Vs, 
                  T_diff_counter, save_text_file, TLSlines2print, Tlines2print, Tbench_Vs_list, 
                  T_Vs, T_diffVs, rejected_elementsLS, rejected_eleNsig):
    """ This function writes to text the results from the Test performed. """
    # unfold variables
    Tstdev_V2_3, Tstdev_V3_3, Tstdev_V2_5, Tstdev_V3_5, Tstdev_V2_7, Tstdev_V3_7 = Tstdev_Vs
    Tmean_V2_3, Tmean_V3_3, Tmean_V2_5, Tmean_V3_5, Tmean_V2_7, Tmean_V3_7 = Tmean_Vs
    T_min_diff, T_counter = T_diff_counter
    TLSlines2print_3, TLSlines2print_5, TLSlines2print_7 = TLSlines2print
    Tlines2print_3, Tlines2print_5, Tlines2print_7 = Tlines2print
    Tbench_V2_list, Tbench_V3_list = Tbench_Vs_list
    T_V2_3, T_V3_3, T_V2_5, T_V3_5, T_V2_7, T_V3_7 = T_Vs
    #T_diffV2_3, T_diffV3_3, T_diffV2_5, T_diffV3_5, T_diffV2_7, T_diffV3_7 = T_diffVs 
    rejected_elements_idx3, rejected_elements_idx5, rejected_elements_idx7 = rejected_elementsLS 
    Nsigrej_elements_idx3, Nsigrej_elements_idx5, Nsigrej_elements_idx7 = rejected_elementsLS 
    # define lines to print
    line0 = "{}".format("Differences = diffs = True_Positions - Measured_Positions")
    if diffs_in_arcsecs:
        line0bis = "{}".format("*** diffs are in units of arcsecs")
    else:
        line0bis = "{}".format("*** diffs are in units of degrees")
    if test2perform == "T1":
        line1 = "{}\n {}".format("Test1: average P1 and P2, transform to V2-V3, calculate differences",
                                 "  * Standard deviations and means ")
    if test2perform == "T2":
        line1 = "{}".format("Test2: P1 P2, average positions in V2-V3, calculate differences")
    if test2perform == "T3":
        line1 = "{}".format("Test3: P1 and P2, transform to V2-V3 space individually, calculate differences position to position")
    # print regular standard deviations and means
    line2a = "std_dev_V2_3 = {:<20}    std_dev_V3_3 = {:<20}".format(Tstdev_V2_3, Tstdev_V3_3)
    line2b = "std_dev_V2_5 = {:<20}    std_dev_V3_5 = {:<20}".format(Tstdev_V2_5, Tstdev_V3_5)
    line2c = "std_dev_V2_7 = {:<20}    std_dev_V3_7 = {:<20}".format(Tstdev_V2_7, Tstdev_V3_7)
    line3a = "   mean_V2_3 = {:<22}     mean_V3_3 = {:<22}".format(Tmean_V2_3, Tmean_V3_3)
    line3b = "   mean_V2_5 = {:<22}     mean_V3_5 = {:<22}".format(Tmean_V2_5, Tmean_V3_5)
    line3c = "   mean_V2_7 = {:<22}     mean_V3_7 = {:<22}".format(Tmean_V2_7, Tmean_V3_7)
    # Print rejected stars for least squares and N-sigma rejection
    stars_samplex2 = []
    for _ in range(2):
        for st_i in stars_sample:
            stars_samplex2.append(st_i) 
    rejected_elements_3 = get_rejected_stars(stars_samplex2, rejected_elements_idx3)
    rejected_elements_5 = get_rejected_stars(stars_samplex2, rejected_elements_idx5)
    rejected_elements_7 = get_rejected_stars(stars_samplex2, rejected_elements_idx7)
    Nsig_rej_elements_3 = get_rejected_stars(stars_samplex2, Nsigrej_elements_idx3)
    Nsig_rej_elements_5 = get_rejected_stars(stars_samplex2, Nsigrej_elements_idx5)
    Nsig_rej_elements_7 = get_rejected_stars(stars_samplex2, Nsigrej_elements_idx7)
    line3bisAa = "- Rejected stars -"
    line3bisAb = "   checkbox3: {} ".format(rejected_elements_3)
    line3bisAc = "   checkbox5: {} ".format(rejected_elements_5)
    line3bisAd = "   checkbox7: {} ".format(rejected_elements_7)
    line3bisAe = "   checkbox3: {} ".format(Nsig_rej_elements_3)
    line3bisAf = "   checkbox5: {} ".format(Nsig_rej_elements_5)
    line3bisAg = "   checkbox7: {} ".format(Nsig_rej_elements_7)
    # Print number of repetitions to find best checkbox
    line3bisB = "\n *** Repetitions Diffs: {}".format(T_counter)
    line4 = "{:<5} {:<20} {:<40} {:<40} {:<38} {:<28} {:<7}".format(
                    "Star", "BG_value", "Pos_Checkbox_3", "Pos_Checkbox_5", "PosCheckbox_7",
                    "True_Pos", "MinDiff")
    line5 = "{:>10} {:>15} {:>17} {:>22} {:>17} {:>22} {:>22} {:>17} {:>17}".format(background_method,
                    "V2", "V3", "V2", "V3", "V2", "V3", "V2", "V3")
    print (line0)
    print (line0bis)
    print (line1)
    print (line2a)
    print (line2b)
    print (line2c)
    print (line3a)
    print (line3b)
    print (line3c)
    print (line3bisAa)
    print (" From least square routine: ")
    print (line3bisAb)
    print (line3bisAc)
    print (line3bisAd)
    print (" From N-sigma rejection routine: ")
    print (line3bisAe)
    print (line3bisAf)
    print (line3bisAg)
    print (line3bisB)
    print (line4)
    print (line5)
    if save_text_file:
        txt_out = path4results+test2perform+"_results_"+case+".txt"
        to = open(txt_out, "w+")
        to.write(line0+"\n")
        to.write(line0bis+"\n")
        to.write(line1+"\n")
        to.write(line2a+"\n")
        to.write(line2b+"\n")
        to.write(line2c+"\n")
        to.write(line3a+"\n")
        to.write(line3b+"\n")
        to.write(line3c+"\n")
        # print standard deviations from least squares routine
        to.write("\n * From least squares routine:  \n")
        to.write(line3bisAa+"\n")
        to.write(line3bisAb+"\n")
        to.write(line3bisAc+"\n")
        to.write(line3bisAd+"\n")        
        to.write("       Checkbox 3:  \n")
        for line2print in TLSlines2print_3:
            to.write(line2print+"\n")
        to.write("       Checkbox 5:  \n")
        for line2print in TLSlines2print_5:
            to.write(line2print+"\n")
        to.write("       Checkbox 7:  \n")
        for line2print in TLSlines2print_7:
            to.write(line2print+"\n")
        # print standard deviations and means after n-sigma rejection
        to.write("\n * From N-sigma rejection routine:  \n")
        to.write(line3bisAa+"\n")
        to.write(line3bisAe+"\n")
        to.write(line3bisAf+"\n")
        to.write(line3bisAg+"\n")        
        to.write(" Checkbox 3:  \n")
        for line2print in Tlines2print_3:
            to.write(line2print+"\n")
        to.write(" Checkbox 5:  \n")
        for line2print in Tlines2print_5:
            to.write(line2print+"\n")
        to.write(" Checkbox 7:  \n")
        for line2print in Tlines2print_7:
            to.write(line2print+"\n")
        to.write(line3bisB+"\n")
        to.write(line4+"\n")
        to.write(line5+"\n")
    j = 0
    for i, _ in enumerate(T_V2_3):
        st = int(stars_sample[j])
        line6 = "{:<5} {:<5} {:>20}  {:<20} {:>18}  {:<20} {:>18}  {:<20} {:>17}  {:<17}  {:>5}".format(
                    st, background2use, 
                    T_V2_3[i], T_V3_3[i], T_V2_5[i], T_V3_5[i], T_V2_7[i], T_V3_7[i], 
                    Tbench_V2_list[i], Tbench_V3_list[i],
                    T_min_diff[i])
        print (line6)
        if save_text_file:
            to.write(line6+"\n")
        j = j + 1
        if st == stars_sample[-1]:
            j = 0
    if save_text_file:
        to.close()
        print (" * Results saved in file: ", txt_out)


def remove_bad_stars(stars_sample):
        """ This function reads the text files of bad stars, compares the sample data, removes 
        the bad stars, and returns the sample without bad stars. """
        # paths to files
        scene1_bad_stars_file = os.path.abspath("../bad_stars/scene1_bad_stars.txt")
        scene2_bad_stars_file = os.path.abspath("../bad_stars/scene2_bad_stars.txt")
        # read files ad get lists
        scene1_bad_stars = np.loadtxt(scene1_bad_stars_file, comments="#", skiprows=2, unpack=True)
        scene2_bad_stars = np.loadtxt(scene2_bad_stars_file, comments="#", skiprows=2, unpack=True)
        print ("There are %i bad stars in Scenario 1 and %i bad stars in Scenario 2." % (len(scene1_bad_stars), 
                                                                                         len(scene2_bad_stars)))
        # compare to stars_sample
        for st_sam in stars_sample:
            if st_sam in scene1_bad_stars or st_sam in scene2_bad_stars:
                stars_sample.pop()
        return stars_sample
        

#######################################################################################################################

#  --> CODE

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
        if keep_bad_stars == False:
            remove_bad_stars(stars_sample)
# remove the bad stars
if keep_bad_stars == False:
    remove_bad_stars(stars_sample)
    print (" * Sample has %i stars left." % len(stars_sample))

# order the star list 
stars_sample.sort(key=lambda xx: xx)
print ("stars_sample =", stars_sample)
raw_input("\nPress enter to continue...")

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
benchmark_data, magnitudes = tf.read_star_param_files(scene2study)
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
        f.write(line0a+"\n")
        f.write(line0b+"\n")
        f.close()

# get the star files to run the TA algorithm on
dir2test_list = tf.get_raw_star_directory(path4starfiles, scene, shutters, noise)

# run centroid algorithm on each position and save them into a text file
x13, x15, x17 = np.array([]), np.array([]), np.array([])
y13, y15, y17 = np.array([]), np.array([]), np.array([])
x23, x25, x27 = np.array([]), np.array([]), np.array([])
y23, y25, y27 = np.array([]), np.array([]), np.array([])
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
                if star_exists == False:
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
                # Do background correction on each of 3 ramp images
                master_img_bgcorr = jtl.bg_correction(master_img, bg_method=background_method, 
                                                      bg_value=bg_value, bg_frac=bg_frac)
                # Obtain the combined FITS image that combines all frames into one image AND
                # check if all image is zeros, take the image that still has a max value
                psf = tu.readimage(master_img_bgcorr, debug=debug)                
                cb_centroid_list = tf.run_recursive_centroids(psf, bg_frac, xwidth_list, ywidth_list, 
                                                           checkbox_size, max_iter, threshold, 
                                                           determine_moments, debug, display_master_img, vlim=vlim)
                cb_centroid_list, loleftcoords, true_center32x32, differences_true_TA = tf.centroid2fulldetector(cb_centroid_list, 
                                                                                                    true_center)
                if output_full_detector == False:
                    true_center = true_center32x32
                # Correct true centers for average value given by Pier 
                if Pier_corr:
                    corr_cb_centroid_list = tf.do_Piers_correction(detector, cb_centroid_list)
                else:
                    corr_cb_centroid_list = cb_centroid_list 
                if show_centroids:
                    print ('***** Measured centroids for checkbox sizes 3, 5, and 7, respectively:')
                    print ('      cb_centroid_list = ', corr_cb_centroid_list)
                    print ('           True center = ', true_center)
                # Show the display with the measured and true positions
                fig_name = os.path.join("../resultsXrandomstars", "centroid_displays/Star"+repr(st)+"_Scene"+repr(scene)+bg_choice+pos+".jpg")
                tf.display_centroids(detector, st, case, psf, true_center32x32, cb_centroid_list, 
                                     show_disp, vlim, savefile=save_centroid_disp, fig_name=fig_name)  
                # Find the best checkbox size = minimum difference with true values
                min_diff, _ = tf.get_mindiff(differences_true_TA[0][0], differences_true_TA[0][1], differences_true_TA[0][2])
                # Write output into text file
                if save_text_file:
                    position = "_Position1"
                    if pos == "_Position2":
                        position = "_Position2"
                    output_file = os.path.join(output_file_path, "centroids_Scene"+repr(scene)+bg_choice+position+".txt")
                    data2write = [save_text_file, output_file, st, background2use, corr_cb_centroid_list, true_center, loleftcoords, mag_i, min_diff]
                    writedatafile(show_centroids, data2write, lines4screenandfile) 
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
                                                                  "Checkbox: 3", "5", "7", 
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
    line1 = "{:<5} {:<10} {:<14} {:<16} {:<14} {:<16} {:<14} {:<16} {:<14} {:<16} {:<8} {:<12} {:<10.2}".format(
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
    resultsTEST1 = runTEST(test2perform, detectors, transf_direction, case, stars_sample, P1P2data, bench_starP1, trueVsP1, trueVsP2)
    T1P1P2data, T1_transformations, T1_diffs, T1_benchVs_list = resultsTEST1
    T1_V2_3, T1_V3_3, T1_V2_5, T1_V3_5, T1_V2_7, T1_V3_7 = T1_transformations
    T1_diffV2_3, T1_diffV3_3, T1_diffV2_5, T1_diffV3_5, T1_diffV2_7, T1_diffV3_7 = T1_diffs
    T1bench_V2_list, T1bench_V3_list = T1_benchVs_list
    # Get the statistics
    results_stats = get_stats(case, T1_transformations, T1_diffs, T1_benchVs_list, Nsigma, max_iters_Nsig)
    # unfold results
    T1_st_devsAndMeans, T1_diff_counter, T1_bench_values, T1_sigmas_deltas, T1_sigma_reject, rejected_elementsLS, rejected_eleNsig = results_stats
    T1stdev_V2_3, T1mean_V2_3, T1stdev_V2_5, T1mean_V2_5, T1stdev_V2_7, T1mean_V2_7, T1stdev_V3_3, T1mean_V3_3, T1stdev_V3_5, T1mean_V3_5, T1stdev_V3_7, T1mean_V3_7 = T1_st_devsAndMeans
    T1_min_diff, T1_counter = T1_diff_counter
    T1bench_V2, T1bench_V3 = T1_bench_values
    T1LSdeltas_3, T1LSsigmas_3, T1LSlines2print_3, T1LSdeltas_5, T1LSsigmas_5, T1LSlines2print_5, T1LSdeltas_7, T1LSsigmas_7, T1LSlines2print_7 = T1_sigmas_deltas
    T1sigmaV2_3, T1meanV2_3, T1sigmaV3_3, T1meanV3_3, T1newV2_3, T1newV3_3, T1niter_3, T1lines2print_3, T1sigmaV2_5, T1meanV2_5, T1sigmaV3_5, T1meanV3_5, T1newV2_5, T1newV3_5, T1niter_5, T1lines2print_5, T1sigmaV2_7, T1meanV2_7, T1sigmaV3_7, T1meanV3_7, T1newV2_7, T1newV3_7, T1niter_7, T1lines2print_7 = T1_sigma_reject

# TEST 2: (a) Transform individual P1 and P2 to V2-V3, (b) avg V2-V3 space positions, (c) compare to avg reference positions
if test2perform == "T2":
    resultsTEST2 = runTEST(test2perform, detectors, transf_direction, case, stars_sample, P1P2data, bench_starP1, trueVsP1, trueVsP2)
    T2P1P2data, T2_transformations, T2_diffs, T2_benchVs_list = resultsTEST2
    T2_V2_3, T2_V3_3, T2_V2_5, T2_V3_5, T2_V2_7, T2_V3_7 = T2_transformations
    T2_diffV2_3, T2_diffV3_3, T2_diffV2_5, T2_diffV3_5, T2_diffV2_7, T2_diffV3_7 = T2_diffs
    T2bench_V2_list, T2bench_V3_list = T2_benchVs_list
    # Get the statistics
    results_stats = get_stats(case, T2_transformations, T2_diffs, T2_benchVs_list, Nsigma, max_iters_Nsig)
    # unfold results
    T2_st_devsAndMeans, T2_diff_counter, T2_bench_values, T2_sigmas_deltas, T2_sigma_reject, rejected_elementsLS, rejected_eleNsig = results_stats
    T2stdev_V2_3, T2mean_V2_3, T2stdev_V2_5, T2mean_V2_5, T2stdev_V2_7, T2mean_V2_7, T2stdev_V3_3, T2mean_V3_3, T2stdev_V3_5, T2mean_V3_5, T2stdev_V3_7, T2mean_V3_7 = T2_st_devsAndMeans
    T2_min_diff, T2_counter = T2_diff_counter
    T2bench_V2, T2bench_V3 = T2_bench_values
    T2LSdeltas_3, T2LSsigmas_3, T2LSlines2print_3, T2LSdeltas_5, T2LSsigmas_5, T2LSlines2print_5, T2LSdeltas_7, T2LSsigmas_7, T2LSlines2print_7 = T2_sigmas_deltas
    T2sigmaV2_3, T2meanV2_3, T2sigmaV3_3, T2meanV3_3, T2newV2_3, T2newV3_3, T2niter_3, T2lines2print_3, T2sigmaV2_5, T2meanV2_5, T2sigmaV3_5, T2meanV3_5, T2newV2_5, T2newV3_5, T2niter_5, T2lines2print_5, T2sigmaV2_7, T2meanV2_7, T2sigmaV3_7, T2meanV3_7, T2newV2_7, T2newV3_7, T2niter_7, T2lines2print_7 = T2_sigma_reject
    
# TEST 3: (a) Transform P1 and P2 individually to V2-V3 (b) compare star by star and position by position
if test2perform == "T3":
    resultsTEST3 = runTEST(test2perform, detectors, transf_direction, case, stars_sample, P1P2data, bench_starP1, trueVsP1, trueVsP2)
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
    T3_V2_3 = combine2arrays(T3_V2_13, T3_V2_23, T3_V2_3)
    T3_V2_5 = combine2arrays(T3_V2_15, T3_V2_25, T3_V2_5)
    T3_V2_7 = combine2arrays(T3_V2_17, T3_V2_27, T3_V2_7)
    T3_V3_3, T3_V3_5, T3_V3_7 = np.array([]), np.array([]), np.array([])
    T3_V3_3 = combine2arrays(T3_V3_13, T3_V3_23, T3_V3_3)
    T3_V3_5 = combine2arrays(T3_V3_15, T3_V3_25, T3_V3_5)
    T3_V3_7 = combine2arrays(T3_V3_17, T3_V3_27, T3_V3_7)
    T3_diffV2_3, T3_diffV2_5, T3_diffV2_7 = np.array([]), np.array([]), np.array([])
    T3_diffV2_3 = combine2arrays(T3_diffV2_13, T3_diffV2_23, T3_diffV2_3)
    T3_diffV2_5 = combine2arrays(T3_diffV2_15, T3_diffV2_25, T3_diffV2_5)
    T3_diffV2_7 = combine2arrays(T3_diffV2_17, T3_diffV2_27, T3_diffV2_7)
    T3_diffV3_3, T3_diffV3_5, T3_diffV3_7 = np.array([]), np.array([]), np.array([])
    T3_diffV3_3 = combine2arrays(T3_diffV3_13, T3_diffV3_23, T3_diffV3_3)
    T3_diffV3_5 = combine2arrays(T3_diffV3_15, T3_diffV3_25, T3_diffV3_5)
    T3_diffV3_7 = combine2arrays(T3_diffV3_17, T3_diffV3_27, T3_diffV3_7)
    T3bench_V2_list, T3bench_V3_list = np.array([]), np.array([])
    T3bench_V2_list = combine2arrays(np.array(T3bench_V2_listP1), np.array(T3bench_V2_listP2), T3bench_V2_list)
    T3bench_V3_list = combine2arrays(np.array(T3bench_V3_listP1), np.array(T3bench_V3_listP2), T3bench_V3_list)
    T3bench_V2_list.tolist()
    T3bench_V3_list.tolist()
    # Get the statistics
    T3_transformations = [T3_V2_3, T3_V3_3, T3_V2_5, T3_V3_5, T3_V2_7, T3_V3_7]
    T3_diffs = [T3_diffV2_3, T3_diffV3_3, T3_diffV2_5, T3_diffV3_5, T3_diffV2_7, T3_diffV3_7]
    T3_benchVs_list = [T3bench_V2_list, T3bench_V3_list] 
    results_stats = get_stats(case, T3_transformations, T3_diffs, T3_benchVs_list, Nsigma, max_iters_Nsig)
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

print_results(stars_sample, case, test2perform, diffs_in_arcsecs, Tstdev_Vs, Tmean_Vs, T_diff_counter, 
              save_text_file, TLSlines2print, Tlines2print, Tbench_Vs_list, T_Vs, T_diffVs, 
              rejected_elementsLS, rejected_eleNsig)

    
print ("\n Script 'testXrandom_stars.py' finished! Took  %s  seconds to finish. \n" % ((time.time() - start_time)) )
