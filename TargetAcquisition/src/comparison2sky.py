from __future__ import print_function, division
from glob import glob
import numpy as np
import os
import collections
# other code
import coords_transform as ct
import testing_functions as tf 
print("Modules correctly imported! \n")

"""
This script tests which of the following options obtains the smaller difference to the true sky positions (V2-V3 space):
1. Average positions P1 and P2, transform to V2-V3 space, and compare to average reference positions (V2-V3 space)
2. Transform individual positions P1 and P2 to V2-V3 space, average V2-V3 space positions, and compare to 
   average reference positions.
3. Transform P1 and P2 individually to V2-V3 space and compare star by star and position by position.

It outputs 3 text file with results of test per case into directory TargetAcquisition/Coords_transforms/results/
* case depends on scene, noise, and shutter velocity; we have 24 cases 
"""

#######################################################################################################################

# general settings
detector = 491           # detector, integer: for now only 491 available
bkgd_method = "all"      # background to test, string: all, None, fixed, frac  
filter_input = "F140X"   # Filter, string: for now only test case is F140X
show_positions = True    # Print positions on file and screen: True or False
tilt = False             # tilt angle: True or False
debug = False            # See screen print statements for intermediate answers: True or False 
save_txt_file = False    # Save text file with resulting transformations: True or False

#######################################################################################################################

#  --> FUNCTIONS

def convert2fulldetector(stars, P1P2data, bench_star, benchmark_xLyL):
    """ This function simply converts from 32x32 pixel to full detector coordinates according to 
    background method - lengths are different for the fractional case. """
    x13,y13, x23,y23, x15,y15, x25,y25, x17,y17, x27,y27 = P1P2data
    benchxL, benchyL = benchmark_xLyL
    if len(stars.tolist()) == len(bench_star.tolist()):   # for the fixed and None background case
        x13, y13 = x13+benchxL, y13+benchyL
        x23, y23 = x23+benchxL, y23+benchyL
        x15, y15 = x15+benchxL, y15+benchyL
        x25, y25 = x25+benchxL, y25+benchyL
        x17, y17 = x17+benchxL, y17+benchyL
        x27, y27 = x27+benchxL, y27+benchyL
    else:                               # for the fractional background case
        for i, s in enumerate(stars):
            if s in bench_star:
                j = bench_star.tolist().index(s)
                x13[i], y13[i] = x13[i]+benchxL[j], y13[i]+benchyL[j]
                x23[i], y23[i] = x23[i]+benchxL[j], y23[i]+benchyL[j]
                x15[i], y15[i] = x15[i]+benchxL[j], y15[i]+benchyL[j]
                x25[i], y25[i] = x25[i]+benchxL[j], y25[i]+benchyL[j]
                x17[i], y17[i] = x17[i]+benchxL[j], y17[i]+benchyL[j]
                x27[i], y27[i] = x27[i]+benchxL[j], y27[i]+benchyL[j]
    return x13,y13, x23,y23, x15,y15, x25,y25, x17,y17, x27,y27
    
    
def read_star_param_file(case, path4starfiles, paths_list):
    """ This function reads the corresponding star parameters file and returns the data. """
    if "Scene1_slow_real" in case:
        dir2test = path4starfiles+paths_list[0]
    elif "Scene1_slow_nonoise" in case:
        dir2test = path4starfiles+paths_list[1]
    elif "Scene1_rapid_real" in case:
        dir2test = path4starfiles+paths_list[2]
    elif "Scene1_rapid_nonoise" in case:
        dir2test = path4starfiles+paths_list[3]
    elif "Scene1_slow_real_shifted" in case:
        dir2test = path4starfiles+paths_list[4]
    elif "Scene1_slow_nonoise_shifted" in case:
        dir2test = path4starfiles+paths_list[5]
    elif "Scene1_rapid_real_shifted" in case:
        dir2test = path4starfiles+paths_list[6]
    elif "Scene1_rapid_nonoise_shifted" in case:
        dir2test = path4starfiles+paths_list[7]
    elif "Scene2_slow_real" in case:
        dir2test = path4starfiles+paths_list[8]
    elif "Scene2_slow_nonoise" in case:
        dir2test = path4starfiles+paths_list[9]
    elif "Scene2_rapid_real" in case:
        dir2test = path4starfiles+paths_list[10]
    elif "Scene2_rapid_nonoise" in case:
        dir2test = path4starfiles+paths_list[11]
    elif "Scene2_slow_real_shifted" in case:
        dir2test = path4starfiles+paths_list[12]
    elif "Scene2_slow_nonoise_shifted" in case:
        dir2test = path4starfiles+paths_list[13]
    elif "Scene2_rapid_real_shifted" in case:
        dir2test = path4starfiles+paths_list[14]
    elif "Scene2_rapid_nonoise_shifted" in case:
        dir2test = path4starfiles+paths_list[15]
    # NOTE: These positions are all 0-indexed!
    #    xL:  x-coordinate of the left edge of the postge stamp in the full image (range 0-2047)
    #    xR: x-coord of right edge of the postage stamp
    #    yL: y-coord of the lower edge of the postage stamp
    #    yU:  y-coord of the upper edge of the postage stamp
    star_param_txt = os.path.join(dir2test,"star parameters.txt")
    benchmark_data = np.loadtxt(star_param_txt, skiprows=3, unpack=True)
    bench_star, quadrant, star_in_quad, x_491, y_491, x_492, y_492, V2, V3, xL, xR, yL, yU = benchmark_data
    return bench_star, quadrant, star_in_quad, x_491, y_491, x_492, y_492, V2, V3, xL, xR, yL, yU
    
    
def compare2ref(case, path4starfiles, paths_list, benchV2, benchV3, stars, V2in, V3in):
    """ This function obtains the differences of the input arrays with the reference or benchmark data. """
    # calculate the differences with respect to the benchmark data
    if len(stars.tolist()) == len(bench_star.tolist()):   # for the fixed and None background case
        diffV2 = benchV2 - V2in
        diffV3 = benchV3 - V2in
        bench_V2_list = benchV2.tolist()
        bench_V3_list = benchV3.tolist()
    else:                               # for the fractional background case
        bench_V2_list, bench_V3_list = [], []
        diffV2, diffV3 = [], []
        for i, s in enumerate(stars):
            if s in bench_star:
                j = bench_star.tolist().index(s)
                dsV2 = benchV2[j] - V2in[i] 
                dsV3 = benchV3[j] - V3in[i] 
                diffV2.append(dsV2)
                diffV3.append(dsV3)
                bench_V2_list.append(benchV2[j])
                bench_V3_list.append(benchV3[j])
        diffV2 = np.array(diffV2)
        diffV3 = np.array(diffV3)
    return diffV2, diffV3, bench_V2_list, bench_V3_list


def get_mindiff(d1, d2, d3):
    """ This function determines the minimum difference from checkboxes 3, 5, and 7,
    and counts the number of repetitions. """
    min_diff = []
    for i, _ in enumerate(d1):
        diffs_list = [d1[i], d2[i], d3[i]]
        md = min(diffs_list)
        if md == d1[i]:
            m_diff = 3
        elif md == d2[i]:
            m_diff = 5
        elif md == d3[i]:
            m_diff = 7
        min_diff.append(m_diff)
    counter=collections.Counter(min_diff)
    return min_diff, counter


#######################################################################################################################

#  --> CODE

single_case = None       # test only a particular case: integer number of index from test_file, else set to None 
if single_case != None:
    save_txt_file = False   # just in case, do not overwrite text file with only one case.

# Define the paths for results and inputs
path4results = "../Coords_transforms/results/"
path4inputP1P2 = "../PFforMaria/comparison_txt_positions/"

# Get true positions from Pierre's position files from Tony's star parameters file 
# Paths to Scenes 1 and 2 local directories: /Users/pena/Documents/AptanaStudio3/NIRSpec/TargetAcquisition/
path4starfiles = "../PFforMaria/"
path_scene1_slow = "Scene_1_AB23/NIRSpec_TA_Sim_AB23 first NRS/postage_redo"
path_scene1_slow_nonoise = "Scene_1_AB23/NIRSpec_TA_Sim_AB23 first NRS no_noise/postage_redo"
path_scene1_rapid = "Scene_1_AB23/NIRSpec_TA_Sim_AB23 first NRSRAPID/postage_redo"
path_scene1_rapid_nonoise = "Scene_1_AB23/NIRSpec_TA_Sim_AB23 first NRS no_noise/postage_redo"
path_scene1_slow_shifted = "Scene_1_AB23/NIRSpec_TA_Sim_AB23 shifted NRS/postage_redo"
path_scene1_slow_shifted_nonoise = "Scene_1_AB23/NIRSpec_TA_Sim_AB23 shifted NRS no_noise/postage_redo"
path_scene1_rapid_shifted = "Scene_1_AB23/NIRSpec_TA_Sim_AB23 shifted NRSRAPID/postage_redo"
path_scene1_rapid_shifted_nonoise = "Scene_1_AB23/NIRSpec_TA_Sim_AB23 shifted NRS no_noise/postage_redo"
path_scene2_slow = "Scene_2_AB1823/NIRSpec_TA_Sim_AB1823 first NRS/postage_redo"
path_scene2_slow_nonoise = "Scene_2_AB1823/NIRSpec_TA_Sim_AB1823 first NRS no_noise/postage_redo"
path_scene2_rapid = "Scene_2_AB1823/NIRSpec_TA_Sim_AB1823 first NRSRAPID/postage_redo"
path_scene2_rapid_nonoise = "Scene_2_AB1823/NIRSpec_TA_Sim_AB1823 first NRSRAPID no_noise/postage_redo"
path_scene2_slow_shifted = "Scene_2_AB1823/NIRSpec_TA_Sim_AB1823 shifted NRS/postage_redo"
path_scene2_slow_shifted_nonoise = "Scene_2_AB1823/NIRSpec_TA_Sim_AB1823 shifted NRS no_noise/postage_redo"
path_scene2_rapid_shifted = "Scene_2_AB1823/NIRSpec_TA_Sim_AB1823 shifted NRSRAPID/postage_redo"
path_scene2_rapid_shifted_nonoise = "Scene_2_AB1823/NIRSpec_TA_Sim_AB1823 shifted NRSRAPID no_noise/postage_redo"

paths_list = [path_scene1_slow, path_scene1_slow_nonoise, 
              path_scene1_rapid, path_scene1_rapid_nonoise,
              path_scene1_slow_shifted, path_scene1_slow_shifted_nonoise, 
              path_scene1_rapid_shifted, path_scene1_rapid_shifted_nonoise, 
              path_scene2_slow, path_scene2_slow_nonoise, 
              path_scene2_rapid, path_scene2_rapid_nonoise,
              path_scene2_slow_shifted, path_scene2_slow_shifted_nonoise, 
              path_scene2_rapid_shifted, path_scene2_rapid_shifted_nonoise]

# Get list of all input files of the background method given to test
if bkgd_method == "all":
    input_files_list = glob(path4inputP1P2+"*.txt")
else:
    input_files_list = glob(path4inputP1P2+"*"+bkgd_method+"*.txt")

# Start TESTS the loop
for infile in input_files_list:
    # define the case to work with 
    case = os.path.basename(infile).replace(".txt", "")
    case = case.replace("_centroids", "")
    
    # determine the background method used
    if "None" in case:
        bg_method = "_None"
    elif "fixed" in case:
        bg_method = "_fixed"
    elif "frac" in case:
        bg_method = "_frac"
    print ("* For background method: ", bg_method)
    
    # get the benchmark data
    bench_star, _, _, _, _, _, _, benchV2, benchV3, benchxL, _, benchyL, _ = read_star_param_file(case, path4starfiles, paths_list)

    # read the measured detector centroids
    data = np.loadtxt(infile, skiprows=2, usecols=(0,1,2,3,4,5,6,7,8,9,10,11,12,13), unpack=True)
    stars, bg_value, x13, y13, x15, y15, x17, y17, x23, y23, x25, y25, x27, y27 = data
    
    # convert from 32x32 pixel to full detector coordinates
    P1P2data = [x13,y13, x23,y23, x15,y15, x25,y25, x17,y17, x27,y27]
    benchmark_xLyL = [benchxL, benchyL]
    x13,y13, x23,y23, x15,y15, x25,y25, x17,y17, x27,y27 = convert2fulldetector(stars, P1P2data, bench_star, benchmark_xLyL)
    
    # TEST 1: (a) Avg P1 and P2, (b) transform to V2-V3, (c) compare to avg reference positions (V2-V3 space)
    transf_direction = "forward"
    # Step (a) - averages
    avgx3 = (x13+x23)/2.0
    avgy3 = (y13+y23)/2.0
    avgx5 = (x15+x25)/2.0
    avgy5 = (y15+y25)/2.0
    avgx7 = (x17+x27)/2.0
    avgy7 = (y17+y27)/2.0
    # Step (b) - transformations
    T1_V2_3, T1_V3_3 = ct.coords_transf(transf_direction, detector, filter_input, avgx3, avgy3, tilt, debug)
    T1_V2_5, T1_V3_5 = ct.coords_transf(transf_direction, detector, filter_input, avgx5, avgy5, tilt, debug)
    T1_V2_7, T1_V3_7 = ct.coords_transf(transf_direction, detector, filter_input, avgx7, avgy7, tilt, debug)
    # Step (c) - comparison
    T1_diffV2_3, T1_diffV3_3, bench_V2_list, bench_V3_list = compare2ref(case, path4starfiles, paths_list, benchV2, benchV3, stars, T1_V2_3, T1_V3_3)
    T1_diffV2_5, T1_diffV3_5, _, _ = compare2ref(case, path4starfiles, paths_list, benchV2, benchV3, stars, T1_V2_5, T1_V3_5)
    T1_diffV2_7, T1_diffV3_7, _, _ = compare2ref(case, path4starfiles, paths_list, benchV2, benchV3, stars, T1_V2_7, T1_V3_7)
    # get the minimum of the differences
    T1_min_diff, T1_counter = get_mindiff(T1_diffV2_3, T1_diffV2_5, T1_diffV2_7)
    # calculate standard deviations and means
    T1stdev_V2_3, T1mean_V2_3 = tf.find_std(T1_diffV2_3)
    T1stdev_V2_5, T1mean_V2_5 = tf.find_std(T1_diffV2_5)
    T1stdev_V2_7, T1mean_V2_7 = tf.find_std(T1_diffV2_7)
    T1stdev_V3_3, T1mean_V3_3 = tf.find_std(T1_diffV3_3)
    T1stdev_V3_5, T1mean_V3_5 = tf.find_std(T1_diffV3_5)
    T1stdev_V3_7, T1mean_V3_7 = tf.find_std(T1_diffV3_7)
    
    # TEST 2: (a) Transform individual P1 and P2 to V2-V3, (b) avg V2-V3 space positions, (c) compare to avg reference positions
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
    T2_diffV2_3, T2_diffV3_3, _, _ = compare2ref(case, path4starfiles, paths_list, benchV2, benchV3, stars, T2_V2_3, T2_V3_3)
    T2_diffV2_5, T2_diffV3_5, _, _ = compare2ref(case, path4starfiles, paths_list, benchV2, benchV3, stars, T2_V2_5, T2_V3_5)
    T2_diffV2_7, T2_diffV3_7, _, _ = compare2ref(case, path4starfiles, paths_list, benchV2, benchV3, stars, T2_V2_7, T2_V3_7)
    # get the minimum of the differences
    T2_min_diff, T2_counter = get_mindiff(T2_diffV2_3, T2_diffV2_5, T2_diffV2_7)
    # calculate standard deviations and means
    T2stdev_V2_3, T2mean_V2_3 = tf.find_std(T2_diffV2_3)
    T2stdev_V2_5, T2mean_V2_5 = tf.find_std(T2_diffV2_5)
    T2stdev_V2_7, T2mean_V2_7 = tf.find_std(T2_diffV2_7)
    T2stdev_V3_3, T2mean_V3_3 = tf.find_std(T2_diffV3_3)
    T2stdev_V3_5, T2mean_V3_5 = tf.find_std(T2_diffV3_5)
    T2stdev_V3_7, T2mean_V3_7 = tf.find_std(T2_diffV3_7)
    
    # TEST 3: (a) Transform P1 and P2 individually to V2-V3 (b) compare star by star and position by position
    # Step (a) - transformations
    T3_V2_13, T3_V3_13 = ct.coords_transf(transf_direction, detector, filter_input, x13, y13, tilt, debug)
    T3_V2_15, T3_V3_15 = ct.coords_transf(transf_direction, detector, filter_input, x15, y15, tilt, debug)
    T3_V2_17, T3_V3_17 = ct.coords_transf(transf_direction, detector, filter_input, x17, y17, tilt, debug)
    T3_V2_23, T3_V3_23 = ct.coords_transf(transf_direction, detector, filter_input, x23, y23, tilt, debug)
    T3_V2_25, T3_V3_25 = ct.coords_transf(transf_direction, detector, filter_input, x25, y25, tilt, debug)
    T3_V2_27, T3_V3_27 = ct.coords_transf(transf_direction, detector, filter_input, x27, y27, tilt, debug)
    # Step (b) - comparison
    T3_diffV2_13, T3_diffV3_13, _, _ = compare2ref(case, path4starfiles, paths_list, benchV2, benchV3, stars, T3_V2_13, T3_V3_13)
    T3_diffV2_23, T3_diffV3_23, _, _ = compare2ref(case, path4starfiles, paths_list, benchV2, benchV3, stars, T3_V2_23, T3_V3_23)
    T3_diffV2_15, T3_diffV3_15, _, _ = compare2ref(case, path4starfiles, paths_list, benchV2, benchV3, stars, T3_V2_15, T3_V3_15)
    T3_diffV2_25, T3_diffV3_25, _, _ = compare2ref(case, path4starfiles, paths_list, benchV2, benchV3, stars, T3_V2_25, T3_V3_25)
    T3_diffV2_17, T3_diffV3_17, _, _ = compare2ref(case, path4starfiles, paths_list, benchV2, benchV3, stars, T3_V2_17, T3_V3_17)
    T3_diffV2_27, T3_diffV3_27, _, _ = compare2ref(case, path4starfiles, paths_list, benchV2, benchV3, stars, T3_V2_27, T3_V3_27)
    # get the minimum of the differences
    T3_min_diff1, T3_counter1 = get_mindiff(T3_diffV2_13, T3_diffV2_15, T3_diffV2_17)
    T3_min_diff2, T3_counter2 = get_mindiff(T3_diffV2_23, T3_diffV2_25, T3_diffV2_27)
    # calculate standard deviations and means
    T3stdev_V2_13, T3mean_V2_13 = tf.find_std(T3_diffV2_13)
    T3stdev_V2_15, T3mean_V2_15 = tf.find_std(T3_diffV2_15)
    T3stdev_V2_17, T3mean_V2_17 = tf.find_std(T3_diffV2_17)
    T3stdev_V3_13, T3mean_V3_13 = tf.find_std(T3_diffV3_13)
    T3stdev_V3_15, T3mean_V3_15 = tf.find_std(T3_diffV3_15)
    T3stdev_V3_17, T3mean_V3_17 = tf.find_std(T3_diffV3_17)
    T3stdev_V2_23, T3mean_V2_23 = tf.find_std(T3_diffV2_23)
    T3stdev_V2_25, T3mean_V2_25 = tf.find_std(T3_diffV2_25)
    T3stdev_V2_27, T3mean_V2_27 = tf.find_std(T3_diffV2_27)
    T3stdev_V3_23, T3mean_V3_23 = tf.find_std(T3_diffV3_23)
    T3stdev_V3_25, T3mean_V3_25 = tf.find_std(T3_diffV3_25)
    T3stdev_V3_27, T3mean_V3_27 = tf.find_std(T3_diffV3_27)
    
    # Print results to screen and save into a text file if told so
    # Text file 1
    line0 = "{}".format("Differences = diffs = True_Positions - Measured_Positions")
    line1 = "{}".format("Test1: average P1 and P2, transform to V2-V3, calculate differences")
    line2a = "std_dev_V2_3 = {:<20}           std_dev_V3_3 = {:<20}".format(T1stdev_V2_3, T1stdev_V3_3)
    line2b = "std_dev_V2_5 = {:<20}           std_dev_V3_5 = {:<20}".format(T1stdev_V2_5, T1stdev_V3_5)
    line2c = "std_dev_V2_7 = {:<20}           std_dev_V3_7 = {:<20}".format(T1stdev_V2_7, T1stdev_V3_7)
    line3a = "   mean_V2_3 = {:<22}            mean_V3_3 = {:<22}".format(T1mean_V2_3, T1mean_V3_3)
    line3b = "   mean_V2_5 = {:<22}            mean_V3_5 = {:<22}".format(T1mean_V2_5, T1mean_V3_5)
    line3c = "   mean_V2_7 = {:<22}            mean_V3_7 = {:<22}".format(T1mean_V2_7, T1mean_V3_7)
    line3bisA = "Repetitions Diffs1: {}".format(T1_counter)
    if show_positions:
        line4 = "{:<5} {:<20} {:<40} {:<40} {:<35} {:<30} {:<40} {:<40} {:<24} {:<6}".format(
                        "Star", "BG_value", "Avg_Pos_Checkbox_3", "Avg_Pos_Checkbox_5", "Avg_PosCheckbox_7",
                        "True_Pos", "Diff_Chbx_3", "Diff_Chbx_5", "Diff_Chbx_7", "MinDiff")
        line5 = "{:>25} {:>17} {:>22} {:>17} {:>22} {:>17} {:>17} {:>11} {:>18} {:>17} {:>22} {:>17} {:>22} {:>17}".format(
                        "x", "y", "x", "y", "x", "y", "x", "y", "x", "y", "x", "y", "x", "y")
    else:
        line4 = "{:<5} {:<20} {:<40} {:<40} {:<23} {:<7}".format("Star", "BG_value", "Checkbox=3", "Checkbox=5", "Checkbox=7", "MinDiff")
        line5 = "{:>25} {:>17} {:>22} {:>17} {:>22} {:>17}".format("x", "y", "x", "y", "x", "y")        
    print (line0)
    print (line1)
    print (line2a)
    print (line2b)
    print (line2c)
    print (line3a)
    print (line3b)
    print (line3c)
    print (line3bisA)
    print (line4)
    print (line5)
    if save_txt_file:
        txt_out = path4results+"Test1_results"+bg_method+".txt"
        to = open(txt_out, "w+")
        to.write(line0+"\n")
        to.write(line1+"\n")
        to.write(line2a+"\n")
        to.write(line2b+"\n")
        to.write(line2c+"\n")
        to.write(line3a+"\n")
        to.write(line3b+"\n")
        to.write(line3c+"\n")
        to.write(line3bisA+"\n")
        to.write(line4+"\n")
        to.write(line5+"\n")
    for i, st in enumerate(stars):
        st = int(st)
        if show_positions:
            line6 = "{:<5} {:<5} {:>20}  {:<20} {:>18}  {:<20} {:>18}  {:<20} {:>10}  {:<14} {:>18}  {:<20} {:>18}  {:<20} {:>18}  {:<17} {:<7}".format(
                        st, bg_value[i], 
                        T1_V2_3[i], T1_V3_3[i], T1_V2_5[i], T1_V3_5[i], T1_V2_7[i], T1_V3_7[i], 
                        bench_V2_list[i], bench_V3_list[i],
                        T1_diffV2_3[i], T1_diffV3_3[i], T1_diffV2_5[i], T1_diffV3_5[i], T1_diffV2_7[i], T1_diffV3_7[i], T1_min_diff[i])
        else:
            line6 = "{:<5} {:<5} {:>20}  {:<20} {:>18}  {:<20} {:>18}  {:<17} {:<7}".format(st, bg_value[i], 
                    T1_diffV2_3[i], T1_diffV3_3[i], T1_diffV2_5[i], T1_diffV3_5[i], T1_diffV2_7[i], T1_diffV3_7[i], T1_min_diff[i])
        print (line6)
        if save_txt_file:
            to.write(line6+"\n")
    if save_txt_file:
        to.close()
        print (" * Results saved in file: ", txt_out)
    #raw_input(" * Press enter to continue... \n")

    # Text file 2
    line0 = "{}".format("Differences = True_Positions - Measured_Positions")
    line1 = "{}".format("Test2: P1 P2, average positions in V2-V3, calculate differences")
    line2a = "std_dev_V2_3 = {:<20}   std_dev_V3_3 = {:<20}".format(T2stdev_V2_3, T2stdev_V3_3)
    line2b = "std_dev_V2_5 = {:<20}   std_dev_V3_5 = {:<20}".format(T2stdev_V2_5, T2stdev_V3_5)
    line2c = "std_dev_V2_7 = {:<20}   std_dev_V3_7 = {:<20}".format(T2stdev_V2_7, T2stdev_V3_7)
    line3a = "   mean_V2_3 = {:<22}    mean_V3_3 = {:<22}".format(T2mean_V2_3, T2mean_V3_3)
    line3b = "   mean_V2_5 = {:<22}    mean_V3_5 = {:<22}".format(T2mean_V2_5, T2mean_V3_5)
    line3c = "   mean_V2_7 = {:<22}    mean_V3_7 = {:<22}".format(T2mean_V2_7, T2mean_V3_7)
    line3bisA = "Repetitions Diffs1: {}".format(T2_counter)
    if show_positions:
        line4 = "{:<5} {:<20} {:<40} {:<40} {:<35} {:<30} {:<40} {:<40} {:<23} {:<7}".format(
                        "Star", "BG_value", "Avg_Pos_Checkbox_3", "Avg_Pos_Checkbox_5", "Avg_PosCheckbox_7",
                        "True_Pos", "Diff_Chbx_3", "Diff_Chbx_5", "Diff_Chbx_7", "MinDiff")
        line5 = "{:>25} {:>17} {:>22} {:>17} {:>22} {:>17} {:>17} {:>11} {:>18} {:>17} {:>22} {:>17} {:>22} {:>17}".format(
                        "x", "y", "x", "y", "x", "y", "x", "y", "x", "y", "x", "y", "x", "y")
    else:
        line4 = "{:<5} {:<20} {:<40} {:<40} {:<23} {:<7}".format("Star", "BG_value", "Checkbox=3", "Checkbox=5", "Checkbox=7", "MinDiff")
        line5 = "{:>25} {:>17} {:>22} {:>17} {:>22} {:>17}".format("x", "y", "x", "y", "x", "y")        
    print (line0)
    print (line1)
    print (line2a)
    print (line2b)
    print (line2c)
    print (line3a)
    print (line3b)
    print (line3c)
    print (line3bisA)
    print (line4)
    print (line5)
    if save_txt_file:
        txt_out = path4results+"Test2_results"+bg_method+".txt"
        to = open(txt_out, "w+")
        to.write(line0+"\n")
        to.write(line1+"\n")
        to.write(line2a+"\n")
        to.write(line2b+"\n")
        to.write(line2c+"\n")
        to.write(line3a+"\n")
        to.write(line3b+"\n")
        to.write(line3c+"\n")
        to.write(line3bisA+"\n")
        to.write(line4+"\n")
        to.write(line5+"\n")
    for i, st in enumerate(stars):
        st = int(st)
        if show_positions:
            line6 = "{:<5} {:<5} {:>20}  {:<20} {:>18}  {:<20} {:>18}  {:<20} {:>10}  {:<14} {:>18}  {:<20} {:>18}  {:<20} {:>18}  {:<20} {:<7}".format(
                        st, bg_value[i], 
                        T2_V2_3[i], T2_V3_3[i], T2_V2_5[i], T2_V3_5[i], T2_V2_7[i], T2_V3_7[i], 
                        bench_V2_list[i], bench_V3_list[i],
                        T2_diffV2_3[i], T2_diffV3_3[i], T2_diffV2_5[i], T2_diffV3_5[i], T2_diffV2_7[i], T2_diffV3_7[i], T2_min_diff[i])
        else:
            line6 = "{:<5} {:<5} {:>20}  {:<20} {:>18}  {:<20} {:>18}  {:<20} {:<7}".format(st, bg_value[i], 
                    T2_diffV2_3[i], T2_diffV3_3[i], T2_diffV2_5[i], T2_diffV3_5[i], T2_diffV2_7[i], T2_diffV3_7[i], T2_min_diff[i])            
        print (line6)
        if save_txt_file:
            to.write(line6+"\n")
    if save_txt_file:
        to.close()
        print (" * Results saved in file: ", txt_out)
    #raw_input(" * Press enter to continue... \n")

    # Text file 3
    line1 = "{}".format("Test3: P1 and P2, transform to V2-V3 space individually, calculate differences position to position")
    line2a = "std_dev_V2_P1_3 = {:<20}           std_dev_V3_P1_3 = {:<20}".format(T3stdev_V2_13, T3stdev_V3_13)
    line2b = "std_dev_V2_P1_5 = {:<20}           std_dev_V3_P1_5 = {:<20}".format(T3stdev_V2_15, T3stdev_V3_15)
    line2c = "std_dev_V2_P1_7 = {:<20}           std_dev_V3_P1_7 = {:<20}".format(T3stdev_V2_17, T3stdev_V3_17)
    line3a = "   mean_V2_P1_3 = {:<22}            mean_V3_P1_3 = {:<22}".format(T3mean_V2_13, T3mean_V3_13)
    line3b = "   mean_V2_P1_5 = {:<22}            mean_V3_P1_5 = {:<22}".format(T3mean_V2_15, T3mean_V3_15)
    line3c = "   mean_V2_P1_7 = {:<22}            mean_V3_P1_7 = {:<22}".format(T3mean_V2_17, T3mean_V3_17)
    line4a = "std_dev_V2_P2_3 = {:<20}           std_dev_V3_P2_3 = {:<20}".format(T3stdev_V2_23, T3stdev_V3_23)
    line4b = "std_dev_V2_P2_5 = {:<20}           std_dev_V3_P2_5 = {:<20}".format(T3stdev_V2_25, T3stdev_V3_25)
    line4c = "std_dev_V2_P2_7 = {:<20}           std_dev_V3_P2_7 = {:<20}".format(T3stdev_V2_27, T3stdev_V3_27)
    line5a = "   mean_V2_P2_3 = {:<22}            mean_V3_P2_3 = {:<22}".format(T3mean_V2_23, T3mean_V3_23)
    line5b = "   mean_V2_P2_5 = {:<22}            mean_V3_P2_5 = {:<22}".format(T3mean_V2_25, T3mean_V3_25)
    line5c = "   mean_V2_P2_7 = {:<22}            mean_V3_P2_7 = {:<22}".format(T3mean_V2_27, T3mean_V3_27)
    line5bisA = "Repetitions Diffs1: {}".format(T3_counter1)
    line5bisB = "Repetitions Diffs2: {}".format(T3_counter2)
    if show_positions:
        line6 = "{:<5} {:<15} {:<35} {:<39} {:<37} {:<38} {:<36} {:<38} {:<30} {:<40} {:<40} {:<40} {:<36} {:<38} {:<27} {:<4} {:<4}".format(
                    "Star", "BG_value", "Pos1_Checkbox_3", "Pos1_Checkbox_5", "Pos1_Checkbox_7",
                    "Pos2_Checkbox_3", "Pos2_Checkbox_5", "Pos2_Checkbox_7", 
                    "True_Pos", "Diff1_Chbx_3", "Diff1_Chbx_5", "Diff1_Chbx_7", 
                    "Diff2_Chbx_3", "Diff2_Chbx_5", "Diff2_Chbx_7", "MinDiff1", "MinDiff2")
        line7 = "{:>22} {:>14} {:>22} {:>14} {:>22} {:>14} {:>22} {:>14} {:>22} {:>14} {:>22} {:>14} {:>22} {:>10} {:>22} {:>14} {:>22} {:>14} {:>22} {:>14} {:>28} {:>14} {:>22} {:>14} {:>22} {:>14}".format(
                        "x", "y", "x", "y", "x", "y", "x", "y", "x", "y", "x", "y", "x", "y",
                        "x", "y", "x", "y", "x", "y", "x", "y", "x", "y", "x", "y", "x", "y", 
                        "x", "y", "x", "y", "x", "y")
    else:
        line6 = "{:<5} {:<20} {:<40} {:<40} {:<40} {:<40} {:<40} {:<25} {:<9} {:<10}".format(
                    "Star", "BG_value", "Diff1_Chbx_3", "Diff1_Chbx_5", "Diff1_Chbx_7",
                    "Diff2_Chbx_3", "Diff2_Chbx_5", "Diff2_Chbx_7", "MinDiff1", "MinDiff2")
        line7 = "{:>25} {:>17} {:>22} {:>17} {:>22} {:>17} {:>22} {:>17} {:>22} {:>17} {:>22} {:>17}".format("x", "y", "x", "y", "x", "y", "x", "y", "x", "y", "x", "y")        
    print (line0)
    print (line1)
    print (line2a)
    print (line2b)
    print (line2c)
    print (line3a)
    print (line3b)
    print (line3c)
    print (line4a)
    print (line4b)
    print (line4b)
    print (line5a)
    print (line5b)
    print (line5c)
    print (line5bisA)
    print (line5bisB)
    print (line6)
    print (line7)
    if save_txt_file:
        txt_out = path4results+"Test3_results"+bg_method+".txt"
        to = open(txt_out, "w+")
        to.write(line0+"\n")
        to.write(line1+"\n")
        to.write(line2a+"\n")
        to.write(line2b+"\n")
        to.write(line2c+"\n")
        to.write(line3a+"\n")
        to.write(line3b+"\n")
        to.write(line3c+"\n")
        to.write(line4a+"\n")
        to.write(line4b+"\n")
        to.write(line4c+"\n")
        to.write(line5a+"\n")
        to.write(line5b+"\n")
        to.write(line5c+"\n")
        to.write(line5bisA+"\n")
        to.write(line5bisB+"\n")
        to.write(line6+"\n")
        to.write(line7+"\n")
    for i, st in enumerate(stars):
        st = int(st)
        if show_positions:
            line8 = "{:<5} {:<5} {:>16}  {:<19} {:>16}  {:<19} {:>16}  {:<19} {:>16}  {:<19} {:>16}  {:<19} {:>16}  {:<19} {:>14}  {:<14} {:>18}  {:<19} {:>18}  {:<19} {:>18}  {:<19} {:>18}  {:<19} {:>18}  {:<19} {:>18}  {:<19} {:<7} {:<7}".format(
                        st, bg_value[i], 
                        T3_V2_13[i], T3_V3_13[i], T3_V2_15[i], T3_V3_15[i], T3_V2_17[i], T3_V3_17[i], 
                        T3_V2_23[i], T3_V3_23[i], T3_V2_25[i], T3_V3_25[i], T3_V2_27[i], T3_V3_27[i], 
                        bench_V2_list[i], bench_V3_list[i],
                        T3_diffV2_13[i], T3_diffV3_13[i], T3_diffV2_15[i], T3_diffV3_15[i], T3_diffV2_17[i], T3_diffV3_17[i],
                        T3_diffV2_23[i], T3_diffV3_23[i], T3_diffV2_25[i], T3_diffV3_25[i], T3_diffV2_27[i], T3_diffV3_27[i],
                        T3_min_diff1[i], T3_min_diff2[i])
        else:
            line8 = "{:<5} {:<5} {:>20}  {:<20} {:>18}  {:<20} {:>18}  {:<20} {:>18}  {:<20} {:>18}  {:<20} {:>18}  {:<20} {:<7} {:<7}".format(st, bg_value[i], 
                    T3_diffV2_13[i], T3_diffV3_13[i], T3_diffV2_15[i], T3_diffV3_15[i], T3_diffV2_17[i], T3_diffV3_17[i], 
                    T3_diffV2_23[i], T3_diffV3_23[i], T3_diffV2_25[i], T3_diffV3_25[i], T3_diffV2_27[i], T3_diffV3_27[i], 
                    T3_min_diff1[i], T3_min_diff2[i])            
        print (line8)
        if save_txt_file:
            to.write(line8+"\n")
    if save_txt_file:
        to.close()
        print (" * Results saved in file: ", txt_out)
    if single_case:
        exit()
    else:
        raw_input(" * Press enter to continue... \n")
    
print ("\n Script 'comparison2sky.py' finished! ")
