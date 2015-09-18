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
* case depends on scene, noise, background value, and shutter velocity; results in 36 files per scene.
"""

#######################################################################################################################

# general settings
detector = 491           # detector, integer: for now only 491 available
bkgd_method = "frac"     # background to test, string: all, None, fixed, frac  
filter_input = "F140X"   # Filter, string: for now only test case is F140X
Pier_corr = True         # Include Pier's corrections to measured positions
show_positions = False   # Print positions on file and screen: True or False
tilt = False             # tilt angle: True or False
debug = False            # See screen print statements for intermediate answers: True or False 
save_txt_file = False    # Save text file with resulting transformations: True or False
diffs_in_arcsecs = True  # Print the differences in arcsecs? True or False (=degrees) 
single_case = 101        # test only a particular case: integer number of star, else set to None 
# Known bad stars in X and Y: 103, 105, 106, 112, 134, 152, 156, 170, 188

#######################################################################################################################

#  --> FUNCTIONS

def convert2fulldetector(detector, stars, P1P2data, bench_stars, benchmark_xLyL_P1, benchmark_xLyL_P2, Pier_corr=True):
    """ This function simply converts from 32x32 pixel to full detector coordinates according to 
    background method - lengths are different for the fractional case. """
    x13,y13, x23,y23, x15,y15, x25,y25, x17,y17, x27,y27 = P1P2data
    benchxL_P1, benchyL_P1 = benchmark_xLyL_P1
    benchxL_P2, benchyL_P2 = benchmark_xLyL_P2
    if len(stars) == len(bench_stars):   # for the fixed and None background case
        #print ("x13[0]+benchxL_P1=", x13[0], benchxL_P1[0], x13[0]+benchxL_P1)
        #raw_input()
        x13, y13 = x13+benchxL_P1, y13+benchyL_P1
        x23, y23 = x23+benchxL_P2, y23+benchyL_P2
        x15, y15 = x15+benchxL_P1, y15+benchyL_P1
        x25, y25 = x25+benchxL_P2, y25+benchyL_P2
        x17, y17 = x17+benchxL_P1, y17+benchyL_P1
        x27, y27 = x27+benchxL_P2, y27+benchyL_P2
    else:                               # for the fractional background case
        for i, s in enumerate(stars):
            if s in bench_stars:
                j = bench_stars.tolist().index(s)
                x13[i], y13[i] = x13[i]+benchxL_P1[j], y13[i]+benchyL_P1[j]
                x23[i], y23[i] = x23[i]+benchxL_P2[j], y23[i]+benchyL_P2[j]
                x15[i], y15[i] = x15[i]+benchxL_P1[j], y15[i]+benchyL_P1[j]
                x25[i], y25[i] = x25[i]+benchxL_P2[j], y25[i]+benchyL_P2[j]
                x17[i], y17[i] = x17[i]+benchxL_P1[j], y17[i]+benchyL_P1[j]
                x27[i], y27[i] = x27[i]+benchxL_P2[j], y27[i]+benchyL_P2[j]
    # Include Pier's corrections
    x_corr = 0.086
    y_corr = 0.077
    if detector == 491:
        x13 = x13 - x_corr
        x15 = x15 - x_corr
        x17 = x17 - x_corr
        y13 = y13 - y_corr
        y15 = y15 - y_corr
        y17 = y17 - y_corr
        x23 = x23 - x_corr
        x25 = x25 - x_corr
        x27 = x27 - x_corr
        y23 = y23 - y_corr
        y25 = y25 - y_corr
        y27 = y27 - y_corr
    elif detector == 492:
        x13 = x13 + x_corr
        x15 = x15 + x_corr
        x17 = x17 + x_corr
        y13 = y13 + y_corr
        y15 = y15 + y_corr
        y17 = y17 + y_corr
        x23 = x23 + x_corr
        x25 = x25 + x_corr
        x27 = x27 + x_corr
        y23 = y23 + y_corr
        y25 = y25 + y_corr
        y27 = y27 + y_corr        
    return x13,y13, x23,y23, x15,y15, x25,y25, x17,y17, x27,y27
    
    
def read_star_param_files(detector, test_case, path4starfiles, paths_list):
    """ This function reads the corresponding star parameters file and returns the data for P1 and P2. """
    cases_list = ["Scene1_slow_real", "Scene1_slow_nonoise", "Scene1_rapid_real", "Scene1_rapid_nonoise",
                  "Scene2_slow_real", "Scene2_slow_nonoise", "Scene2_rapid_real", "Scene2_rapid_nonoise"]
    bench_dirs = [[path4starfiles+paths_list[0], path4starfiles+paths_list[4]],
                  [path4starfiles+paths_list[1], path4starfiles+paths_list[5]],
                  [path4starfiles+paths_list[2], path4starfiles+paths_list[6]],
                  [path4starfiles+paths_list[3], path4starfiles+paths_list[7]],
                  [path4starfiles+paths_list[8], path4starfiles+paths_list[12]],
                  [path4starfiles+paths_list[9], path4starfiles+paths_list[13]],
                  [path4starfiles+paths_list[10], path4starfiles+paths_list[14]],
                  [path4starfiles+paths_list[11], path4starfiles+paths_list[15]]]
    for i, case in enumerate(cases_list):
        if case in test_case:
            dirs4test = bench_dirs[i]
    # NOTE: These positions are all 0-indexed!
    indexing_fix = 0.0
    #    xL:  x-coordinate of the left edge of the postge stamp in the full image (range 0-2047)
    #    xR: x-coord of right edge of the postage stamp
    #    yL: y-coord of the lower edge of the postage stamp
    #    yU:  y-coord of the upper edge of the postage stamp
    # Load parameters of Position 1
    star_param_txt = os.path.join(dirs4test[0],"star parameters.txt")
    benchmark_dataP1 = np.loadtxt(star_param_txt, skiprows=3, unpack=True)
    #bench_star, quadrant, star_in_quad, x_491, y_491, x_492, y_492, V2, V3, xL, xR, yL, yU = benchmark_data
    bench_starP1, _, _, x_491P1, y_491P1, x_492P1, y_492P1, V2P1, V3P1, xLP1, _, yLP1, _ = benchmark_dataP1
    if detector == 491:
        bench_P1 = [bench_starP1, x_491P1+indexing_fix, y_491P1+indexing_fix, V2P1, V3P1, xLP1+indexing_fix, yLP1+indexing_fix]
    elif detector == 492:
        bench_P1 = [bench_starP1, x_492P1+indexing_fix, y_492P1+indexing_fix, V2P1, V3P1, xLP1+indexing_fix, yLP1+indexing_fix]        
    # Load parameters of Position 2
    star_param_txt = os.path.join(dirs4test[1],"star parameters.txt")
    benchmark_dataP2 = np.loadtxt(star_param_txt, skiprows=3, unpack=True)
    #bench_star, quadrant, star_in_quad, x_491, y_491, x_492, y_492, V2, V3, xL, xR, yL, yU = benchmark_data
    bench_starP2, _, _, x_491P2, y_491P2, x_492P2, y_492P2, V2P2, V3P2, xLP2, _, yLP2, _ = benchmark_dataP2
    if detector == 491:
        bench_P2 = [bench_starP2, x_491P2+indexing_fix, y_491P2+indexing_fix, V2P2, V3P2, xLP2+indexing_fix, yLP2+indexing_fix]
    elif detector == 492:
        bench_P2 = [bench_starP2, x_492P2+indexing_fix, y_492P2+indexing_fix, V2P2, V3P2, xLP2+indexing_fix, yLP2+indexing_fix]        
    benchmark_data = [bench_P1, bench_P2]
    return benchmark_data
    
    
def compare2ref(case, path4starfiles, paths_list, bench_stars, benchV2, benchV3, stars, V2in, V3in, arcsecs=True):
    """ This function obtains the differences of the input arrays with the reference or benchmark data. """
    # calculate the differences with respect to the benchmark data
    multiply_by = 1.0          # keep differences in degrees
    if arcsecs:
        multiply_by = 3600.0   # to convert from degrees to arcsecs
    if len(stars) == len(bench_stars):   # for the fixed and None background case
        diffV2 = (benchV2 - V2in) * multiply_by
        diffV3 = (benchV3 - V3in) * multiply_by
        bench_V2_list = benchV2.tolist()
        bench_V3_list = benchV3.tolist()
    else:                               # for the fractional background case
        bench_V2_list, bench_V3_list = [], []
        diffV2, diffV3 = [], []
        for i, s in enumerate(stars):
            if s in bench_stars:
                j = bench_stars.tolist().index(s)
                dsV2 = (benchV2[j] - V2in[i]) * multiply_by
                dsV3 = (benchV3[j] - V3in[i]) * multiply_by 
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


def get_best_fracvalue(stars, diffV2, diffV3):
    """ This function determines which background fractional value has the smallest difference 
    with respect true sky positions. """
    frac_values = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
    min_diff_fracvalueV2, min_diff_fracvalueV3 = [], []
    best_fracval_V2, best_fracval_V3 = [], []
    fracbg_values = 11
    counter = 0
    print (len(diffV2))
    for i, _ in enumerate(diffV2):
        frac_diffs2analyzeV2 = []
        frac_diffs2analyzeV3 = []
        counter = counter + 1
        if counter <= fracbg_values:
            frac_diffs2analyzeV2.append(diffV2[i])
            frac_diffs2analyzeV3.append(diffV3[i])
        print (stars[i], frac_diffs2analyzeV2, frac_diffs2analyzeV3)
        if counter == fracbg_values:
            counter = 0
            minV2 = min(frac_diffs2analyzeV2)
            minV3 = min(frac_diffs2analyzeV3)
            minV2_idx = frac_diffs2analyzeV2.index(minV2)
            minV3_idx = frac_diffs2analyzeV3.index(minV3)
            for i, fv in enumerate(frac_values):
                if minV2_idx == i:
                    min_diff_fvV2 = fv
                if minV3_idx == i:
                    min_diff_fvV3 = fv
                min_diff_fracvalueV2.append(min_diff_fvV2)
                min_diff_fracvalueV3.append(min_diff_fvV3)
    # now repeat the value 11 times for that star number
    counter = 0
    for mdV2, mdV3 in zip(min_diff_fracvalueV2, min_diff_fracvalueV3):
        counter = counter + 1            
        if counter <= fracbg_values:
            best_fracval_V2.append(mdV2)
            best_fracval_V3.append(mdV3)
        if counter == fracbg_values:
            counter = 0
    return best_fracval_V2, best_fracval_V3


#######################################################################################################################

#  --> CODE

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
    print ("* studying case: ", case)
    
    # get the benchmark data
    benchmark_data = read_star_param_files(detector, case, path4starfiles, paths_list)
    bench_P1, bench_P2 = benchmark_data
    bench_starP1, bench_xP1, bench_yP1, bench_V2P1, bench_V3P1, bench_xLP1, bench_yLP1 = bench_P1
    bench_starP2, bench_xP2, bench_yP2, bench_V2P2, bench_V3P2, bench_xLP2, bench_yLP2 = bench_P2
    if single_case:
        bench_stars = bench_starP1.tolist()
        star_idx = bench_stars.index(single_case)
        bench_starP1 = np.array([bench_starP1[star_idx]])
        bench_xP1 = np.array([bench_xP1[star_idx]])
        bench_yP1 = np.array([bench_yP1[star_idx]])
        bench_V2P1 = np.array([bench_V2P1[star_idx]])
        bench_V3P1 = np.array([bench_V3P1[star_idx]])
        bench_xLP1 = np.array([bench_xLP1[star_idx]])
        bench_yLP1 = np.array([bench_yLP1[star_idx]])
        bench_starP2 = np.array([bench_starP2[star_idx]])
        bench_xP2 = np.array([bench_xP2[star_idx]])
        bench_yP2 = np.array([bench_yP2[star_idx]])
        bench_V2P2 = np.array([bench_V2P2[star_idx]])
        bench_V3P2 = np.array([bench_V3P2[star_idx]])
        bench_xLP2 = np.array([bench_xLP2[star_idx]])
        bench_yLP2 = np.array([bench_yLP2[star_idx]])
    #avg_benchX = (bench_xP1 + bench_xP2)/2.0
    #avg_benchY = (bench_yP1 + bench_yP2)/2.0
    avg_benchV2 = (bench_V2P1 + bench_V2P2)/2.0
    avg_benchV3 = (bench_V3P1 + bench_V3P2)/2.0
        
    # read the measured detector centroids
    data = np.loadtxt(infile, skiprows=2, usecols=(0,1,2,3,4,5,6,7,8,9,10,11,12,13), unpack=True)
    stars, bg_value, x13, y13, x15, y15, x17, y17, x23, y23, x25, y25, x27, y27 = data
    if single_case:
        star_idx = stars.tolist().index(single_case)
        if "frac" in case:
            bg_value = np.array(bg_value[ (stars==stars[star_idx]) & (stars <= stars[star_idx+11]) ])
            x13 = np.array(x13[ (stars==stars[star_idx]) & (stars <= stars[star_idx+11]) ])
            y13 = np.array(y13[ (stars==stars[star_idx]) & (stars <= stars[star_idx+11]) ])
            x15 = np.array(x15[ (stars==stars[star_idx]) & (stars <= stars[star_idx+11]) ])
            y15 = np.array(y15[ (stars==stars[star_idx]) & (stars <= stars[star_idx+11]) ])
            x17 = np.array(x17[ (stars==stars[star_idx]) & (stars <= stars[star_idx+11]) ])
            y17 = np.array(y17[ (stars==stars[star_idx]) & (stars <= stars[star_idx+11]) ])
            x23 = np.array(x23[ (stars==stars[star_idx]) & (stars <= stars[star_idx+11]) ])
            y23 = np.array(y23[ (stars==stars[star_idx]) & (stars <= stars[star_idx+11]) ])
            x25 = np.array(x25[ (stars==stars[star_idx]) & (stars <= stars[star_idx+11]) ])
            y25 = np.array(y25[ (stars==stars[star_idx]) & (stars <= stars[star_idx+11]) ])
            x27 = np.array(x27[ (stars==stars[star_idx]) & (stars <= stars[star_idx+11]) ])
            y27 = np.array(y27[ (stars==stars[star_idx]) & (stars <= stars[star_idx+11]) ])
            stars = np.array(stars[ (stars==stars[star_idx]) & (stars <= stars[star_idx+11]) ])
        else:
            stars = np.array([stars[star_idx]])
            bg_value = np.array([bg_value[star_idx]])
            x13 = np.array([x13[star_idx]])
            y13 = np.array([y13[star_idx]])
            x15 = np.array([x15[star_idx]])
            y15 = np.array([y15[star_idx]])
            x17 = np.array([x17[star_idx]])
            y17 = np.array([y17[star_idx]])
            x23 = np.array([x23[star_idx]])
            y23 = np.array([y23[star_idx]])
            x25 = np.array([x25[star_idx]])
            y25 = np.array([y25[star_idx]])
            x27 = np.array([x27[star_idx]])
            y27 = np.array([y27[star_idx]])
            
    
    if debug or single_case:
        print ("Check that read BENCHMARK values correspond to expected for case: ", case)
        print ("Star, xP1, yP1, V2P1, V3P1, xLP1, yLP1")
        print (bench_starP1[0], bench_xP1[0], bench_yP1[0], bench_V2P1[0], bench_V3P1[0], bench_xLP1[0], bench_yLP1[0])
        print ("Star, xP2, yP2, V2P2, V3P2, xLP2, yLP2")
        print (bench_starP2[0], bench_xP2[0], bench_yP2[0], bench_V2P2[0], bench_V3P2[0], bench_xLP2[0], bench_yLP2[0])
        print ("Check that read MEASURED values correspond to expected for the same case: ", case)
        print ("Star, BG, x13, y13, x15, y15, x17, y17, LoLeftP1 (x, y), TrueP1 (x, y)")
        print (stars[0], bg_value[0], x13[0], y13[0], x15[0], y15[0], x17[0], y17[0], bench_xLP1, bench_yLP1, bench_xP1, bench_yP1)
        print ("Star, BG, x23, y23, x25, y25, x27, y27, LoLeftP2 (x, y), TrueP2 (x, y)")
        print (stars[0], x23[0], y23[0], x25[0], y25[0], x27[0], y27[0], bench_xLP2, bench_yLP2, bench_xP2, bench_yP2)
        raw_input(" * press enter to continue... \n")

    # convert from 32x32 pixel to full detector coordinates
    P1P2data = [x13,y13, x23,y23, x15,y15, x25,y25, x17,y17, x27,y27]
    benchmark_xLyL_P1 = [bench_xLP1, bench_yLP1]
    benchmark_xLyL_P2 = [bench_xLP2, bench_yLP2]
    x13,y13, x23,y23, x15,y15, x25,y25, x17,y17, x27,y27 = convert2fulldetector(detector, stars, P1P2data, bench_starP1, benchmark_xLyL_P1, benchmark_xLyL_P2, Pier_corr=Pier_corr)
    
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
    T1_diffV2_3, T1_diffV3_3, bench_V2_list, bench_V3_list = compare2ref(case, path4starfiles, paths_list, bench_starP1, avg_benchV2, avg_benchV3, stars, T1_V2_3, T1_V3_3, arcsecs=diffs_in_arcsecs)
    T1_diffV2_5, T1_diffV3_5, _, _ = compare2ref(case, path4starfiles, paths_list, bench_starP1, avg_benchV2, avg_benchV3, stars, T1_V2_5, T1_V3_5, arcsecs=diffs_in_arcsecs)
    T1_diffV2_7, T1_diffV3_7, _, _ = compare2ref(case, path4starfiles, paths_list, bench_starP1, avg_benchV2, avg_benchV3, stars, T1_V2_7, T1_V3_7, arcsecs=diffs_in_arcsecs)
    # get the minimum of the differences
    T1_min_diff, T1_counter = get_mindiff(T1_diffV2_3, T1_diffV2_5, T1_diffV2_7)
    # get the fractional value that has the smaller difference
    #if "frac" in case:
    #    T1_best_fracval_V2, T1_best_fracval_V3 = get_best_fracvalue(stars, T1_diffV2_3, T1_diffV3_3)
    #    print (T1_best_fracval_V2, T1_best_fracval_V3)
    #    raw_input()
    # calculate standard deviations and means
    if not single_case:    
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
    T2_diffV2_3, T2_diffV3_3, _, _ = compare2ref(case, path4starfiles, paths_list, bench_starP1, avg_benchV2, avg_benchV3, stars, T2_V2_3, T2_V3_3, arcsecs=diffs_in_arcsecs)
    T2_diffV2_5, T2_diffV3_5, _, _ = compare2ref(case, path4starfiles, paths_list, bench_starP1, avg_benchV2, avg_benchV3, stars, T2_V2_5, T2_V3_5, arcsecs=diffs_in_arcsecs)
    T2_diffV2_7, T2_diffV3_7, _, _ = compare2ref(case, path4starfiles, paths_list, bench_starP1, avg_benchV2, avg_benchV3, stars, T2_V2_7, T2_V3_7, arcsecs=diffs_in_arcsecs)
    # get the minimum of the differences
    T2_min_diff, T2_counter = get_mindiff(T2_diffV2_3, T2_diffV2_5, T2_diffV2_7)
    # get the fractional value that has the smaller difference
    #if "frac" in case:
    #    best_fracval_V2, best_fracval_V3 = get_best_fracvalue(T2_diffV2_3, T2_diffV3_3)
    # calculate standard deviations and means
    if not single_case:    
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
    T3_diffV2_13, T3_diffV3_13, _, _ = compare2ref(case, path4starfiles, paths_list, bench_starP1, bench_V2P1, bench_V3P1, stars, T3_V2_13, T3_V3_13, arcsecs=diffs_in_arcsecs)
    T3_diffV2_23, T3_diffV3_23, _, _ = compare2ref(case, path4starfiles, paths_list, bench_starP1, bench_V2P2, bench_V3P2, stars, T3_V2_23, T3_V3_23, arcsecs=diffs_in_arcsecs)
    T3_diffV2_15, T3_diffV3_15, _, _ = compare2ref(case, path4starfiles, paths_list, bench_starP1, bench_V2P1, bench_V3P1, stars, T3_V2_15, T3_V3_15, arcsecs=diffs_in_arcsecs)
    T3_diffV2_25, T3_diffV3_25, _, _ = compare2ref(case, path4starfiles, paths_list, bench_starP1, bench_V2P2, bench_V3P2, stars, T3_V2_25, T3_V3_25, arcsecs=diffs_in_arcsecs)
    T3_diffV2_17, T3_diffV3_17, _, _ = compare2ref(case, path4starfiles, paths_list, bench_starP1, bench_V2P1, bench_V3P1, stars, T3_V2_17, T3_V3_17, arcsecs=diffs_in_arcsecs)
    T3_diffV2_27, T3_diffV3_27, _, _ = compare2ref(case, path4starfiles, paths_list, bench_starP1, bench_V2P2, bench_V3P2, stars, T3_V2_27, T3_V3_27, arcsecs=diffs_in_arcsecs)
    # get the minimum of the differences
    T3_min_diff1, T3_counter1 = get_mindiff(T3_diffV2_13, T3_diffV2_15, T3_diffV2_17)
    T3_min_diff2, T3_counter2 = get_mindiff(T3_diffV2_23, T3_diffV2_25, T3_diffV2_27)
    # calculate standard deviations and means
    if not single_case:    
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
    
    if debug or single_case:
        print ("TEST 1: ")
        print ("transformations: detector (avgx, avgy),  sky (V2, V3),  true (avgV2, avgV3)")
        print ("            ChBx3: ", avgx3, avgy3, T1_V2_3, T1_V3_3, avg_benchV2, avg_benchV3)
        print ("            ChBx5: ", avgx5, avgy5, T1_V2_5, T1_V3_5, avg_benchV2, avg_benchV3)
        print ("            ChBx7: ", avgx7, avgy7, T1_V2_7, T1_V3_7, avg_benchV2, avg_benchV3)
        print ("TEST 2: ")
        print ("transformations: detector P1 and P2 (x, y),  sky (avgV2, avgV3),  true (avgV2, avgV3)")
        print ("            ChBx3: ", x13, y13, x23, y23, T2_V2_3, T2_V3_3, avg_benchV2, avg_benchV3)
        print ("            ChBx5: ", x15, y15, x25, y25, T2_V2_5, T2_V3_5, avg_benchV2, avg_benchV3)
        print ("            ChBx7: ", x17, y17, x27, y27, T2_V2_7, T2_V3_7, avg_benchV2, avg_benchV3)
        print ("TEST 3: ")
        print ("transformations: detector P1 and P2 (x, y),  sky P1 and P2 (V2, V3),  true P1 and P2 (V2, V3)")
        print ("            ChBx3: ", x13, y13, x23, y23, T3_V2_13, T3_V3_13, T3_V2_23, T3_V3_23, bench_V2P1, bench_V3P1, bench_V2P2, bench_V3P2)
        print ("            ChBx5: ", x15, y15, x25, y25, T3_V2_13, T3_V3_13, T3_V2_23, T3_V3_23, bench_V2P1, bench_V3P1, bench_V2P2, bench_V3P2)
        print ("            ChBx7: ", x17, y17, x27, y27, T3_V2_13, T3_V3_13, T3_V2_23, T3_V3_23, bench_V2P1, bench_V3P1, bench_V2P2, bench_V3P2)
        raw_input(" * press enter to continue... \n")

    # Print results to screen and save into a text file if told so
    # Text file 1
    line0 = "{}".format("Differences = diffs = True_Positions - Measured_Positions")
    if diffs_in_arcsecs:
        line0bis = "{}".format("*** diffs are in units of arcsecs")
    else:
        line0bis = "{}".format("*** diffs are in units of degrees")
    line1 = "{}".format("Test1: average P1 and P2, transform to V2-V3, calculate differences")
    if not single_case:    
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
                        "AvgTrue_Pos", "Diff_Chbx_3", "Diff_Chbx_5", "Diff_Chbx_7", "MinDiff")
        line5 = "{:>25} {:>17} {:>22} {:>17} {:>22} {:>17} {:>17} {:>11} {:>18} {:>17} {:>22} {:>17} {:>22} {:>17}".format(
                        "x", "y", "x", "y", "x", "y", "x", "y", "x", "y", "x", "y", "x", "y")
    else:
        line4 = "{:<5} {:<20} {:<40} {:<40} {:<23} {:<7}".format("Star", "BG_value", "Checkbox=3", "Checkbox=5", "Checkbox=7", "MinDiff")
        line5 = "{:>25} {:>17} {:>22} {:>17} {:>22} {:>17}".format("x", "y", "x", "y", "x", "y")        
    print (line0)
    print (line0bis)
    print (line1)
    if not single_case:    
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
        txt_out = path4results+"Test1_results_"+case+".txt"
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
        to.write(line3bisA+"\n")
        to.write(line4+"\n")
        to.write(line5+"\n")
    for i, st in enumerate(stars):
        st = int(st)
        if show_positions:
            line6 = "{:<5} {:<5} {:>20}  {:<20} {:>18}  {:<20} {:>18}  {:<20} {:>10}  {:<14} {:>18}  {:<20} {:>18}  {:<20} {:>18}  {:<20} {:<7}".format(
                        st, bg_value[i], 
                        T1_V2_3[i], T1_V3_3[i], T1_V2_5[i], T1_V3_5[i], T1_V2_7[i], T1_V3_7[i], 
                        bench_V2_list[i], bench_V3_list[i],
                        T1_diffV2_3[i], T1_diffV3_3[i], T1_diffV2_5[i], T1_diffV3_5[i], T1_diffV2_7[i], T1_diffV3_7[i], T1_min_diff[i])
        else:
            line6 = "{:<5} {:<5} {:>20}  {:<20} {:>18}  {:<20} {:>18}  {:<20} {:<7}".format(st, bg_value[i], 
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
    if diffs_in_arcsecs:
        line0bis = "{}".format("*** diffs are in units of arcsecs")
    else:
        line0bis = "{}".format("*** diffs are in units of degrees")
    line1 = "{}".format("Test2: P1 P2, average positions in V2-V3, calculate differences")
    if not single_case:    
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
                        "AvgTrue_Pos", "Diff_Chbx_3", "Diff_Chbx_5", "Diff_Chbx_7", "MinDiff")
        line5 = "{:>25} {:>17} {:>22} {:>17} {:>22} {:>17} {:>17} {:>11} {:>18} {:>17} {:>22} {:>17} {:>22} {:>17}".format(
                        "x", "y", "x", "y", "x", "y", "x", "y", "x", "y", "x", "y", "x", "y")
    else:
        line4 = "{:<5} {:<20} {:<40} {:<40} {:<23} {:<7}".format("Star", "BG_value", "Checkbox=3", "Checkbox=5", "Checkbox=7", "MinDiff")
        line5 = "{:>25} {:>17} {:>22} {:>17} {:>22} {:>17}".format("x", "y", "x", "y", "x", "y")        
    print (line0)
    print (line0bis)
    print (line1)
    if not single_case:    
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
        txt_out = path4results+"Test2_results_"+case+".txt"
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
    line0 = "{}".format("Differences = True_Positions - Measured_Positions")
    if diffs_in_arcsecs:
        line0bis = "{}".format("*** diffs are in units of arcsecs")
    else:
        line0bis = "{}".format("*** diffs are in units of degrees")
    line1 = "{}".format("Test3: P1 and P2, transform to V2-V3 space individually, calculate differences position to position")
    if not single_case:    
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
    print (line0bis)
    print (line1)
    if not single_case:    
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
        txt_out = path4results+"Test3_results_"+case+".txt"
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
    #else:
    #    raw_input(" * Press enter to continue... \n")
    
print ("\n Script 'comparison2sky.py' finished! ")
