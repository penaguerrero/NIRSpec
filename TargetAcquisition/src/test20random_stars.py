from __future__ import print_function, division
from glob import glob
import numpy as np
#import os
import collections
import random
import PIL.Image as Image
# other code
import coords_transform as ct
import testing_functions as tf 
import least_squares_iterate as lsi
print("Modules correctly imported! \n")

"""
This script tests can choose 20 random stars from either detector and runs the following
transformations for the selected (or all) test(s) for the a given set of 20 stars:  

 TEST1 - Average positions P1 and P2, transform to V2-V3 space, and compare to average 
         reference positions (V2-V3 space)
 TEST2 - Transform individual positions P1 and P2 to V2-V3 space, average V2-V3 space 
         positions, and compare to average reference positions.
 TEST3 - Transform P1 and P2 individually to V2-V3 space and compare star by star and 
         position by position.

Outputs:
    - display images with true and calculated centroid
    - text file with 
"""


#######################################################################################################################


# general settings
random_sample = True     # choose a random sample of 20 stars from either detector: True or False
# if wanting a specific sample of stars, input integer numbers into following list
# Known bad stars in X and Y: 103, 105, 106, 112, 134, 152, 156, 170, 188
stars_sample = []
test2perform = "all"     # string, type "all" for 3 tests or "T1", "T2", "T3" for test 1, 2, and 3, respectively
Nsigma = 3               # N-sigma rejection of bad stars, integer or float
max_iterations = 10      # Max number of iterations for N-sigma function, integer
scene = 1                # integer or string, scene=1 is constant Mag 23, scene=2 is stars with Mag 18-23
# Select EITHER a particular shutter velocity to study OR choose "all"
# to study both rapid and slow shutters
shutters = "rapid"       # string, shutter velocity: "rapid", "slow", "all"
bkgd_method = "None"     # background to test, string: "all", "None", "fixed", "frac"  
noise = "real"           # string, noise level: "nonoise" or "real"
filter_input = "F140X"   # Filter, string: for now only test case is "F140X"
show_display = False     # Show display of resulting positions: True or False
save_txt_file = False    # Save text file with resulting transformations: True or False
Pier_corr = True         # Include Pier's corrections to measured positions
show_positions = False   # Print positions on file and screen: True or False
tilt = False             # tilt angle: True or False
debug = False            # See screen print statements for intermediate answers: True or False 
diffs_in_arcsecs = True  # Print the differences in arcsecs? True or False (=degrees) 


#######################################################################################################################

#  --> FUNCTIONS        
    
def get_mindiff(d1, d2, d3):
    """ This function determines the minimum difference from checkboxes 3, 5, and 7,
    and counts the number of repetitions.  """
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


def find_best_fracbgvalue(diffV2, diffV3):
    """ This function uses the next 2 functions to determine the fractional background value that has
    the smallest difference with respect to true values for all 100 stars. 
    * Disregard fractional backgrounds of 0 and 1.
    Input:      - diffV2, diffV3  = numpy arrays of the fractional background True-Measured V2 
                                    and V3 positions for the 100 stars
    Output:     - list_best_frbg_value_V2, list_best_frbg_value_V3 = list of best fractional backgrounds 
                  (the best value is repeated 11 times for each star, so that the length is the same to
                  that of the differences arrays)
    """
    list_best_frbg_value_V2, list_best_frbg_value_V3 = [], []
    # Divide the whole array into a list of arrays according to star number
    arrsV2_list, arrsV3_list = divide_array_by_starNo(diffV2, diffV3)
    # Obtain the list of best fractional background for each star in the list
    for slice_diffV2, slice_diffV3 in zip(arrsV2_list, arrsV3_list):
        best_frbg_value_V2, best_frbg_value_V3 = get_best_fracvalue(slice_diffV2, slice_diffV3)
        list_best_frbg_value_V2.append(best_frbg_value_V2)
        list_best_frbg_value_V3.append(best_frbg_value_V3)
    # Flatten the list of lists into a single list 
    list_best_frbg_value_V2 = [item for sublist in list_best_frbg_value_V2 for item in sublist]
    list_best_frbg_value_V3 = [item for sublist in list_best_frbg_value_V3 for item in sublist]
    counterV2 = collections.Counter(list_best_frbg_value_V2)
    counterV3 = collections.Counter(list_best_frbg_value_V3)
    return list_best_frbg_value_V2, list_best_frbg_value_V3, counterV2, counterV3

def divide_array_by_starNo(arrX, arrY):
    """ This function divides the array of fractional background into slices corresponding to star number. 
    Input:     - arr        =  numpy array to be divided
    Output:    - arrsX_list  =  list of numpy X arrays corresponding to star number
               - arrsY_list  =  list of numpy Y arrays corresponding to star number
    """
    arrsX_list, arrsY_list = [], []
    frac_values = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
    counter = 0
    star_arrayX, star_arrayY = np.array([]), np.array([])
    for itemX, itemY in zip(arrX, arrY):
        counter = counter + 1
        if counter <= len(frac_values):
            star_arrayX = np.append(star_arrayX, itemX)
            star_arrayY = np.append(star_arrayY, itemY)
        if counter == len(frac_values):
            arrsX_list.append(star_arrayX)
            arrsY_list.append(star_arrayY)
            counter = 0
            star_arrayX, star_arrayY = np.array([]), np.array([])
    return arrsX_list, arrsY_list
    
def get_best_fracvalue(slice_diffV2, slice_diffV3):
    """ This function determines which background fractional value has the smallest difference 
    with respect true sky positions. 
    Input:     diffV2, diffV3 = numpy arrays of the fractional background True-Measured V2 
                                and V3 positions for the same star (each array has 11 elements)
    Output:    best_frbg_value_V2 = list of same minimum difference of True-Measured V2
               best_frbg_value_V3 = list of same minimum difference of True-Measured V3
    """
    frac_values = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
    min_diff_V2 = min(slice_diffV2)
    min_diff_V3 = min(slice_diffV3)
    #print (slice_diffV2)
    #print (slice_diffV3)
    #print ("mins for fractional backgrounds in V2 and V3: ", min_diff_V2, min_diff_V3)
    #raw_input()
    best_frbg_value_V2, best_frbg_value_V3 = [], []
    for dv2, dv3 in zip(slice_diffV2, slice_diffV3):
        if dv2 == min_diff_V2:
            idx2 = slice_diffV2.tolist().index(min_diff_V2)
            # make sure the smallest value is not 0.0 or 1.0
            if idx2==0 or idx2==10:
                min_diff_V2 = second_smallest(slice_diffV2)
                idx2 = slice_diffV2.tolist().index(min_diff_V2)
            #print ("best fractional background value for V2 is: ", frac_values[idx])
        if dv3 == min_diff_V3:
            idx3 = slice_diffV3.tolist().index(min_diff_V3)
            # make sure the smallest value is not 0.0 or 1.0
            if idx3==0 or idx3==10:
                min_diff_V3 = second_smallest(slice_diffV3)
                idx3 = slice_diffV3.tolist().index(min_diff_V3)
            #print ("best fractional background value for V3 is: ", frac_values[idx])
    # repeat each value 11 times
    for _ in frac_values:
        best_frbg_value_V2.append(frac_values[idx2])
        best_frbg_value_V3.append(frac_values[idx3])    
    return best_frbg_value_V2, best_frbg_value_V3

def second_smallest(numbers):
    m1, m2 = float('inf'), float('inf')
    for x in numbers:
        if x <= m1:
            m1, m2 = x, m1
        elif x < m2:
            m2 = x
    return m2


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
        print ("            ChBx3: ", avgx3, avgy3, T1_V2_3, T1_V3_3, avg_benchV2, avg_benchV3)
        print ("            ChBx5: ", avgx5, avgy5, T1_V2_5, T1_V3_5, avg_benchV2, avg_benchV3)
        print ("            ChBx7: ", avgx7, avgy7, T1_V2_7, T1_V3_7, avg_benchV2, avg_benchV3)
        raw_input(" * press enter to continue... \n")
    # Organize results
    T1_transformations = [T1_V2_3, T1_V3_3, T1_V2_5, T1_V3_5, T1_V2_7, T1_V3_7]
    T1_diffs = [T1_diffV2_3, T1_diffV3_3, T1_diffV2_5, T1_diffV3_5, T1_diffV2_7, T1_diffV3_7]
    T1_benchVs_list = [T1bench_V2_list, T1bench_V3_list]
    return T1_transformations, T1_diffs, T1_benchVs_list
    
def runTest1_and_append_results(data4test1, Vs, diffs, benchVs):
    """ This function runs the test for the specified detector and sliced arrays, and appends it to the results. """
    detector, transf_direction, case, stars, P1P2data, bench_starP1, avg_benchV23, LoLeftCornersP1, LoLeftCornersP2, Pier_corr = data4test1
    T1_V2_3, T1_V3_3, T1_V2_5, T1_V3_5, T1_V2_7, T1_V3_7 = Vs
    T1_diffV2_3, T1_diffV3_3, T1_diffV2_5, T1_diffV3_5, T1_diffV2_7, T1_diffV3_7 = diffs
    T1bench_V2_list, T1bench_V3_list = benchVs
    # convert from 32x32 pixel to full detector coordinates
    x13,y13, x23,y23, x15,y15, x25,y25, x17,y17, x27,y27 = tf.convert2fulldetector(detector, stars, P1P2data, bench_starP1, LoLeftCornersP1, LoLeftCornersP2, Pier_corr=Pier_corr)
    P1P2data = [x13,y13, x23,y23, x15,y15, x25,y25, x17,y17, x27,y27]
    transformations, diffs, benchVs_list = TEST1(detector, transf_direction, stars, case, bench_starP1, avg_benchV23, P1P2data)
    V2_3, V3_3, V2_5, V3_5, V2_7, V3_7 = transformations
    diffV2_3, diffV3_3, diffV2_5, diffV3_5, diffV2_7, diffV3_7 = diffs
    bench_V2_list, bench_V3_list = benchVs_list
    for v23, v33, v25, v35, v27, v37 in zip(V2_3, V3_3, V2_5, V3_5, V2_7, V3_7):
        T1_V2_3.append(v23)
        T1_V3_3.append(v33)
        T1_V2_5.append(v25)
        T1_V3_5.append(v35)
        T1_V2_7.append(v27)
        T1_V3_7.append(v37)
    for dv23, dv33, dv25, dv35, dv27, dv37 in zip(diffV2_3, diffV3_3, diffV2_5, diffV3_5, diffV2_7, diffV3_7):
        T1_diffV2_3.append(dv23)
        T1_diffV3_3.append(dv33)
        T1_diffV2_5.append(dv25)
        T1_diffV3_5.append(dv35)
        T1_diffV2_7.append(dv27)
        T1_diffV3_7.append(dv37)
    for bv2, bv3 in zip(bench_V2_list, bench_V3_list):
        T1bench_V2_list.append(bv2)
        T1bench_V3_list.append(bv3)
    Vs = [T1_V2_3, T1_V3_3, T1_V2_5, T1_V3_5, T1_V2_7, T1_V3_7]
    diffs = [T1_diffV2_3, T1_diffV3_3, T1_diffV2_5, T1_diffV3_5, T1_diffV2_7, T1_diffV3_7]
    benchVs = [T1bench_V2_list, T1bench_V3_list]
    return Vs, diffs, benchVs

def runTEST1(detectors, transf_direction, case, stars, P1P2data, bench_starP1, trueVsP1, trueVsP2, LoLeftCornersP1, LoLeftCornersP2, Pier_corr):
    """ This function runs the test 1 for both detectors and returns the results for the 20 star sample """
    x13,y13, x23,y23, x15,y15, x25,y25, x17,y17, x27,y27 = P1P2data
    bench_V2P1, bench_V3P1 = trueVsP1
    bench_V2P2, bench_V3P2 = trueVsP2
    avg_benchV2 = (bench_V2P1 + bench_V2P2)/2.0
    avg_benchV3 = (bench_V3P1 + bench_V3P2)/2.0
    T1_V2_3, T1_V3_3, T1_V2_5, T1_V3_5, T1_V2_7, T1_V3_7 = [], [], [], [], [], [] 
    T1_diffV2_3, T1_diffV3_3, T1_diffV2_5, T1_diffV3_5, T1_diffV2_7, T1_diffV3_7 = [], [], [], [], [], []
    T1bench_V2_list, T1bench_V3_list = [], []
    Vs = [T1_V2_3, T1_V3_3, T1_V2_5, T1_V3_5, T1_V2_7, T1_V3_7]
    diffs = [T1_diffV2_3, T1_diffV3_3, T1_diffV2_5, T1_diffV3_5, T1_diffV2_7, T1_diffV3_7]
    benchVs = [T1bench_V2_list, T1bench_V3_list]
    # Find the index at which to change detector
    for st in stars:
        if st >= 100:
            change_detector_idx = stars.tolist().index(st)
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
    d2avg_benchV23 = [avg_benchV2[:change_detector_idx], avg_benchV3[:change_detector_idx]]
    d2stars = stars[:change_detector_idx]
    d2LoLeftCornersP1 = [LoLeftCornersP1[0][:change_detector_idx], LoLeftCornersP1[1][:change_detector_idx]]
    d2LoLeftCornersP2 = [LoLeftCornersP2[0][:change_detector_idx], LoLeftCornersP2[1][:change_detector_idx]]
    data4test1 = [detector, transf_direction, case, d2stars, P1P2data, d2bench_starP1, d2avg_benchV23, d2LoLeftCornersP1, d2LoLeftCornersP2, Pier_corr]
    Vs, diffs, benchVs = runTest1_and_append_results(data4test1, Vs, diffs, benchVs)
    # detector 491
    detector = detectors[0]  
    d1x13, d1y13 = x13[change_detector_idx:], y13[change_detector_idx:]
    d1x23, d1y23 = x23[change_detector_idx:], y23[change_detector_idx:]
    d1x15, d1y15 = x15[change_detector_idx:], y15[change_detector_idx:]
    d1x25, d1y25 = x25[change_detector_idx:], y25[change_detector_idx:]
    d1x17, d1y17 = x17[change_detector_idx:], y17[change_detector_idx:]
    d1x27, d1y27 = x27[change_detector_idx:], y27[change_detector_idx:]
    P1P2data = [d1x13, d1y13, d1x23, d1y23, d1x15, d1y15, d1x25, d1y25, d1x17, d1y17, d1x27, d1y27]
    d1bench_starP1 = bench_starP1[change_detector_idx:]
    d1avg_benchV23 = [avg_benchV2[change_detector_idx:], avg_benchV3[change_detector_idx:]]
    d1stars = stars[change_detector_idx:]
    d1LoLeftCornersP1 = [LoLeftCornersP1[0][change_detector_idx:], LoLeftCornersP1[1][change_detector_idx:]]
    d1LoLeftCornersP2 = [LoLeftCornersP2[0][change_detector_idx:], LoLeftCornersP2[1][change_detector_idx:]]
    data4test1 = [detector, transf_direction, case, d1stars, P1P2data, d1bench_starP1, d1avg_benchV23, d1LoLeftCornersP1, d1LoLeftCornersP2, Pier_corr]
    Vs, diffs, benchVs = runTest1_and_append_results(data4test1, Vs, diffs, benchVs)
    resultsTEST1 = [Vs, diffs, benchVs]
    return resultsTEST1


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
        print ("            ChBx3: ", x13, y13, x23, y23, T2_V2_3, T2_V3_3, avg_benchV2, avg_benchV3)
        print ("            ChBx5: ", x15, y15, x25, y25, T2_V2_5, T2_V3_5, avg_benchV2, avg_benchV3)
        print ("            ChBx7: ", x17, y17, x27, y27, T2_V2_7, T2_V3_7, avg_benchV2, avg_benchV3)
        raw_input(" * press enter to continue... \n")
    # Organize results
    T2_transformations = [T2_V2_3, T2_V3_3, T2_V2_5, T2_V3_5, T2_V2_7, T2_V3_7]
    T2_diffs = [T2_diffV2_3, T2_diffV3_3, T2_diffV2_5, T2_diffV3_5, T2_diffV2_7, T2_diffV3_7]
    T2_benchVs_list = [T2bench_V2_list, T2bench_V3_list]
    return T2_transformations, T2_diffs, T2_benchVs_list
    

def TEST3(detector, transf_direction, stars, bench_Vs, case, x13, x23, y13, y23, x15, y15, x25, y25, x17, y17, x27, y27):
    # TEST 3: (a) Transform P1 and P2 individually to V2-V3 (b) compare star by star and position by position
    bench_V2P1, bench_V3P1, bench_V2P2, bench_V3P2 = bench_Vs
    # Step (a) - transformations
    T3_V2_13, T3_V3_13 = ct.coords_transf(transf_direction, detector, filter_input, x13, y13, tilt, debug)
    T3_V2_15, T3_V3_15 = ct.coords_transf(transf_direction, detector, filter_input, x15, y15, tilt, debug)
    T3_V2_17, T3_V3_17 = ct.coords_transf(transf_direction, detector, filter_input, x17, y17, tilt, debug)
    T3_V2_23, T3_V3_23 = ct.coords_transf(transf_direction, detector, filter_input, x23, y23, tilt, debug)
    T3_V2_25, T3_V3_25 = ct.coords_transf(transf_direction, detector, filter_input, x25, y25, tilt, debug)
    T3_V2_27, T3_V3_27 = ct.coords_transf(transf_direction, detector, filter_input, x27, y27, tilt, debug)
    # Step (b) - comparison
    T3_diffV2_13, T3_diffV3_13, T3bench_V2_listP1, T3bench_V3_listP1 = tf.compare2ref(case, path4starfiles, paths_list, bench_starP1, bench_V2P1, bench_V3P1, stars, T3_V2_13, T3_V3_13, arcsecs=diffs_in_arcsecs)
    T3_diffV2_23, T3_diffV3_23, T3bench_V2_listP2, T3bench_V3_listP2 = tf.compare2ref(case, path4starfiles, paths_list, bench_starP1, bench_V2P2, bench_V3P2, stars, T3_V2_23, T3_V3_23, arcsecs=diffs_in_arcsecs)
    T3_diffV2_15, T3_diffV3_15, _, _ = tf.compare2ref(case, path4starfiles, paths_list, bench_starP1, bench_V2P1, bench_V3P1, stars, T3_V2_15, T3_V3_15, arcsecs=diffs_in_arcsecs)
    T3_diffV2_25, T3_diffV3_25, _, _ = tf.compare2ref(case, path4starfiles, paths_list, bench_starP1, bench_V2P2, bench_V3P2, stars, T3_V2_25, T3_V3_25, arcsecs=diffs_in_arcsecs)
    T3_diffV2_17, T3_diffV3_17, _, _ = tf.compare2ref(case, path4starfiles, paths_list, bench_starP1, bench_V2P1, bench_V3P1, stars, T3_V2_17, T3_V3_17, arcsecs=diffs_in_arcsecs)
    T3_diffV2_27, T3_diffV3_27, _, _ = tf.compare2ref(case, path4starfiles, paths_list, bench_starP1, bench_V2P2, bench_V3P2, stars, T3_V2_27, T3_V3_27, arcsecs=diffs_in_arcsecs)
    if debug:
        print ("TEST 3: ")
        print ("transformations: detector P1 and P2 (x, y),  sky P1 and P2 (V2, V3),  true P1 and P2 (V2, V3)")
        print ("            ChBx3: ", x13, y13, x23, y23, T3_V2_13, T3_V3_13, T3_V2_23, T3_V3_23, bench_V2P1, bench_V3P1, bench_V2P2, bench_V3P2)
        print ("            ChBx5: ", x15, y15, x25, y25, T3_V2_13, T3_V3_13, T3_V2_23, T3_V3_23, bench_V2P1, bench_V3P1, bench_V2P2, bench_V3P2)
        print ("            ChBx7: ", x17, y17, x27, y27, T3_V2_13, T3_V3_13, T3_V2_23, T3_V3_23, bench_V2P1, bench_V3P1, bench_V2P2, bench_V3P2)
        raw_input(" * press enter to continue... \n")
    # Organize results
    T3_transformationsP1 = [T3_V2_13, T3_V3_13, T3_V2_15, T3_V3_15, T3_V2_17, T3_V3_17]
    T3_transformationsP2 = [T3_V2_23, T3_V3_23, T3_V2_25, T3_V3_25, T3_V2_27, T3_V3_27]
    T3_diffsP1 = [T3_diffV2_13, T3_diffV3_13, T3_diffV2_15, T3_diffV3_15, T3_diffV2_17, T3_diffV3_17]
    T3_diffsP2 = [T3_diffV2_23, T3_diffV3_23, T3_diffV2_25, T3_diffV3_25, T3_diffV2_27, T3_diffV3_27]
    T3_benchVs_list = [T3bench_V2_listP1, T3bench_V3_listP1, T3bench_V2_listP2, T3bench_V3_listP2]
    return T3_transformationsP1, T3_diffsP1, T3_transformationsP2, T3_diffsP2, T3_benchVs_list


def get_stats(case, T_transformations, T_diffs, T_benchVs_list, Nsigma, max_iterations):
    T_V2_3, T_V3_3, T_V2_5, T_V3_5, T_V2_7, T_V3_7 = T_transformations
    T_diffV2_3, T_diffV3_3, T_diffV2_5, T_diffV3_5, T_diffV2_7, T_diffV3_7 = T_diffs
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
    T_min_diff, T_counter = get_mindiff(T_diffV2_3, T_diffV2_5, T_diffV2_7)
    # get the fractional value that has the smaller difference
    if "frac" in case:
        Tlist_best_frbg_value_V2_3, Tlist_best_frbg_value_V3_3, TcounterV2_3, TcounterV3_3 = find_best_fracbgvalue(T_diffV2_3, T_diffV3_3)
        Tlist_best_frbg_value_V2_5, Tlist_best_frbg_value_V3_5, TcounterV2_5, TcounterV3_5 = find_best_fracbgvalue(T_diffV2_5, T_diffV3_5)
        Tlist_best_frbg_value_V2_7, Tlist_best_frbg_value_V3_7, TcounterV2_7, TcounterV3_7 = find_best_fracbgvalue(T_diffV2_7, T_diffV3_7)
    # to express in arcsecs multiply by 3600.0
    TLSdeltas_3, TLSsigmas_3, TLSlines2print_3 = lsi.ls_fit_iter(max_iterations, T_V2_3*3600.0, T_V3_3*3600.0, Tbench_V2*3600.0, Tbench_V3*3600.0)
    TLSdeltas_5, TLSsigmas_5, TLSlines2print_5 = lsi.ls_fit_iter(max_iterations, T_V2_5*3600.0, T_V3_5*3600.0, Tbench_V2*3600.0, Tbench_V3*3600.0)
    TLSdeltas_7, TLSsigmas_7, TLSlines2print_7 = lsi.ls_fit_iter(max_iterations, T_V2_7*3600.0, T_V3_7*3600.0, Tbench_V2*3600.0, Tbench_V3*3600.0)
    # Do N-sigma rejection
    TsigmaV2_3, TmeanV2_3, TsigmaV3_3, TmeanV3_3, TnewV2_3, TnewV3_3, Tniter_3, Tlines2print_3 = tf.Nsigma_rejection(Nsigma, T_diffV2_3, T_diffV3_3, max_iterations)
    TsigmaV2_5, TmeanV2_5, TsigmaV3_5, TmeanV3_5, TnewV2_5, TnewV3_5, Tniter_5, Tlines2print_5 = tf.Nsigma_rejection(Nsigma, T_diffV2_5, T_diffV3_5, max_iterations)
    TsigmaV2_7, TmeanV2_7, TsigmaV3_7, TmeanV3_7, TnewV2_7, TnewV3_7, Tniter_7, Tlines2print_7 = tf.Nsigma_rejection(Nsigma, T_diffV2_7, T_diffV3_7, max_iterations)
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
    results_stats = [st_devsAndMeans, diff_counter, bench_values, sigmas_deltas, sigma_reject]
    if "frac" in case:
        best_frac_values = [Tlist_best_frbg_value_V2_3, Tlist_best_frbg_value_V3_3, TcounterV2_3, TcounterV3_3,
                            Tlist_best_frbg_value_V2_5, Tlist_best_frbg_value_V3_5, TcounterV2_5, TcounterV3_5,
                            Tlist_best_frbg_value_V2_7, Tlist_best_frbg_value_V3_7, TcounterV2_7, TcounterV3_7]
        results_stats = [st_devsAndMeans, diff_counter, bench_values, sigmas_deltas, sigma_reject, best_frac_values]
    return results_stats


def dat4case(case, measured_centroids491, measured_centroids492):
    scene, shutters, noise, bkgd_method = case
    scene = "scene"+str(scene)
    # Create the 200 element arrays
    stars1, stars2 = np.array([]), np.array([])
    bg1, x13, y13, x15, y15, x17, y17, factor1 = np.array([]), np.array([]), np.array([]), np.array([]), np.array([]), np.array([]), np.array([]), np.array([])
    bg2, x23, y23, x25, y25, x27, y27, factor2 = np.array([]), np.array([]), np.array([]), np.array([]), np.array([]), np.array([]), np.array([]), np.array([])
    for mc1, mc2 in zip(measured_centroids491, measured_centroids492):
        if scene in mc1 and shutters in mc1 and noise in mc1 and bkgd_method in mc1:   # detector 491
            # Position 1
            d1stars1, d1bg1, d1x13, d1y13, d1x15, d1y15, d1x17, d1y17, d1factor1 = np.loadtxt(mc1, skiprows=3, usecols=(0,1,2,3,4,5,6,7,12), unpack=True)
            if scene in mc1 and shutters in mc1 and noise in mc1 and bkgd_method in mc1 and "shifted" in mc1:    # Position 2 (or shifted position)
                # Position 2
                d1stars2, d1bg2, d1x23, d1y23, d1x25, d1y25, d1x27, d1y27, d1factor2 = np.loadtxt(mc1, skiprows=3, usecols=(0,1,2,3,4,5,6,7,12), unpack=True)
        if scene in mc2 and shutters in mc2 and noise in mc2 and bkgd_method in mc2:   # detector 492
            # Position 1
            d2stars1, d2bg1, d2x13, d2y13, d2x15, d2y15, d2x17, d2y17, d2factor1 = np.loadtxt(mc2, skiprows=3, usecols=(0,1,2,3,4,5,6,7,12), unpack=True)
            if scene in mc2 and shutters in mc2 and noise in mc2 and bkgd_method in mc2 and "shifted" in mc2:   # detector 492
                # Position 2
                d2stars2, d2bg2, d2x23, d2y23, d2x25, d2y25, d2x27, d2y27, d2factor2 = np.loadtxt(mc2, skiprows=3, usecols=(0,1,2,3,4,5,6,7,12), unpack=True)
    # Position 1
    stars1 = np.append(stars1, d2stars1)  # append fist stars of detector 492, then 491
    stars1 = np.append(stars1, d1stars1)
    bg1 = np.append(bg1, d2bg1)
    bg1 = np.append(bg1, d1bg1)
    factor1 = np.append(factor1, d2factor1)
    factor1 = np.append(factor1, d1factor1)
    x13 = np.append(x13, d2x13)
    x13 = np.append(x13, d1x13)
    y13 = np.append(y13, d2y13)
    y13 = np.append(y13, d1y13)
    x15 = np.append(x15, d2x15)
    x15 = np.append(x15, d1x15)
    y15 = np.append(y15, d2y15)
    y15 = np.append(y15, d1y15)
    x17 = np.append(x17, d2x17)
    x17 = np.append(x17, d1x17)
    y17 = np.append(y17, d2y17)
    y17 = np.append(y17, d1y17)
    print ("len(stars1): ", len(stars1))
    # Position 2
    stars2 = np.append(stars2, d2stars2)  # append fist stars of detector 492, then 491
    stars2 = np.append(stars2, d1stars2)
    bg2 = np.append(bg2, d2bg2)
    bg2 = np.append(bg2, d1bg2)
    factor2 = np.append(factor2, d2factor2)
    factor2 = np.append(factor2, d1factor2)
    x23 = np.append(x23, d2x23)
    x23 = np.append(x23, d1x23)
    y23 = np.append(y23, d2y23)
    y23 = np.append(y23, d1y23)
    x25 = np.append(x25, d2x25)
    x25 = np.append(x25, d1x25)
    y25 = np.append(y25, d2y25)
    y25 = np.append(y25, d1y25)
    x27 = np.append(x27, d2x27)
    x27 = np.append(x27, d1x27)
    y27 = np.append(y27, d2y27)
    y27 = np.append(y27, d1y27)
    # organize return arrays
    measured_positions1 = [stars1, bg1, factor1, x13, y13, x15, y15, x17, y17]
    measured_positions2 = [stars2, bg2, factor2, x23, y23, x25, y25, x27, y27]
    return measured_positions1, measured_positions2


def get_sample_data4case(star_idx_list, case2study, measured_centroids491, measured_centroids492):    
    # get the 200 element arrays for measured centroid positions
    measured_positions1, measured_positions2 = dat4case(case2study, measured_centroids491, measured_centroids492)
    allstars1, allbg1, allfactor1, allx13, ally13, allx15, ally15, allx17, ally17 = measured_positions1
    allstars2, allbg2, allfactor2, allx23, ally23, allx25, ally25, allx27, ally27 = measured_positions2
    # get the data for the sample of stars
    stars1, stars2 = np.array([]), np.array([])
    bg1, x13, y13, x15, y15, x17, y17, factor1 = np.array([]), np.array([]), np.array([]), np.array([]), np.array([]), np.array([]), np.array([]), np.array([])
    bg2, x23, y23, x25, y25, x27, y27, factor2 = np.array([]), np.array([]), np.array([]), np.array([]), np.array([]), np.array([]), np.array([]), np.array([])
    for i in star_idx_list:
        # Positions 1 and 2
        stars1 = np.append(stars1, allstars1[i])
        bg1 = np.append(bg1, allbg1[i])
        factor1 = np.append(factor1, allfactor1[i])
        x13 = np.append(x13, allx13[i])
        y13 = np.append(y13, ally13[i])
        x15 = np.append(x15, allx15[i])
        y15 = np.append(y15, ally15[i])
        x17 = np.append(x17, allx17[i])
        y17 = np.append(y17, ally17[i])
        stars2 = np.append(stars2, allstars2[i])
        bg2 = np.append(bg2, allbg2[i])
        factor2 = np.append(factor2, allfactor2[i])
        x23 = np.append(x23, allx23[i])
        y23 = np.append(y23, ally23[i])
        x25 = np.append(x25, allx25[i])
        y25 = np.append(y25, ally25[i])
        x27 = np.append(x27, allx27[i])
        y27 = np.append(y27, ally27[i])
    # organize return arrays
    sample_pos1 = [stars1, bg1, factor1, x13, y13, x15, y15, x17, y17]
    sample_pos2 = [stars2, bg2, factor2, x23, y23, x25, y25, x27, y27]
    return sample_pos1, sample_pos2


def get_figs(star_idx_list, case, centroid_figs491, centroid_figs492):
    scene, shutters, noise, _ = case
    scene = "Scene"+str(scene)
    allfigs1, allfigs2, figs1, figs2 = [], [], [], []
    for subdir1, subdir2 in zip(centroid_figs491, centroid_figs492):
        # detector 492
        if scene in subdir2 and shutters in subdir2 and noise in subdir2:  
            figs = glob(subdir2+"/*")
            for fig in figs:
                # Position 1
                if "shifted" not in fig:
                    allfigs1.append(fig)
                # Position 2
                if "shifted" in fig:
                    allfigs2.append(fig)
        # detector 491
        if scene in subdir1 and shutters in subdir1 and noise in subdir1:  
            figs = glob(subdir1+"/*")
            for fig in figs:
                # Position 1
                if "shifted" not in fig:
                    allfigs1.append(fig)
                # Position 2
                if "shifted" in fig:
                    allfigs2.append(fig)
    # get the sample of stars
    for i in star_idx_list:
        figs1.append(allfigs1[i])
        figs2.append(allfigs2[i])
    return figs1, figs2


def display_figs(figs1, figs2):
    # display the figures
    for f1, f2 in zip(figs1, figs2):
        # Position 1
        image = Image.open(f1)
        image.show()
        image.close()
        # Position 2
        image = Image.open(f2)
        image.show()
        image.close()
        raw_input("press enter to continue... ")


def show_star_displays(star_idx_list, case, centroid_figs491, centroid_figs492):
    figs1, figs2 = get_figs(star_idx_list, case, centroid_figs491, centroid_figs492)
    display_figs(figs1, figs2)


def put_lists_in_correct_order(resultsTEST, ):
    if "frac" in case:
        st_devsAndMeans, diff_counter, bench_values, sigmas_deltas, sigma_reject = resultsTEST


#######################################################################################################################

#  --> CODE

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

detectors = [491, 492]

# Stars of detector 491 and 492
stars_detectors = range(1, 201)

if random_sample:
    # select 20 stars from 1 to 200
    stars_sample = []
    for i in range(20):
        random_star = random.choice(stars_detectors)
        stars_sample.append(random_star)
    # make sure that there are no repetitions
    list(set(stars_sample))
    while len(stars_sample) != 20:
        random_star = random.choice(stars_detectors)
        stars_sample.append(random_star)  
        list(set(stars_sample))      
else:
    stars_sample = stars_sample
# order the star list 
stars_sample.sort(key=lambda xx: xx)

# Define the paths for results 
path4results = "../results20randomstars/"

# Set the case to study according to the selected scene
case2study = "Scene"+str(scene)+"_"

# get the benchmark data according to Scene selected
benchmark_data = tf.read_star_param_files(case2study)
bench_P1, bench_P2 = benchmark_data
allbench_starP1, allbench_xP1, allbench_yP1, allbench_V2P1, allbench_V3P1, allbench_xLP1, allbench_yLP1 = bench_P1
allbench_starP2, allbench_xP2, allbench_yP2, allbench_V2P2, allbench_V3P2, allbench_xLP2, allbench_yLP2 = bench_P2
allbench_stars = allbench_starP1.tolist()

# get the index for the sample stars
star_idx_list = []
print (stars_sample)
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

# get the lists of the text files for measured centroid positions
measured_centroids491 = glob(path4starfiles+"detector_491_resulting_centroid_txt_files_redo/*")
measured_centroids492 = glob(path4starfiles+"detector_492_resulting_centroid_txt_files_redo/*")

# get the lists of the paths for centroid displays
centroid_figs491 = glob(path4starfiles+"detector_491_centroid_figs_redo/*")
centroid_figs492 = glob(path4starfiles+"detector_492_centroid_figs_redo/*")

# specify the case a bit more
shutters_list = ["rapid", "slow"]
bkgd_method_list = ["None", "fixed", "frac"]
noise_list = ["nonoise", "real"]

# for a single case
if shutters != "all" and bkgd_method != "all":
    case2study = [scene, shutters, noise, bkgd_method]
    case = "Scene"+str(scene)+"_"+shutters+"_"+noise+"_"+bkgd_method
    print ("Analyzing case: ", case)
    sample_pos1, sample_pos2 = get_sample_data4case(star_idx_list, case2study, measured_centroids491, measured_centroids492)
    stars1, bg1, factor1, x13, y13, x15, y15, x17, y17 = sample_pos1
    stars2, bg2, factor2, x23, y23, x25, y25, x27, y27 = sample_pos2
    line0 = "Centroid indexing starting at 1 !"
    line0a = "{:<5} {:<15} {:<16} {:>23} {:>30} {:>44} {:>17} {:>15}".format("Star", "Background", 
                                                                      "Checkbox: 3", "5", "7", 
                                                                      "TruePositions", "LoLeftCoords",
                                                                      "Factor")
    line0b = "{:>25} {:>12} {:>16} {:>14} {:>16} {:>14} {:>16} {:>18} {:>12} {:>10}".format(
                                                                           "x", "y", "x", "y", "x", "y", 
                                                                           "TrueX", "TrueY", "LoLeftX", "LoLeftY")
    print (line0)
    print (line0a)
    print (line0b)
    for i, st in enumerate(stars1): 
        line1 = "{:<5} {:<10} {:<14} {:<16} {:<14} {:<16} {:<14} {:<16} {:<14} {:<16} {:<8} {:<12} {:<10}".format(
                                                                    int(st), bg1[i], 
                                                                    x13[i], y13[i], x15[i], y15[i], x17[i], y17[i],
                                                                    bench_xP1[i]-bench_xLP1[i], bench_yP1[i]-bench_yLP1[i],
                                                                    bench_xLP1[i], bench_yLP1[i],
                                                                    factor1[i])
        print (line1)
            
    # compact results for functions
    P1P2data = [x13,y13, x23,y23, x15,y15, x25,y25, x17,y17, x27,y27]

    # show the displays with positions for the sample of stars
    if show_display:
        show_star_displays(star_idx_list, case, centroid_figs491, centroid_figs492)
    
    # Now run the tests
    transf_direction = "forward"
    
    if test2perform == "T1" or "all":
        # TEST 1: (a) Avg P1 and P2, (b) transform to V2-V3, (c) compare to avg reference positions (V2-V3 space)
        resultsTEST1 = runTEST1(detectors, transf_direction, case, stars1, P1P2data, bench_starP1, trueVsP1, trueVsP2, LoLeftCornersP1, LoLeftCornersP2, Pier_corr)
        T1_transformations, T1_diffs, T1_benchVs_list = resultsTEST1
        T1_V2_3, T1_V3_3, T1_V2_5, T1_V3_5, T1_V2_7, T1_V3_7 = T1_transformations
        T1_diffV2_3, T1_diffV3_3, T1_diffV2_5, T1_diffV3_5, T1_diffV2_7, T1_diffV3_7 = T1_diffs
        T1bench_V2_list, T1bench_V3_list = T1_benchVs_list
        print(len(resultsTEST1))
    exit()
    
"""
    if test2perform == "T2" or "all":
        # TEST 2: (a) Transform individual P1 and P2 to V2-V3, (b) avg V2-V3 space positions, (c) compare to avg reference positions
        resTEST2 = TEST2(detector, transf_direction, avg_benchV23, case, x13, x23, y13, y23, x15, y15, y25, y25, x17, y17, x27, y27)
        resultsTEST2.append(resTEST2)
    
    if test2perform == "T3" or "all":
        # TEST 3: (a) Transform P1 and P2 individually to V2-V3 (b) compare star by star and position by position
        resTEST3 = TEST3(detector, transf_direction, bench_V23, case, x13, x23, y13, y23, x15, y15, y25, y25, x17, y17, x27, y27)
        resultsTEST3.append(resTEST3)
    
    exit()


# Start the loop over the cases to study of scence 1 and 2
#if shutters == "all" and bkgd_method == "all":
    
        
#elif shutters == "all" and bkgd_method != "all":


bench_starP1_list, bench_xP1_list, bench_yP1_list, bench_V2P1_list = [], [], [], []
bench_V3P1_list, bench_xLP1_list, bench_yLP1_list = [], [], []
bench_starP2_list, bench_xP2_list, bench_yP2_list, bench_V2P2_list = [], [], [], []
bench_V3P2_list, bench_xLP2_list, bench_yLP2_list = [], [], []
avg_benchV2_list, avg_benchV3_list = [], []

x13_list,y13_list, x23_list,y23_list = [], [], [], []
x15_list,y15_list, x25_list,y25_list  = [], [], [], []
x17_list,y17_list, x27_list,y27_list = [], [], [], []

resultsTEST1, resultsTEST2, resultsTEST3 = [], [], []






# Define the paths for inputs
for single_case_idx, single_case in enumerate(stars_sample):
    if single_case <= 100:
        detector = 491
    else:
        detector = 492
    path4inputP1P2 = "../PFforMaria/detector_"+str(detector)+"_comparison_txt_positions/"
    # Get list of all input files of the background method given to test
    if bkgd_method == "all":
        input_files_list = glob(path4inputP1P2+"*.txt")
    else:
        input_files_list = glob(path4inputP1P2+"*"+bkgd_method+"*.txt")
    
    # Find stars and their values in the files
    for infile in input_files_list:
        # define the case to work with 
        case = os.path.basename(infile).replace(".txt", "")
        case = case.replace("_centroids", "")
        print ("* studying case: ", case)
        
        # get the benchmark data
        benchmark_data = tf.read_star_param_files(detector, case, path4starfiles, paths_list)
        bench_P1, bench_P2 = benchmark_data
        bench_starP1, bench_xP1, bench_yP1, bench_V2P1, bench_V3P1, bench_xLP1, bench_yLP1 = bench_P1
        bench_starP2, bench_xP2, bench_yP2, bench_V2P2, bench_V3P2, bench_xLP2, bench_yLP2 = bench_P2
        bench_stars = bench_starP1.tolist()
        # get the data for the star in question
        star_idx = bench_stars.index(single_case)
        bench_starP1_single = np.array([bench_starP1[star_idx]])
        bench_xP1_single = np.array([bench_xP1[star_idx]])
        bench_yP1_single = np.array([bench_yP1[star_idx]])
        bench_V2P1_single = np.array([bench_V2P1[star_idx]])
        bench_V3P1_single = np.array([bench_V3P1[star_idx]])
        bench_xLP1_single = np.array([bench_xLP1[star_idx]])
        bench_yLP1_single = np.array([bench_yLP1[star_idx]])
        bench_starP2_single = np.array([bench_starP2[star_idx]])
        bench_xP2_single = np.array([bench_xP2[star_idx]])
        bench_yP2_single = np.array([bench_yP2[star_idx]])
        bench_V2P2_single = np.array([bench_V2P2[star_idx]])
        bench_V3P2_single = np.array([bench_V3P2[star_idx]])
        bench_xLP2_single = np.array([bench_xLP2[star_idx]])
        bench_yLP2_single = np.array([bench_yLP2[star_idx]])
        avg_benchV2_single = (bench_V2P1 + bench_V2P2)/2.0
        avg_benchV3_single = (bench_V3P1 + bench_V3P2)/2.0
        bench_starP1_list.append(bench_starP1_single)
        bench_xP1_list.append(bench_xP1_single)
        bench_yP1_list.append(bench_yP1_single)
        bench_V2P1_list.append(bench_V2P1_single)
        bench_V3P1_list.append(bench_V3P1_single)
        bench_xLP1_list.append(bench_xLP1_single)
        bench_yLP1_list.append(bench_yLP1_single)
        bench_starP2_list.append(bench_starP2_single)
        bench_xP2_list.append(bench_xP2_single)
        bench_yP2_list.append(bench_yP2_single)
        bench_V2P2_list.append(bench_V2P2_single)
        bench_V3P2_list.append(bench_V3P2_single)
        bench_xLP2_list.append(bench_xLP2_single)
        bench_yLP2_list.append(bench_yLP2_single)
        avg_benchV2_list.append(avg_benchV2_single)
        avg_benchV3_list.append(avg_benchV2_single)
            
        # read the measured detector centroids
        data = np.loadtxt(infile, skiprows=2, usecols=(0,1,2,3,4,5,6,7,8,9,10,11,12,13), unpack=True)
        stars, bg_value, x13, y13, x15, y15, x17, y17, x23, y23, x25, y25, x27, y27 = data            
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
        x13_list.append(x13)
        y13_list.append(y13)
        x23_list.append(x23)
        y23_list.append(y23)
        x15_list.append(x15)
        y15_list.append(y15)
        x25_list.append(x25)
        y25_list.append(y25)
        x17_list.append(x17)
        y17_list.append(y17)
        x27_list.append(x27)
        y27_list.append(y27)
        
        if debug:
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
            
    bench_starP1 = np.array(bench_starP1_list)
    bench_xP1 = np.array(bench_xP1_list)
    bench_yP1 = np.array(bench_yP1_list)
    bench_V2P1 = np.array(bench_V2P1_list)
    bench_V3P1 = np.array(bench_V3P1_list)
    bench_xLP1 = np.array(bench_xLP1_list)
    bench_yLP1 = np.array(bench_yLP1_list)
    bench_starP2 = np.array(bench_starP2_list)
    bench_xP2 = np.array(bench_xP2_list)
    bench_yP2 = np.array(bench_yP2_list)
    bench_V2P2 = np.array(bench_V2P2_list)
    bench_V3P2 = np.array(bench_V3P2_list)
    bench_xLP2 = np.array(bench_xLP2_list)
    bench_yLP2 = np.array(bench_yLP2_list)
    avg_benchV2 = np.array(avg_benchV2_list)
    avg_benchV3 = np.array(avg_benchV3_list)
    x13 = np.array(x13_list)
    y13 = np.array(y13_list)
    x23 = np.array(x23_list)
    y23 = np.array(y23_list)
    x15 = np.array(x15_list)
    y15 = np.array(y15_list)
    x25 = np.array(x25_list)
    y25 = np.array(y25_list)
    x17 = np.array(x17_list)
    y17 = np.array(y17_list)
    x27 = np.array(x27_list)
    y27 = np.array(y27_list)
    
    # divide the arrays according to detector
    for i, single_case in enumerate(stars_sample):
        if single_case <= 100:
        else:
            break
    x13_det1, y13_det1, x23_det1, y23_det1 = x13[:i], y13[:i], x23[:i], y23[:i]
    x15_det1, y15_det1, x25_det1, y25_det1 = x15[:i], y15[:i], x25[:i], y25[:i]
    x17_det1, y17_det1, x27_det1, y27_det1 = x17[:i], y17[:i], x27[:i], y27[:i]
    bench_V23_det1 = [bench_V2P1[:i], bench_V3P1[:i], bench_V2P2[:i], bench_V3P2[:i]]
    avg_bench_det1_V23 = [avg_benchV2[:i], avg_benchV3[:i]]
    x13_det2, y13_det2, x23_det2, y23_det2 = x13[i:], y13[i:], x23[i:], y23[i:]
    x15_det2, y15_det2, x25_det2, y25_det2 = x15[i:], y15[i:], x25[i:], y25[i:]
    x17_det2, y17_det2, x27_det2, y27_det2 = x17[i:], y17[i:], x27[i:], y27[i:]
    bench_V23_det1 = [bench_V2P1[i:], bench_V3P1[i:], bench_V2P2[i:], bench_V3P2[i:]]
    avg_bench_det1_V23 = [avg_benchV2[i:], avg_benchV3[i:]]
    
    # Start TESTS the loop
    for infile in input_files_list:
        # define the case to work with 
        case = os.path.basename(infile).replace(".txt", "")
        case = case.replace("_centroids", "")
        print ("* studying case: ", case)
        # convert from 32x32 pixel to full detector coordinates
        P1P2data = [x13,y13, x23,y23, x15,y15, x25,y25, x17,y17, x27,y27]
        benchmark_xLyL_P1 = [bench_xLP1, bench_yLP1]
        benchmark_xLyL_P2 = [bench_xLP2, bench_yLP2]
        x13,y13, x23,y23, x15,y15, x25,y25, x17,y17, x27,y27 = tf.convert2fulldetector(detector, stars, P1P2data, bench_starP1, benchmark_xLyL_P1, benchmark_xLyL_P2, Pier_corr=Pier_corr)
        
        transf_direction = "forward"
        
        if test2perform == "T1" or "all":
            # TEST 1: (a) Avg P1 and P2, (b) transform to V2-V3, (c) compare to avg reference positions (V2-V3 space)
            resTEST1 = TEST1(detector, transf_direction, avg_benchV23, case, x13, x23, y13, y23, x15, y15, y25, y25, x17, y17, x27, y27)
            resultsTEST1.append(resTEST1)
        if test2perform == "T2" or "all":
            # TEST 2: (a) Transform individual P1 and P2 to V2-V3, (b) avg V2-V3 space positions, (c) compare to avg reference positions
            resTEST2 = TEST2(detector, transf_direction, avg_benchV23, case, x13, x23, y13, y23, x15, y15, y25, y25, x17, y17, x27, y27)
            resultsTEST2.append(resTEST2)
        if test2perform == "T3" or "all":
            # TEST 3: (a) Transform P1 and P2 individually to V2-V3 (b) compare star by star and position by position
            resTEST3 = TEST3(detector, transf_direction, bench_V23, case, x13, x23, y13, y23, x15, y15, y25, y25, x17, y17, x27, y27)
            resultsTEST3.append(resTEST3)
            
            
    if test2perform == "T1" or "all":
        st_devsAndMeans, diff_counter, bench_values, sigmas_deltas, sigma_reject = resultsTEST1
        T1stdev_V2_3, T1mean_V2_3, T1stdev_V2_5, T1mean_V2_5, T1stdev_V2_7, T1mean_V2_7, T1stdev_V3_3, T1mean_V3_3, T1stdev_V3_5, T1mean_V3_5, T1stdev_V3_7, T1mean_V3_7 = st_devsAndMeans
        T1_min_diff, T1_counter = diff_counter
        T1bench_V2, T1bench_V3 = bench_values
        T1LSdeltas_3, T1LSsigmas_3, T1LSlines2print_3, T1LSdeltas_5, T1LSsigmas_5, T1LSlines2print_5, T1LSdeltas_7, T1LSsigmas_7, T1LSlines2print_7 = sigmas_deltas
        T1sigmaV2_3, T1meanV2_3, T1sigmaV3_3, T1meanV3_3, T1newV2_3, T1newV3_3, T1niter_3, T1lines2print_3, T1sigmaV2_5, T1meanV2_5, T1sigmaV3_5, T1meanV3_5, T1newV2_5, T1newV3_5, T1niter_5, T1lines2print_5, T1sigmaV2_7, T1meanV2_7, T1sigmaV3_7, T1meanV3_7, T1newV2_7, T1newV3_7, T1niter_7, T1lines2print_7 = sigma_reject
        if "frac" in case:
            st_devsAndMeans, diff_counter, bench_values, sigmas_deltas, sigma_reject, best_frac_values = resultsTEST1
            T1list_best_frbg_value_V2_3, T1list_best_frbg_value_V3_3, T1counterV2_3, T1counterV3_3, T1list_best_frbg_value_V2_5, T1list_best_frbg_value_V3_5, T1counterV2_5, T1counterV3_5, T1list_best_frbg_value_V2_7, T1list_best_frbg_value_V3_7, T1counterV2_7, T1counterV3_7 = best_frac_values
    
    if test2perform == "T2" or "all":
        st_devsAndMeans, diff_counter, bench_values, sigmas_deltas, sigma_reject = resultsTEST2
        T2stdev_V2_3, T2mean_V2_3, T2stdev_V2_5, T2mean_V2_5, T2stdev_V2_7, T2mean_V2_7, T2stdev_V3_3, T2mean_V3_3, T2stdev_V3_5, T2mean_V3_5, T2stdev_V3_7, T2mean_V3_7 = st_devsAndMeans
        T2_min_diff, T2_counter = diff_counter
        T2bench_V2, T2bench_V3 = bench_values
        T2LSdeltas_3, T2LSsigmas_3, T2LSlines2print_3, T2LSdeltas_5, T2LSsigmas_5, T2LSlines2print_5, T2LSdeltas_7, T2LSsigmas_7, T2LSlines2print_7 = sigmas_deltas
        T2sigmaV2_3, T2meanV2_3, T2sigmaV3_3, T2meanV3_3, T2newV2_3, T2newV3_3, T2niter_3, T2lines2print_3, T2sigmaV2_5, T2meanV2_5, T2sigmaV3_5, T2meanV3_5, T2newV2_5, T2newV3_5, T2niter_5, T2lines2print_5, T2sigmaV2_7, T2meanV2_7, T2sigmaV3_7, T2meanV3_7, T2newV2_7, T2newV3_7, T2niter_7, T2lines2print_7 = sigma_reject
        if "frac" in case:
            st_devsAndMeans, diff_counter, bench_values, sigmas_deltas, sigma_reject, best_frac_values = resultsTEST1
            T2list_best_frbg_value_V2_3, T2list_best_frbg_value_V3_3, T2counterV2_3, T2counterV3_3, T2list_best_frbg_value_V2_5, T2list_best_frbg_value_V3_5, T2counterV2_5, T2counterV3_5, T2list_best_frbg_value_V2_7, T2list_best_frbg_value_V3_7, T2counterV2_7, T2counterV3_7 = best_frac_values
    
    if test2perform == "T3" or "all":
        st_devsAndMeans, diff_counter, sigmas_deltas, sigma_reject = resultsTEST3
        T3stdev_V2_13, T3mean_V2_13, T3stdev_V2_15, T3mean_V2_15, T3stdev_V2_17, T3mean_V2_17, T3stdev_V3_13, T3mean_V3_13, T3stdev_V3_15, T3mean_V3_15, T3stdev_V3_17, T3mean_V3_17, T3stdev_V2_23, T3mean_V2_23, T3stdev_V2_25, T3mean_V2_25, T3stdev_V2_27, T3mean_V2_27, T3stdev_V3_23, T3mean_V3_23, T3stdev_V3_25, T3mean_V3_25, T3stdev_V3_27, T3mean_V3_27 = st_devsAndMeans
        T3_min_diff1, T3_counter1, T3_min_diff2, T3_counter2 = diff_counter
        T3LSdeltas_13, T3LSsigmas_13, T3LSlines2print_13, T3LSdeltas_15, T3LSsigmas_15, T3LSlines2print_15, T3LSdeltas_17, T3LSsigmas_17, T3LSlines2print_17, T3LSdeltas_23, T3LSsigmas_23, T3LSlines2print_23,  T3LSdeltas_25, T3LSsigmas_25, T3LSlines2print_25, T3LSdeltas_27, T3LSsigmas_27, T3LSlines2print_27 = sigmas_deltas
        T3sigmaV2_13, T3meanV2_13, T3sigmaV3_13, T3meanV3_13, T3newV2_13, T3newV3_13, T3niter_13, T3lines2print_13, T3sigmaV2_15, T3meanV2_15, T3sigmaV3_15, T3meanV3_15, T3newV2_15, T3newV3_15, T3niter_15, T3lines2print_15, T3sigmaV2_17, T3meanV2_17, T3sigmaV3_17, T3meanV3_17, T3newV2_17, T3newV3_17, T3niter_17, T3lines2print_17, T3sigmaV2_23, T3meanV2_23, T3sigmaV3_23, T3meanV3_23, T3newV2_23, T3newV3_23, T3niter_23, T3lines2print_23, T3sigmaV2_25, T3meanV2_25, T3sigmaV3_25, T3meanV3_25, T3newV2_25, T3newV3_25, T3niter_25, T3lines2print_25, T3sigmaV2_27, T3meanV2_27, T3sigmaV3_27, T3meanV3_27, T3newV2_27, T3newV3_27, T3niter_27, T3lines2print_27 = sigma_reject
        if "frac" in case:
            st_devsAndMeans, diff_counter, sigmas_deltas, sigma_reject, best_frac_values = resultsTEST3
            T3list_best_frbg_value_V2_13, T3list_best_frbg_value_V3_13, T3counterV2_13, T3counterV3_13, T3list_best_frbg_value_V2_15, T3list_best_frbg_value_V3_15, T3counterV2_15, T3counterV3_15, T3list_best_frbg_value_V2_17, T3list_best_frbg_value_V3_17, T3counterV2_17, T3counterV3_17, T3list_best_frbg_value_V2_23, T3list_best_frbg_value_V3_23, T3counterV2_23, T3counterV3_23, T3list_best_frbg_value_V2_25, T3list_best_frbg_value_V3_25, T3counterV2_25, T3counterV3_25, T3list_best_frbg_value_V2_27, T3list_best_frbg_value_V3_27, T3counterV2_27, T3counterV3_27 = best_frac_values

    
    # Print results to screen and save into a text file if told so
    # Text file 1
    line0 = "{}".format("Differences = diffs = True_Positions - Measured_Positions")
    if diffs_in_arcsecs:
        line0bis = "{}".format("*** diffs are in units of arcsecs")
    else:
        line0bis = "{}".format("*** diffs are in units of degrees")
    line1 = "{}\n {}".format("Test1: average P1 and P2, transform to V2-V3, calculate differences",
                             "  * Standard deviations and means ")
    # print regular standard deviations and means
    line2a = "std_dev_V2_3 = {:<20}    std_dev_V3_3 = {:<20}".format(T1stdev_V2_3, T1stdev_V3_3)
    line2b = "std_dev_V2_5 = {:<20}    std_dev_V3_5 = {:<20}".format(T1stdev_V2_5, T1stdev_V3_5)
    line2c = "std_dev_V2_7 = {:<20}    std_dev_V3_7 = {:<20}".format(T1stdev_V2_7, T1stdev_V3_7)
    line3a = "   mean_V2_3 = {:<22}     mean_V3_3 = {:<22}".format(T1mean_V2_3, T1mean_V3_3)
    line3b = "   mean_V2_5 = {:<22}     mean_V3_5 = {:<22}".format(T1mean_V2_5, T1mean_V3_5)
    line3c = "   mean_V2_7 = {:<22}     mean_V3_7 = {:<22}".format(T1mean_V2_7, T1mean_V3_7)
    # Print number of repetitions to find best checkbox
    line3bisA = "Repetitions Diffs: {}".format(T1_counter)
    if "frac" in case:
        line3bisB = "Repetitions Fractional BG values ChBx3: V2-{}".format(T1counterV2_3)
        line3bisC = "                                        V3-{}".format(T1counterV3_3)
        line3bisD = "                                 ChBx5: V2-{}".format(T1counterV2_5)
        line3bisE = "                                        V3-{}".format(T1counterV3_5)
        line3bisF = "                                 ChBx7: V2-{}".format(T1counterV2_7)
        line3bisG = "                                        V3-{}".format(T1counterV3_7)
    if show_positions:
        line4 = "{:<5} {:<20} {:<40} {:<40} {:<38} {:<35} {:<40} {:<40} {:<24} {:<7}".format(
                        "Star", "BG_value", "Avg_Pos_Checkbox_3", "Avg_Pos_Checkbox_5", "Avg_PosCheckbox_7",
                        "AvgTrue_Pos", "Diff_Chbx_3", "Diff_Chbx_5", "Diff_Chbx_7", "MinDiff")
        line5 = "{:>25} {:>17} {:>22} {:>17} {:>22} {:>22} {:>17} {:>11} {:>18} {:>17} {:>22} {:>17} {:>22} {:>17}".format(
                        "x", "y", "x", "y", "x", "y", "x", "y", "x", "y", "x", "y", "x", "y")
    else:
        line4 = "{:<5} {:<20} {:<40} {:<40} {:<23} {:<9}".format("Star", "BG_value", "Checkbox=3", "Checkbox=5", "Checkbox=7", "MinDiff")
        line5 = "{:>25} {:>17} {:>22} {:>17} {:>22} {:>17}".format("x", "y", "x", "y", "x", "y")        
    if "frac" in case:
        line4 = line4 + " {:>23}".format("BestFracVal")
        line5 = line5 + " {:>25} {:>10} {:>10}".format("Cbx3", "Cbx5", "Cbx7")
    print (line0)
    print (line0bis)
    print (line1)
    print (line2a)
    print (line2b)
    print (line2c)
    print (line3a)
    print (line3b)
    print (line3c)
    print (line3bisA)
    if "frac" in case:
        print (line3bisB)
        print (line3bisC)
        print (line3bisD)
        print (line3bisE)
        print (line3bisF)
        print (line3bisG)
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
        # print standard deviations from least squares routine
        to.write("  * From least squares routine:  \n")
        to.write("       Checkbox 3:  \n")
        for line2print in T1LSlines2print_3:
            to.write(line2print+"\n")
        to.write("       Checkbox 5:  \n")
        for line2print in T1LSlines2print_5:
            to.write(line2print+"\n")
        to.write("       Checkbox 7:  \n")
        for line2print in T1LSlines2print_7:
            to.write(line2print+"\n")
        # print standard deviations and means after n-sigma rejection
        to.write(" Checkbox 3:  \n")
        for line2print in T1lines2print_3:
            to.write(line2print+"\n")
        to.write(" Checkbox 5:  \n")
        for line2print in T1lines2print_5:
            to.write(line2print+"\n")
        to.write(" Checkbox 7:  \n")
        for line2print in T1lines2print_7:
            to.write(line2print+"\n")
        to.write(line3bisA+"\n")
        if "frac" in case:
            to.write(line3bisB+"\n")
            to.write(line3bisC+"\n")
            to.write(line3bisD+"\n")
            to.write(line3bisE+"\n")
            to.write(line3bisF+"\n")
            to.write(line3bisG+"\n")
        to.write(line4+"\n")
        to.write(line5+"\n")
    for i, st in enumerate(stars):
        st = int(st)
        if show_positions:
            line6 = "{:<5} {:<5} {:>20}  {:<20} {:>18}  {:<20} {:>18}  {:<20} {:>10}  {:<14} {:>18}  {:<20} {:>18}  {:<20} {:>18}  {:<20} {:<7}".format(
                        st, bg_value[i], 
                        T1_V2_3[i], T1_V3_3[i], T1_V2_5[i], T1_V3_5[i], T1_V2_7[i], T1_V3_7[i], 
                        T1bench_V2_list[i], T1bench_V3_list[i],
                        T1_diffV2_3[i], T1_diffV3_3[i], T1_diffV2_5[i], T1_diffV3_5[i], T1_diffV2_7[i], T1_diffV3_7[i], T1_min_diff[i])
        else:
            line6 = "{:<5} {:<5} {:>20}  {:<20} {:>18}  {:<20} {:>18}  {:<20} {:<7}".format(st, bg_value[i], 
                    T1_diffV2_3[i], T1_diffV3_3[i], T1_diffV2_5[i], T1_diffV3_5[i], T1_diffV2_7[i], T1_diffV3_7[i], T1_min_diff[i])
        if "frac" in case:
            line6 = line6 + " {:<4} {:<5} {:<4} {:<5} {:<4} {:<5}".format(
                            T1list_best_frbg_value_V2_3[i], T1list_best_frbg_value_V3_3[i],
                            T1list_best_frbg_value_V2_5[i], T1list_best_frbg_value_V3_5[i],
                            T1list_best_frbg_value_V2_7[i], T1list_best_frbg_value_V3_7[i])
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
    line2a = "std_dev_V2_3 = {:<20}   std_dev_V3_3 = {:<20}".format(T2stdev_V2_3, T2stdev_V3_3)
    line2b = "std_dev_V2_5 = {:<20}   std_dev_V3_5 = {:<20}".format(T2stdev_V2_5, T2stdev_V3_5)
    line2c = "std_dev_V2_7 = {:<20}   std_dev_V3_7 = {:<20}".format(T2stdev_V2_7, T2stdev_V3_7)
    line3a = "   mean_V2_3 = {:<22}    mean_V3_3 = {:<22}".format(T2mean_V2_3, T2mean_V3_3)
    line3b = "   mean_V2_5 = {:<22}    mean_V3_5 = {:<22}".format(T2mean_V2_5, T2mean_V3_5)
    line3c = "   mean_V2_7 = {:<22}    mean_V3_7 = {:<22}".format(T2mean_V2_7, T2mean_V3_7)
    line3bisA = "Repetitions Diffs1: {}".format(T2_counter)
    if "frac" in case:
        line3bisB = "Repetitions Fractional BG values ChBx3: V2-{}".format(T2counterV2_3)
        line3bisC = "                                        V3-{}".format(T2counterV3_3)
        line3bisD = "                                 ChBx5: V2-{}".format(T2counterV2_5)
        line3bisE = "                                        V3-{}".format(T2counterV3_5)
        line3bisF = "                                 ChBx7: V2-{}".format(T2counterV2_7)
        line3bisG = "                                        V3-{}".format(T2counterV3_7)
    if show_positions:
        line4 = "{:<5} {:<20} {:<40} {:<40} {:<35} {:<30} {:<40} {:<40} {:<23} {:<7}".format(
                        "Star", "BG_value", "Avg_Pos_Checkbox_3", "Avg_Pos_Checkbox_5", "Avg_PosCheckbox_7",
                        "AvgTrue_Pos", "Diff_Chbx_3", "Diff_Chbx_5", "Diff_Chbx_7", "MinDiff")
        line5 = "{:>25} {:>17} {:>22} {:>17} {:>22} {:>17} {:>17} {:>11} {:>18} {:>17} {:>22} {:>17} {:>22} {:>17}".format(
                        "x", "y", "x", "y", "x", "y", "x", "y", "x", "y", "x", "y", "x", "y")
    else:
        line4 = "{:<5} {:<20} {:<40} {:<40} {:<23} {:<9}".format("Star", "BG_value", "Checkbox=3", "Checkbox=5", "Checkbox=7", "MinDiff")
        line5 = "{:>25} {:>17} {:>22} {:>17} {:>22} {:>17}".format("x", "y", "x", "y", "x", "y")        
    if "frac" in case:
        line4 = line4 + " {:>23}".format("BestFracVal")
        line5 = line5 + " {:>25} {:>10} {:>10}".format("Cbx3", "Cbx5", "Cbx7")
    print (line0)
    print (line0bis)
    print (line1)
    print (line2a)
    print (line2b)
    print (line2c)
    print (line3a)
    print (line3b)
    print (line3c)
    print (line3bisA)
    if "frac" in case:
        print (line3bisB)
        print (line3bisC)
        print (line3bisD)
        print (line3bisE)
        print (line3bisF)
        print (line3bisG)
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
        # print standard deviations from least squares routine
        to.write("  * From least squares routine:  \n")
        to.write("       Checkbox 3:  \n")
        for line2print in T2LSlines2print_3:
            to.write(line2print+"\n")
        to.write("       Checkbox 5:  \n")
        for line2print in T2LSlines2print_5:
            to.write(line2print+"\n")
        to.write("       Checkbox 7:  \n")
        for line2print in T2LSlines2print_7:
            to.write(line2print+"\n")
        # print standard deviations and means after n-sigma rejection
        to.write(" Checkbox 3:  \n")
        for line2print in T2lines2print_3:
            to.write(line2print+"\n")
        to.write(" Checkbox 5:  \n")
        for line2print in T2lines2print_5:
            to.write(line2print+"\n")
        to.write(" Checkbox 7:  \n")
        for line2print in T2lines2print_7:
            to.write(line2print+"\n")
        to.write(line3bisA+"\n")
        if "frac" in case:
            to.write(line3bisB+"\n")
            to.write(line3bisC+"\n")
            to.write(line3bisD+"\n")
            to.write(line3bisE+"\n")
            to.write(line3bisF+"\n")
            to.write(line3bisG+"\n")
        to.write(line4+"\n")
        to.write(line5+"\n")
    for i, st in enumerate(stars):
        st = int(st)
        if show_positions:
            line6 = "{:<5} {:<5} {:>20}  {:<20} {:>18}  {:<20} {:>18}  {:<20} {:>10}  {:<14} {:>18}  {:<20} {:>18}  {:<20} {:>18}  {:<20} {:<7}".format(
                        st, bg_value[i], 
                        T2_V2_3[i], T2_V3_3[i], T2_V2_5[i], T2_V3_5[i], T2_V2_7[i], T2_V3_7[i], 
                        T2bench_V2_list[i], T2bench_V3_list[i],
                        T2_diffV2_3[i], T2_diffV3_3[i], T2_diffV2_5[i], T2_diffV3_5[i], T2_diffV2_7[i], T2_diffV3_7[i], T2_min_diff[i])
        else:
            line6 = "{:<5} {:<5} {:>20}  {:<20} {:>18}  {:<20} {:>18}  {:<20} {:<7}".format(st, bg_value[i], 
                    T2_diffV2_3[i], T2_diffV3_3[i], T2_diffV2_5[i], T2_diffV3_5[i], T2_diffV2_7[i], T2_diffV3_7[i], T2_min_diff[i])            
        if "frac" in case:
            line6 = line6 + " {:<4} {:<5} {:<4} {:<5} {:<4} {:<5}".format(
                            T2list_best_frbg_value_V2_3[i], T2list_best_frbg_value_V3_3[i],
                            T2list_best_frbg_value_V2_5[i], T2list_best_frbg_value_V3_5[i],
                            T2list_best_frbg_value_V2_7[i], T2list_best_frbg_value_V3_7[i])
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
    line2a = "std_dev_V2_P1_3 = {:<20}    std_dev_V3_P1_3 = {:<20}".format(T3stdev_V2_13, T3stdev_V3_13)
    line2b = "std_dev_V2_P1_5 = {:<20}    std_dev_V3_P1_5 = {:<20}".format(T3stdev_V2_15, T3stdev_V3_15)
    line2c = "std_dev_V2_P1_7 = {:<20}    std_dev_V3_P1_7 = {:<20}".format(T3stdev_V2_17, T3stdev_V3_17)
    line3a = "   mean_V2_P1_3 = {:<22}     mean_V3_P1_3 = {:<22}".format(T3mean_V2_13, T3mean_V3_13)
    line3b = "   mean_V2_P1_5 = {:<22}     mean_V3_P1_5 = {:<22}".format(T3mean_V2_15, T3mean_V3_15)
    line3c = "   mean_V2_P1_7 = {:<22}     mean_V3_P1_7 = {:<22}".format(T3mean_V2_17, T3mean_V3_17)
    line4a = "std_dev_V2_P2_3 = {:<20}    std_dev_V3_P2_3 = {:<20}".format(T3stdev_V2_23, T3stdev_V3_23)
    line4b = "std_dev_V2_P2_5 = {:<20}    std_dev_V3_P2_5 = {:<20}".format(T3stdev_V2_25, T3stdev_V3_25)
    line4c = "std_dev_V2_P2_7 = {:<20}    std_dev_V3_P2_7 = {:<20}".format(T3stdev_V2_27, T3stdev_V3_27)
    line5a = "   mean_V2_P2_3 = {:<22}     mean_V3_P2_3 = {:<22}".format(T3mean_V2_23, T3mean_V3_23)
    line5b = "   mean_V2_P2_5 = {:<22}     mean_V3_P2_5 = {:<22}".format(T3mean_V2_25, T3mean_V3_25)
    line5c = "   mean_V2_P2_7 = {:<22}     mean_V3_P2_7 = {:<22}".format(T3mean_V2_27, T3mean_V3_27)
    line5bisA = "Repetitions Diffs1: {}".format(T3_counter1)
    line5bisB = "Repetitions Diffs2: {}".format(T3_counter2)
    if "frac" in case:
        line3bisC = "Repetitions Fractional BG values Positions 1 and 2 ChBx3: V2 {}  {}".format(T3counterV2_13, T3counterV2_23)
        line3bisD = "                                                          V3 {}  {}".format(T3counterV3_13, T3counterV3_23)
        line3bisE = "                                                   ChBx5: V2 {}  {}".format(T3counterV2_15, T3counterV2_25)
        line3bisF = "                                                          V3 {}  {}".format(T3counterV3_15, T3counterV3_25)
        line3bisG = "                                                   ChBx7: V2 {}  {}".format(T3counterV2_17, T3counterV2_27)
        line3bisH = "                                                          V3 {}  {}".format(T3counterV3_17, T3counterV3_27)
    if show_positions:
        line6 = "{:<5} {:<15} {:<35} {:<39} {:<37} {:<38} {:<36} {:<38} {:<30} {:<40} {:<40} {:<40} {:<36} {:<38} {:<35} {:<25} {:<7}".format(
                    "Star", "BG_value", "Pos1_Checkbox_3", "Pos1_Checkbox_5", "Pos1_Checkbox_7",
                    "Pos2_Checkbox_3", "Pos2_Checkbox_5", "Pos2_Checkbox_7", 
                    "True_Pos1", "True_Pos2", "Diff1_Chbx_3", "Diff1_Chbx_5", "Diff1_Chbx_7", 
                    "Diff2_Chbx_3", "Diff2_Chbx_5", "Diff2_Chbx_7", "MinDiff1", "MinDiff2")
        line7 = "{:>22} {:>14} {:>22} {:>14} {:>22} {:>14} {:>22} {:>14} {:>22} {:>14} {:>22} {:>14} {:>22} {:>10} {:>22} {:>10} {:>22} {:>14} {:>22} {:>14} {:>22} {:>14} {:>28} {:>14} {:>22} {:>14} {:>22} {:>14}".format(
                        "x", "y", "x", "y", "x", "y", "x", "y", "x", "y", "x", "y", "x", "y",
                        "x", "y", "x", "y", "x", "y", "x", "y", "x", "y", "x", "y", "x", "y", 
                        "x", "y", "x", "y", "x", "y", "x", "y")
    else:
        line6 = "{:<5} {:<20} {:<40} {:<40} {:<40} {:<40} {:<40} {:<25} {:<9} {:<9}".format(
                    "Star", "BG_value", "Diff1_Chbx_3", "Diff1_Chbx_5", "Diff1_Chbx_7",
                    "Diff2_Chbx_3", "Diff2_Chbx_5", "Diff2_Chbx_7", "MinDiff1", "MinDiff2")
        line7 = "{:>25} {:>17} {:>22} {:>17} {:>22} {:>17} {:>22} {:>17} {:>22} {:>17} {:>22} {:>17}".format("x", "y", "x", "y", "x", "y", "x", "y", "x", "y", "x", "y")        
    if "frac" in case:
        line6 = line6 + " {:>33}".format("BestFracVal")
        line7 = line7 + " {:>33} {:>10} {:>10} {:>11} {:>10} {:>10}".format("P13", "P15", "P17", "P23", "P25", "P27")
    print (line0)
    print (line0bis)
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
    if "frac" in case:
        print (line3bisC)
        print (line3bisD)
        print (line3bisE)
        print (line3bisF)
        print (line3bisG)
        print (line3bisH)
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
        # print standard deviations from least squares routine
        to.write("  * From least squares routine:  \n")
        to.write("     Position 1:  \n")
        to.write("       Checkbox 3:  \n")
        for line2print in T3LSlines2print_13:
            to.write(line2print+"\n")
        to.write("       Checkbox 5:  \n")
        for line2print in T3LSlines2print_15:
            to.write(line2print+"\n")
        to.write("       Checkbox 7:  \n")
        for line2print in T3LSlines2print_17:
            to.write(line2print+"\n")
        to.write("     Position 2:  \n")
        to.write("       Checkbox 3:  \n")
        for line2print in T3LSlines2print_23:
            to.write(line2print+"\n")
        to.write("       Checkbox 5:  \n")
        for line2print in T3LSlines2print_25:
            to.write(line2print+"\n")
        to.write("       Checkbox 7:  \n")
        for line2print in T3LSlines2print_27:
            to.write(line2print+"\n")
        # print standard deviations and means after n-sigma rejection
        to.write("     Position 1:  \n")
        to.write("       Checkbox 3:  \n")
        for l2p in T3lines2print_13:
            to.write(l2p+"\n")
        to.write("       Checkbox 5:  \n")
        for l2p in T3lines2print_15:
            to.write(l2p+"\n")
        to.write("       Checkbox 7:  \n")
        for l2p in T3lines2print_17:
            to.write(l2p+"\n")
        to.write("     Position 2:  \n")
        to.write("       Checkbox 3:  \n")
        for l2p in T3lines2print_23:
            to.write(l2p+"\n")
        to.write("       Checkbox 5:  \n")
        for l2p in T3lines2print_25:
            to.write(l2p+"\n")
        to.write("       Checkbox 7:  \n")
        for l2p in T3lines2print_27:
            to.write(l2p+"\n")
        to.write(line3bisA+"\n")
        if "frac" in case:
            to.write(line3bisB+"\n")
            to.write(line3bisC+"\n")
            to.write(line3bisD+"\n")
            to.write(line3bisE+"\n")
            to.write(line3bisF+"\n")
            to.write(line3bisG+"\n")
        to.write(line5bisA+"\n")
        to.write(line5bisB+"\n")
        to.write(line6+"\n")
        to.write(line7+"\n")
    for i, st in enumerate(stars):
        st = int(st)
        if show_positions:
            line8 = "{:<5} {:<5} {:>16}  {:<19} {:>16}  {:<19} {:>16}  {:<19} {:>16}  {:<19} {:>16}  {:<19} {:>16}  {:<19} {:>14}  {:<14}  {:>14}  {:<14} {:>18}  {:<19} {:>18}  {:<19} {:>18}  {:<19} {:>18}  {:<19} {:>18}  {:<19} {:>18}  {:<19} {:<7} {:<7}".format(
                        st, bg_value[i], 
                        T3_V2_13[i], T3_V3_13[i], T3_V2_15[i], T3_V3_15[i], T3_V2_17[i], T3_V3_17[i], 
                        T3_V2_23[i], T3_V3_23[i], T3_V2_25[i], T3_V3_25[i], T3_V2_27[i], T3_V3_27[i], 
                        T3bench_V2_listP1[i], T3bench_V3_listP1[i], T3bench_V2_listP2[i], T3bench_V3_listP2[i],
                        T3_diffV2_13[i], T3_diffV3_13[i], T3_diffV2_15[i], T3_diffV3_15[i], T3_diffV2_17[i], T3_diffV3_17[i],
                        T3_diffV2_23[i], T3_diffV3_23[i], T3_diffV2_25[i], T3_diffV3_25[i], T3_diffV2_27[i], T3_diffV3_27[i],
                        T3_min_diff1[i], T3_min_diff2[i])
        else:
            line8 = "{:<5} {:<5} {:>20}  {:<20} {:>18}  {:<20} {:>18}  {:<20} {:>18}  {:<20} {:>18}  {:<20} {:>18}  {:<20} {:<7} {:<7}".format(st, bg_value[i], 
                    T3_diffV2_13[i], T3_diffV3_13[i], T3_diffV2_15[i], T3_diffV3_15[i], T3_diffV2_17[i], T3_diffV3_17[i], 
                    T3_diffV2_23[i], T3_diffV3_23[i], T3_diffV2_25[i], T3_diffV3_25[i], T3_diffV2_27[i], T3_diffV3_27[i], 
                    T3_min_diff1[i], T3_min_diff2[i])            
        if "frac" in case:
            line8 = line8 + " {:<4} {:<5} {:<4} {:<5} {:<4} {:<6} {:<4} {:<5} {:<4} {:<5} {:<4} {:<5}".format(
                            T3list_best_frbg_value_V2_13[i], T3list_best_frbg_value_V3_13[i],
                            T3list_best_frbg_value_V2_15[i], T3list_best_frbg_value_V3_15[i],
                            T3list_best_frbg_value_V2_17[i], T3list_best_frbg_value_V3_17[i],
                            T3list_best_frbg_value_V2_23[i], T3list_best_frbg_value_V3_23[i],
                            T3list_best_frbg_value_V2_25[i], T3list_best_frbg_value_V3_25[i],
                            T3list_best_frbg_value_V2_27[i], T3list_best_frbg_value_V3_27[i])
        print (line8)
        if save_txt_file:
            to.write(line8+"\n")
    if save_txt_file:
        to.close()
        print (" * Results saved in file: ", txt_out)
    else:
        raw_input(" * Press enter to continue... \n")
"""

print ("\n Script 'comparison2sky.py' finished! ")
