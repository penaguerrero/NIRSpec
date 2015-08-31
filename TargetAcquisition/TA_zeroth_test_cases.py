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

# extra coding used
import testing_functions as tf

print()

# Header
__author__ = "Maria A. Pena-Guerrero"
__version__ = "1.0"

"""

DESCRIPTION:
    This script runs 4 'zeroth' test cases for the Target Locate code jwst_targloc.py 
    written by T. Le Blanc, and modified by M. A. Pena-Guerrero. These fits files are
    a perfect MSA (all shutters open without hot pixels) over 200 stars randomly distributed, 
    100 over detector 491 (extension 1) and 100 over detector 492 (extension 2).
    Test cases are:
        1) All stars with magnitude 23.0
        2) Magnitude range from 18.0 to 23.0
    Tasks performed are:
        - Performing centroiding for widths of 3, 5, and 7 pixels
        - Fractional backgrounds from 0.0 to 1.0 in increments of 0.1, and a fixed 
        or None background=0 case
        - Moments are included in the code though not currently being used
        
"""

# Paths
perfect_scene1 = "PFforMaria/electron_rate_maps/SKY-F140X-MIRROR_MOS_simuTA20150528-F140X-S50-K-AB23_004.erm"
perfect_scene2 = "PFforMaria/electron_rate_maps/SKY-F140X-MIRROR_MOS_simuTA20150528-F140X-S50-K-AB18to23_002.erm"
noisy_scene1 = "PFforMaria/electron_rate_maps/SKY-F140X-MIRROR_MOS_simuTA20150528-F140X-S50-K-AB23_005.erm"
noisy_scene2 = "PFforMaria/electron_rate_maps/SKY-F140X-MIRROR_MOS_simuTA20150528-F140X-S50-K-AB18to23_003.erm"
#                    1               2             3             4
paths_list = [perfect_scene1, perfect_scene2, noisy_scene1, noisy_scene2]

###########################################################################################################

# Set test parameters
path_number = 1             # Select 0 through 4 from paths_list above 
perform_cutouts = True
detector = 491
background_method = None#'frac'  # Select either 'fractional', 'fixed', or None   
use_list_files = False

###########################################################################################################

# --> FUNCTIONS

###########################################################################################################

# ---> CODE

# start the timer to compute the whole running time
start_time = time.time()

# Star numbers of detectors 492 and 491
stars_492 = range(1, 101)
stars_491 = [x+100 for x in range(1, 101)]

# Read fits table with benchmark data
if path_number==1 or path_number==3:
    path2listfile = "PFforMaria/Scene_1_AB23"
    list_file = "simuTA20150528-F140X-S50-K-AB23.list"
    path_scene = "PFforMaria/Scene_1_AB23/NIRSpec_TA_Sim_AB23 first NRS no_noise/postage"
if path_number==2 or path_number==4:
    # Read the text file just written to get the offsets from the "real" positions of the fake stars
    path2listfile = "PFforMaria/Scene_2_AB1823"
    list_file = "simuTA20150528-F140X-S50-K-AB18to23.list"
    path_scene = "PFforMaria/Scene_2_AB1823/NIRSpec_TA_Sim_AB1823 first NRS no_noise/postage"
if use_list_files:   
    lf = os.path.join(path2listfile,list_file)
    star_number, xpos, ypos, factor, mag, bg_method = tf.read_listfile(lf, detector, background_method)
    #print (xpos)
else:
    # Read the star parameters file to compare results
    star_param_txt = os.path.join(path_scene,"star parameters.txt")
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
    #print (true_x)

if perform_cutouts:
    # Read FITS image and get data
    fimg = fits.open(paths_list[path_number])
    fimg.info()
    hdr = fimg[0].header
    data_detector491 = fimg[1].data
    data_detector492 = fimg[2].data
    fimg.close()
    #print("** HEADER: \n", hdr)
    print("** shape data_detector491: ", np.shape(data_detector491))
    print("** shape data_detector492: ", np.shape(data_detector492))
    detector_x, detector_y = data_detector491[0], data_detector491[1]
    for i, v in enumerate(data_detector491[0]):#.ravel()):
        if i != 0:
            print (i, v)
    exit()
    if detector == 492:
        detector_x, detector_y = data_detector492[0], data_detector492[1]
    # With the positions of each reference star, get the cutouts
    for st in stars_detector:
        if st in bench_star:
            idx = np.where(bench_star==st)[0][0]
            xi = true_x[idx]
            yi = true_y[idx]
            ref_star_position = [xi, yi]
            #tf.get_cutouts(ref_star_position, detector_x, detector_y)
            

print ("\n Controlled test script finished. Took  %s  seconds to finish. \n" % ((time.time() - start_time)) )
