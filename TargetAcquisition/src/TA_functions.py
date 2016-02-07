from __future__ import print_function, division
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
import copy
import os
import collections

# Tommy's code
import jwst_targloc as jtl

# Header
__author__ = "Maria A. Pena-Guerrero"
__version__ = "1.0"

"""
This script has Target Acquisition functions frequently used. They are ordered alphabetically.
"""

# FUNCTIONS

def bg_correction(img, bg_method=None, bg_value=None, bg_frac=None, debug=False):
    """
    Subtract a background value from every pixel in the image, based on
    the background method (None, Fixed, or Fraction):
        - If None, the image is used as-is.
        - If Fixed, the given background level (bg_value) is the value
        to be subtracted from each pixel in the image.
        - If Fraction, the given background level is the fraction (e.g. if
        bg_fraction = 0.5, the background is set to the median pixel value
        in the image; if bg_fraction = 0.4, 40% of the pixels have data
        values less than background, while 60% have data values larger than 
        background, and implicitly, the top 20% of the data values are 
        assumed to contain significant counts from astronomical sources, 
        cosmic rays, or hot pixels. See code).

    Keyword arguments:
    img        -- Image 
    bg_method  -- Either None value or string: "fixed" or "frac"
    bg_value   -- Fixed value to subtract from each pixel (this has to 
                  be set if bg_method = "fixed")
    bg_frac    -- Fractional value to subtract from image (this has to 
                  be set if bg_method = "frac")
    
    Output(s):
    img_bgcorr -- The group of 3 background subtracted images 

    Example usage:
    
        >> img_bgcorr = bg_correction(master_img, bg_method='frac', bg_value=0.4)

        Correct each ramped image from background.
    """
    if bg_method is None:
        return img
    
    elif "fix" in bg_method:
        # Check that bg_value is defined
        if bg_value is None:
            print ("ERROR - Background_method set to 'fixed': bg_value needs to be a float number, got None.")
            exit()
        master_img_bgcorr = img - bg_value
        return master_img_bgcorr
    
    elif "frac" in bg_method:
        # Check that bg_value is defined
        if bg_frac is None:
            print ("ERROR - Background_method set to 'fractional': bg_frac needs to be a float number, got None.")
            exit()
        #if bg_frac == 0.0:   #?
        #    return master_img
        # Find the pixel value (bg) that represents that fraction of the population
        img_original = copy.deepcopy(img)
        sorted_img = np.sort(np.ravel(img))   # flatten the image and sort it
        xsize = np.shape(img)[1]
        ysize = np.shape(img)[0]
        idx_bg = np.floor(bg_frac * xsize * ysize)
        # If at the edge, correct
        if idx_bg == np.shape(sorted_img)[0]:
            idx_bg = idx_bg - 1
        bg = sorted_img[idx_bg]
        img_bgcorr = img_original - bg
        # Debugging messages
        if debug:
            print("(bg_correction): xsize = {},  ysize= {}".format(xsize, ysize))
            print("(bg_correction): sorted_img = {}".format(sorted_img))
            print("(bg_correction): idx_bg = {}".format(idx_bg))
            print("(bg_correction): bg = {}".format(bg))
        return img_bgcorr


def centroid2fulldetector(cb_centroid_list, true_center):
    """
    Transform centroid coordinates into full detector coordinates.
    
    KEYWORD ARGUMENTS:
    cb_centroid_list           -- Checkbox based centroid determined by target acquisition (TA) algorithm in 
                                  terms of 32 by 32 pixels for checkbox sizes 3, 5, and 7
    true_center                -- Actual (true) position of star in terms of full detector  
    
    OUTPUT:
    cb_centroid_list_fulldetector  -- List of centroid locations determined with the TA algorithm in 
                                      terms of full detector. List is for positions determined with
                                      3, 5, and 7 checkbox sizes. 
    loleftcoords                   -- Coordinates of the lower left corner of the 32x32 pixel box
    true_center32x32               -- True center given in coordinates of 32x32 pix
    differences_true_TA            -- Difference of true-observed positions       
    """
        
    # Get the lower left corner coordinates in terms of full detector. We subtract 16.0 because indexing
    # from centroid function starts with 1
    corrected_x = true_center[0]
    corrected_y = true_center[1]
    loleft_x = np.floor(corrected_x) - 16.0
    loleft_y = np.floor(corrected_y) - 16.0
    loleftcoords = [loleft_x, loleft_y]
    #print(loleft_x, loleft_y)
    
    # get center in term of 32x32 checkbox
    true_center32x32 = [corrected_x-loleft_x, corrected_y-loleft_y]
    
    # Add lower left corner to centroid location to get it in terms of full detector
    cb_centroid_list_fulldetector = []
    for centroid_location in cb_centroid_list:
        centroid_fulldetector_x = centroid_location[0] + loleft_x
        centroid_fulldetector_y = centroid_location[1] + loleft_y
        centroid_fulldetector = [centroid_fulldetector_x, centroid_fulldetector_y]
        cb_centroid_list_fulldetector.append(centroid_fulldetector)
    corr_cb_centroid_list = cb_centroid_list_fulldetector
    
    # Determine difference between center locations
    differences_true_TA = []
    d3_x = true_center[0] - corr_cb_centroid_list[0][0]
    d3_y = true_center[1] - corr_cb_centroid_list[0][1]
    d3 = [d3_x, d3_y]
    if len(corr_cb_centroid_list) != 1:   # make sure this function works even for one checkbox
        d5_x = true_center[0] - corr_cb_centroid_list[1][0]
        d5_y = true_center[1] - corr_cb_centroid_list[1][1]
        d7_x = true_center[0] - corr_cb_centroid_list[2][0]
        d7_y = true_center[1] - corr_cb_centroid_list[2][1]
        d5 = [d5_x, d5_y]
        d7 = [d7_x, d7_y]
        diffs = [d3, d5, d7]
    else:
        diffs = d3
    differences_true_TA.append(diffs)
    return corr_cb_centroid_list, loleftcoords, true_center32x32, differences_true_TA



def compare2ref(case, bench_stars, benchV2, benchV3, stars, V2in, V3in, arcsecs=True):
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


def display_centroids(detector, st, case, psf, corr_true_center_centroid, 
                      corr_cb_centroid_list, show_disp, vlims=None, savefile=False, 
                      fig_name=None, redos=False, display_master_img=False):  
    if isinstance(st, int): 
        fig_title = "star_"+str(st)+"_"+case
    else:
        fig_title = st
    if vlims is None:
        vlims = (10, 50)
    if display_master_img is not False:
        # Display original image.
        _, ax = plt.subplots(figsize=(8, 8))
        ax.set_title(fig_title+"_original")
        ax.autoscale(enable=False, axis='both')
        ax.imshow(display_master_img, cmap='gray', interpolation='nearest')
        ax.set_ylim(1.0, np.shape(display_master_img)[0])
        ax.set_xlim(1.0, np.shape(display_master_img)[1])
        ax.imshow(display_master_img, cmap='gray', interpolation='nearest', vmin=vlims[0], vmax=vlims[1])
    # Add plot of measured centroids
    fig, ax = plt.subplots(figsize=(8, 8))
    ax.set_title(fig_title)
    ax.autoscale(enable=False, axis='both')
    ax.imshow(psf, cmap='gray', interpolation='nearest')
    ax.set_ylim(1.0, np.shape(psf)[0])
    ax.set_xlim(1.0, np.shape(psf)[1])
    ax.plot(corr_cb_centroid_list[0][0], corr_cb_centroid_list[0][1], marker='*', ms=20, mec='black', mfc='blue', ls='', label='Checkbox=3')
    if len(corr_cb_centroid_list) != 1:
        ax.plot(corr_cb_centroid_list[1][0], corr_cb_centroid_list[1][1], marker='*', ms=17, mec='black', mfc='green', ls='', label='Checkbox=5')
        ax.plot(corr_cb_centroid_list[2][0], corr_cb_centroid_list[2][1], marker='*', ms=15, mec='black', mfc='red', ls='', label='Checkbox=7')
        ax.plot(corr_true_center_centroid[0], corr_true_center_centroid[1], marker='o', ms=8, mec='black', mfc='yellow', ls='', label='True Centroid')
    # Shrink current axis by 10%
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.9, box.height])
    ax.legend(loc='upper right', bbox_to_anchor=(1.26, 1.0), prop={"size":"small"})   # put legend out of the plot box   
    ax.imshow(psf, cmap='gray', interpolation='nearest', vmin=vlims[0], vmax=vlims[1])
    if show_disp:
        plt.show()
    else:
        plt.close('all')
    if savefile:
        path4fig = "../PFforMaria/detector_"+str(detector)+"_centroid_figs"
        if "scene1" in fig_title:
            if "slow" in fig_title:
                if "real" in fig_title:
                    in_dir = "Scene1_slow_real"
                else:
                    in_dir = "Scene1_slow_nonoise"
            elif "rapid" in fig_title:
                if "real" in fig_title:
                    in_dir = "Scene1_rapid_real"
                else:
                    in_dir = "Scene1_rapid_nonoise"
        if "scene2" in fig_title:
            if "slow" in fig_title:
                if "real" in fig_title:
                    in_dir = "Scene2_slow_real"
                else:
                    in_dir = "Scene2_slow_nonoise"
            elif "rapid" in fig_title:
                if "real" in fig_title:
                    in_dir = "Scene2_rapid_real"
                else:
                    in_dir = "Scene2_rapid_nonoise"
        if fig_name is None:
            fig_name = path4fig+in_dir+"/"+fig_title+".jpg"
        if redos:
            fig_name = path4fig+"_redo/"+in_dir+"_redo/"+fig_title+"_redo.jpg"
        fig.savefig(fig_name)
        print ("Figure ", fig_name, " was saved!")
    
    
def display_ns_psf(image, vlim=(), fsize=(8, 8), interp='nearest', \
    title='', cmap='gray', extent=None, savefile=None, cb=False):
    """
    Custom display a PSF generated with WEBBPSF or similar tool.
    A quick tool for displaying NIRSpec images in native size 
    (2048x2048) with additional options for display.
    Keyword arguments:
    image    --  A 2D image to display
    vlim     --  The image range (in terms of image counts) to display.
                 Defaults to empty (), displaying full spectrum.
    fsize    --  Figure image size (in cm?)
    interp   --  Interpolation type. Defaults to 'nearest'.
    title    --  Title for plot. Defaults to ''.
    cmap     --  Color map for plot. Defaults to 'gray'.
    cb       --  Color bar toggle. Defaults to 'False'.
    savefile --  Figure save toggle. Defaults to 'None'. A filename
                 (with directory structure) can be entered to save to
                 given filename.
    """

    # Display PSF (oversampled and detector levels)
    fig, ax = plt.subplots(figsize=fsize)
    ax.set_title(title)
    ax.set_aspect('equal')
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_ylim(0.0, np.shape(image)[0])

    if vlim == ():
        vlim = (image.min(), image.max())
    
    if extent is not None:     
        cax = ax.imshow(image, cmap=cmap, interpolation=interp, vmin=vlim[0], \
              extent=extent-.5, vmax=vlim[1])
    else:
        cax = ax.imshow(image, cmap=cmap, interpolation=interp, vmin=vlim[0], vmax=vlim[1])       
    
    if cb: fig.colorbar(cax, ax=ax, shrink=0.8)
    
    # Added by Maria to see plots when not in Notebook environment
    plt.show()
    
    if savefile is not None:
        fig.savefig(savefile)


def do_Piers_correction(detector, cb_centroid_list):
    xy3, xy5, xy7 = cb_centroid_list
    xy3corr = Pier_correction(detector, xy3)
    xy5corr = Pier_correction(detector, xy5)
    xy7corr = Pier_correction(detector, xy7)
    corr_cb_centroid_list = [xy3corr, xy5corr, xy7corr]
    return corr_cb_centroid_list 
    
    
def find_std(arr):
    """ This function determines the standard deviation of the given array. """
    N = float(len(arr))
    mean = sum(arr) / N
    diff2meansq_list = []
    for a in arr:
        diff = a - mean
        diffsq = diff * diff
        diff2meansq_list.append(diffsq)
    std = ( 1.0/(N) * sum(diff2meansq_list) )**(0.5)
    #print ('sigma = ', std, '    mean = ', mean)
    return std, mean


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


def get_raw_star_directory(path4starfiles, scene, shutters, noise, redo=True):
    """
    This function returns a list of the directories (positions 1 and 2) to be studied.
    # Paths to Scenes 1 and 2 local directories: /Users/pena/Documents/AptanaStudio3/NIRSpec/TargetAcquisition/PFforMaria
    path_scene1_slow = "Scene_1_AB23/NIRSpec_TA_Sim_AB23 first NRS/postage"
    path_scene1_slow_nonoise = "Scene_1_AB23/NIRSpec_TA_Sim_AB23 first NRS no_noise/postage"
    path_scene1_rapid = "Scene_1_AB23/NIRSpec_TA_Sim_AB23 first NRSRAPID/postage"
    path_scene1_rapid_nonoise = "Scene_1_AB23/NIRSpec_TA_Sim_AB23 first NRS no_noise/postage"
    path_scene1_slow_shifted = "Scene_1_AB23/NIRSpec_TA_Sim_AB23 shifted NRS/postage"
    path_scene1_slow_shifted_nonoise = "Scene_1_AB23/NIRSpec_TA_Sim_AB23 shifted NRS no_noise/postage"
    path_scene1_rapid_shifted = "Scene_1_AB23/NIRSpec_TA_Sim_AB23 shifted NRSRAPID/postage"
    path_scene1_rapid_shifted_nonoise = "Scene_1_AB23/NIRSpec_TA_Sim_AB23 shifted NRS no_noise/postage"
    path_scene2_slow = "Scene_2_AB1823/NIRSpec_TA_Sim_AB1823 first NRS/postage"
    path_scene2_slow_nonoise = "Scene_2_AB1823/NIRSpec_TA_Sim_AB1823 first NRS no_noise/postage"
    path_scene2_rapid = "Scene_2_AB1823/NIRSpec_TA_Sim_AB1823 first NRSRAPID/postage"
    path_scene2_rapid_nonoise = "Scene_2_AB1823/NIRSpec_TA_Sim_AB1823 first NRSRAPID no_noise/postage"
    path_scene2_slow_shifted = "Scene_2_AB1823/NIRSpec_TA_Sim_AB1823 shifted NRS/postage"
    path_scene2_slow_shifted_nonoise = "Scene_2_AB1823/NIRSpec_TA_Sim_AB1823 shifted NRS no_noise/postage"
    path_scene2_rapid_shifted = "Scene_2_AB1823/NIRSpec_TA_Sim_AB1823 shifted NRSRAPID/postage"
    path_scene2_rapid_shifted_nonoise = "Scene_2_AB1823/NIRSpec_TA_Sim_AB1823 shifted NRSRAPID no_noise/postage"
    """
    # define shutter velocity to be used
    shutter_vel = "NRS"   # for slow case
    if shutters == "rapid":
        shutter_vel = "NRSRAPID"
    # define noise level
    noise_level = " no_noise"
    if noise == "real":
        noise_level = ""
    # define directory path for scenario 1
    position1 = path4starfiles+"Scene_"+repr(scene)+"_AB23/NIRSpec_TA_Sim_AB23 first "+shutter_vel+noise_level+"/postage"
    position2 = path4starfiles+"Scene_"+repr(scene)+"_AB23/NIRSpec_TA_Sim_AB23 shifted "+shutter_vel+noise_level+"/postage"
    if scene == 2:        
        position1 = path4starfiles+"Scene_"+repr(scene)+"_AB1823/NIRSpec_TA_Sim_AB1823 first "+shutter_vel+noise_level+"/postage"
        position2 = path4starfiles+"Scene_"+repr(scene)+"_AB1823/NIRSpec_TA_Sim_AB1823 shifted "+shutter_vel+noise_level+"/postage"
    if redo:
        position1 = position1+"_redo"
        position2 = position2+"_redo"
    dir2test_list = [position1, position2]
    return dir2test_list    
    
    
def Pier_correction(detector, XandYarr):
    """
    KEYWORD ARGUMENTS:
        Pier_corr                  -- Perform average correction suggested by Pier: True or False
        
    OUTPUT:
    cb_centroid_list               -- Values corrected for Pier's values 
    """
    # Corrections for offsets in positions (see section 2.5 of Technical Notes in Documentation directory)
    offset_491 = (-0.086, -0.077)
    offset_492 = (0.086, 0.077)
    corrected_x = XandYarr[0]
    corrected_y = XandYarr[1]
    if detector == 491:
        corrected_x = XandYarr[0] + offset_491[0]
        corrected_y = XandYarr[1] + offset_491[1]
    elif detector == 492:
        corrected_x = XandYarr[0] + offset_492[0]
        corrected_y = XandYarr[1] + offset_492[1]
    corr_XandYarr = [corrected_x, corrected_y]
    return corr_XandYarr


def Nsigma_rejection(N, x, y, max_iterations=10):
    """ This function will reject any residuals that are not within N*sigma in EITHER coordinate. 
        Input: 
                 - x and y must be the arrays of the differences with respect to true values: True-Measured 
                 - N is the factor (integer or float) by which sigma will be multiplied
                 - max_iterations is the maximum integer allowed iterations 
        Output:
                 - sigma_x = the standard deviation of the new array x 
                 - mean_x  = the mean of the new array x 
                 - sigma_y = the standard deviation of the new array y
                 - mean_y  = the mean of the new array y
                 - x_new   = the new array x (with rejections) 
                 - y_new   = the new array y (with rejections) 
                 - niter   = the number of iterations to reach a convergence (no more rejections)
        Usage:
             import testing_finctions as tf
             sigma_x, mean_x, sigma_y, mean_y, x_new, y_new, niter = tf.Nsigma_rejection(N, x, y, max_iterations=10)
    """
    N = float(N)
    or_sigma_x, or_mean_x = find_std(x)
    or_sigma_y, or_mean_y = find_std(y)
    x_new = copy.deepcopy(x)
    y_new = copy.deepcopy(y)
    original_diffs = copy.deepcopy(x)

    for nit in range(max_iterations):
        # Determine the standard deviation for each array
        sigma_x, mean_x = find_std(x_new)
        sigma_y, mean_y = find_std(y_new)
        thres_x = N*sigma_x
        thres_y = N*sigma_y
        xdiff = np.abs(x_new - mean_x) 
        ydiff = np.abs(y_new - mean_y)
        xn = x_new[(np.where((xdiff<=thres_x) & (ydiff<=thres_y)))]
        yn = y_new[(np.where((xdiff<=thres_x) & (ydiff<=thres_y)))]
        if len(xn) == len(x_new): 
            niter = nit
            break   # exit the loop since no additional rejections on this iteration
        else:
            x_new, y_new = xn, yn
            niter = nit
    line0 = "N-sigma rejection function:  values calculated from differences"
    line1 = "                             - stopped at {} iterations".format(niter)
    line2 = "                             - arrays have {} elements left out of {} initial".format(len(x_new), len(x))
    line3 = "                             - original sigma and mean in X = {}  {}".format(or_sigma_x, or_mean_x)
    line4 = "                             - original sigma and mean in Y = {}  {}".format(or_sigma_y, or_mean_y)
    line5 = "                             - new sigma and mean in X = {}  {}".format(sigma_x, mean_x)
    line6 = "                             - new sigma and mean in Y {}  {}".format(sigma_y, mean_y)
    lines2print = [line0, line1, line2, line3, line4, line5, line6]
    print (line0)
    print (line1)
    print (line2)
    print (line3)
    print (line4)
    print (line5)
    print (line6)
    # find what elements got rejected            
    rejected_elements_idx = []
    for i, centroid in enumerate(original_diffs):
        if centroid not in x:
            rejected_elements_idx.append(i)
    return sigma_x, mean_x, sigma_y, mean_y, x_new, y_new, niter, lines2print, rejected_elements_idx


def read_listfile(list_file_name, detector=None, background_method=None):    
    """ This function reads the fits table that contains the flux and converts to magnitude for the 
    simulated stars. """
    listfiledata = fits.getdata(list_file_name)
    star_number, xpos, ypos, orient, factor = np.array([]), np.array([]), np.array([]), np.array([]), np.array([])  
    for row in listfiledata:
        #print ("row: ", row)
        star_number = np.append(star_number, row[0]) 
        xpos = np.append(xpos, row[1]) 
        ypos = np.append(ypos, row[2])
        orient = np.append(orient, row[3])
        factor = np.append(factor, row[4])
    # convert the flux into magnitude (factor=1.0 is equivalent to magnitude=23.0,
    #  and factor=100.0 is equivalent to magnitude=18.0)
    mag = -2.5*np.log10(factor) + 23.0
    #mag = 2.5*np.log10(factor) + 18.0   # --> this conversion is wrong!
    # Get the correct slices according to detector or return all if no detector was chosen
    if detector is not None:
        if detector == 491:   # slice from star 101 to 200
            star_number, xpos, ypos, factor, mag = star_number[100:], xpos[100:], ypos[100:], factor[100:], mag[100:]
        elif detector == 492:   # slice from star 1 to 100
            star_number, xpos, ypos, factor, mag = star_number[:100], xpos[:100], ypos[:100:], factor[:100], mag[:100]
    bg_method = background_method
    if background_method is None:   # convert the None value to string
        bg_method = 'None'
    return star_number, xpos, ypos, factor, mag, bg_method


def read_positionsfile(positions_file_name, detector=None):
    """ This function reads the fits table that contains the true full detector positions of all simulated stars. """
    posfiledata = fits.getdata(positions_file_name)
    star_number, xpos491, ypos491, xpos492, ypos492 = np.array([]), np.array([]), np.array([]), np.array([]), np.array([])  
    trueV2, trueV3 = np.array([]), np.array([])
    for row in posfiledata:
        star_number = np.append(star_number, row[0]) 
        xpos491 = np.append(xpos491, row[3]) 
        ypos491 = np.append(ypos491, row[4])
        xpos492 = np.append(xpos492, row[5]) 
        ypos492 = np.append(ypos492, row[6])
        trueV2 = np.append(trueV2, row[13])
        trueV3 = np.append(trueV3, row[14])
    # Get the correct slices according to detector or return all if no detector was chosen
    xpos, ypos = np.array([]), np.array([])
    if detector is not None:
        if detector == 491:   # slice from star 101 to 200
            star_number, xpos, ypos, trueV2, trueV3 = star_number[100:], xpos491[100:], ypos491[100:], trueV2[100:], trueV3[100:]
        elif detector == 492:   # slice from star 1 to 100
            star_number, xpos, ypos, trueV2, trueV3 = star_number[:100], xpos492[:100], ypos492[:100], trueV2[:100], trueV3[:100]
    else:
        # return the 200 values
        xpos = np.append(xpos, xpos492[:100])
        xpos = np.append(xpos, xpos491[100:])
        ypos = np.append(ypos, ypos492[:100])
        ypos = np.append(ypos, ypos491[100:])
    #for s, x, y in zip(star_number, xpos, ypos):
    #    print(s, x, y)
    #    raw_input()
    return star_number, xpos, ypos, trueV2, trueV3


def run_recursive_centroids(psf, background, xwidth_list, ywidth_list, checkbox_size, max_iter, 
                            threshold, determine_moments, debug):   
    """
    Determine the centroid location given the that the background is already subtracted. 
    """ 
    # Test checkbox piece
    print ("Centroid measurement for background of: ", background)
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


def read_star_param_files(test_case, detector=None, path4starfiles=None, paths_list=None):
    """ This function reads the corresponding star parameters file and returns the data for P1 and P2. """
    cases_list = ["Scene1_slow_real", "Scene1_slow_nonoise", "Scene1_rapid_real", "Scene1_rapid_nonoise",
                  "Scene2_slow_real", "Scene2_slow_nonoise", "Scene2_rapid_real", "Scene2_rapid_nonoise"]
    if path4starfiles is not None and paths_list is not None:
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
            
        """
        *** WE ARE NOT USING THIS PART RIGHT NOW BECAUSE THE star_parameters FILES HAVE THE SAME DATA FOR 
        BOTH DETECTORS.
        #    xL:  x-coordinate of the left edge of the postge stamp in the full image (range 0-2047)
        #    xR: x-coord of right edge of the postage stamp
        #    yL: y-coord of the lower edge of the postage stamp
        #    yU:  y-coord of the upper edge of the postage stamp
        
        # Load parameters of Position 1
        star_param_txt = os.path.join(dirs4test[0],"star parameters.txt")
        if detector == 492:
            star_param_txt = os.path.join(dirs4test[0],"star parameters_492.txt")
        benchmark_dataP1 = np.loadtxt(star_param_txt, skiprows=3, unpack=True)
        # benchmark_data is list of: bench_star, quadrant, star_in_quad, x_491, y_491, x_492, y_492, V2, V3, xL, xR, yL, yU 
        bench_starP1, _, _, x_491P1, y_491P1, x_492P1, y_492P1, V2P1, V3P1, xLP1, _, yLP1, _ = benchmark_dataP1
        if detector == 491:
            bench_P1 = [bench_starP1, x_491P1, y_491P1, V2P1, V3P1, xLP1, yLP1]
        elif detector == 492:
            bench_P1 = [bench_starP1, x_492P1, y_492P1, V2P1, V3P1, xLP1, yLP1]        
        # Load parameters of Position 2
        star_param_txt = os.path.join(dirs4test[1],"star parameters.txt")
        if detector == 492:
            star_param_txt = os.path.join(dirs4test[1],"star parameters_492.txt")
        benchmark_dataP2 = np.loadtxt(star_param_txt, skiprows=3, unpack=True)
        #bench_star, quadrant, star_in_quad, x_491, y_491, x_492, y_492, V2, V3, xL, xR, yL, yU = benchmark_data
        bench_starP2, _, _, x_491P2, y_491P2, x_492P2, y_492P2, V2P2, V3P2, xLP2, _, yLP2, _ = benchmark_dataP2
        if detector == 491:
            bench_P2 = [bench_starP2, x_491P2, y_491P2, V2P2, V3P2, xLP2, yLP2]
        elif detector == 492:
            bench_P2 = [bench_starP2, x_492P2, y_492P2, V2P2, V3P2, xLP2, yLP2]        
        """
    #else:
    # Read fits table with benchmark data
    main_path_infiles = "../PFforMaria/"
    S1path2listfile = main_path_infiles+"Scene_1_AB23"
    S1list_file1 = "simuTA20150528-F140X-S50-K-AB23.list"
    S1positions_file1 = "simuTA20150528-F140X-S50-K-AB23_positions.fits" 
    S1list_file2 = "simuTA20150528-F140X-S50-K-AB23-shifted.list"
    S1positions_file2 = "simuTA20150528-F140X-S50-K-AB23-shifted_positions.fits"
    S2path2listfile = main_path_infiles+"Scene_2_AB1823"
    S2list_file1 = "simuTA20150528-F140X-S50-K-AB18to23.list"
    S2positions_file1 = "simuTA20150528-F140X-S50-K-AB18to23_positions.fits"
    S2list_file2 = "simuTA20150528-F140X-S50-K-AB18to23-shifted.list"
    S2positions_file2 = "simuTA20150528-F140X-S50-K-AB18to23-shifted_positions.fits"
    if "Scene1" in test_case:
        benchmark_data, magnitudes = read_TruePosFromFits(S1path2listfile, S1list_file1, S1positions_file1, S1list_file2, S1positions_file2)
    if "Scene2" in test_case:
        benchmark_data, magnitudes = read_TruePosFromFits(S2path2listfile, S2list_file1, S2positions_file1, S2list_file2, S2positions_file2)
    return benchmark_data, magnitudes


def read_TruePosFromFits(path2listfile, list_file1, positions_file1, list_file2, positions_file2, test_case=None, detector=None):
    # Read the text file just written to get the offsets from the "real" positions of the fake stars
    lf1 = os.path.join(path2listfile, list_file1)
    pf1 = os.path.join(path2listfile, positions_file1)
    lf2 = os.path.join(path2listfile, list_file2)
    pf2 = os.path.join(path2listfile, positions_file2)
    if test_case is not None:
        if "None" in test_case:
            background_method = None
        elif "fix" in test_case:
            background_method = "fix"
        elif "frac" in test_case:
            background_method = "frac"
    else:
        background_method = None
    bench_starP1, xpos_arcsecP1, ypos_arcsecP1, factorP1, magP1, bg_methodP1 = read_listfile(lf1, detector, background_method)
    _, true_xP1, true_yP1, trueV2P1, trueV3P1 = read_positionsfile(pf1, detector)
    bench_starP2, xpos_arcsecP2, ypos_arcsecP2, factorP2, magP2, bg_methodP2 = read_listfile(lf2, detector, background_method)
    _, true_xP2, true_yP2, trueV2P2, trueV3P2 = read_positionsfile(pf2, detector)
    # Get the lower left corner coordinates in terms of full detector. We subtract 15.0 because indexing
    # starts with 0
    xLP1 = np.floor(true_xP1) - 16.0
    yLP1 = np.floor(true_yP1) - 16.0
    xLP2 = np.floor(true_xP2) - 16.0
    yLP2 = np.floor(true_yP2) - 16.0
    # Organize elements of positions 1 and 2
    bench_P1 = [bench_starP1, true_xP1, true_yP1, trueV2P1, trueV3P1, xLP1, yLP1]
    bench_P2 = [bench_starP2, true_xP2, true_yP2, trueV2P2, trueV3P2, xLP2, yLP2]
    benchmark_data = [bench_P1, bench_P2]
    return benchmark_data, magP1


# Extract an image from a multi-ramp integration FITS file
def readimage(master_img, backgnd_subtraction_method=None, bg_method=None, bg_value=None, bg_frac=None, debug=False):
    """
    Extract am image from a multi-ramp integration FITS file.
    Currently, JWST NIRSpec FITS images consists of a 3-ramp integration, 
    with each succesive image containing more photon counts than the next. 
    Uses a cube-differencing calculation to eliminate random measurements
    such as cosmic rays.
    Keyword arguments:
    master_img                 -- 3-frame image (as per NIRSpec output images)
    backgnd_subtraction_method -- 1 = Do background subtraction on final image (after subtracting 3-2 and 2-1), 
                                      before converting negative values into zeros
                                  2 = Do background subtraction on 3-2 and 2-1 individually
    Output(s):
    omega -- A combined FITS image that combines all frames into one image.
    """
    
    # Read in input file, and generate the alpha and beta images
    # (for subtraction)
    #alpha = master_img[1, :, :] - master_img[0, :, :]
    #beta = master_img[2, :, :] - master_img[1, :, :]
    alpha = master_img[1] - master_img[0]
    beta = master_img[2] - master_img[1]
    #fits.writeto("/Users/pena/Documents/AptanaStudio3/NIRSpec/TargetAcquisition/alpha.fits", alpha)
    #fits.writeto("/Users/pena/Documents/AptanaStudio3/NIRSpec/TargetAcquisition/beta.fits", beta)
    
    # Perform background subtraction if backgnd_subtraction_method=1
    if backgnd_subtraction_method == 2:
        print ("*  Background subtraction being done on 3-2 and 2-1 individually...")
        alpha = bg_correction(alpha, bg_method=bg_method, bg_value=bg_value, bg_frac=bg_frac, debug=debug)
        beta = bg_correction(beta, bg_method=bg_method, bg_value=bg_value, bg_frac=bg_frac, debug=debug)
    
    # Generate a final image by doing a pixel to pixel check 
    # between alpha and beta images, storing lower value
    omega = np.where(alpha < beta, alpha, beta)
    
    # Perform background subtraction if backgnd_subtraction_method=1
    if backgnd_subtraction_method == 1:
        print ("*  Background subtraction being done on image of min between 3-2 and 2-1...")
        omega = bg_correction(omega, bg_method=bg_method, bg_value=bg_value, bg_frac=bg_frac, debug=debug)

    # Convert negative pixel values to zero
    negative_idx = np.where(omega < 0.0)
    omega[negative_idx] = 0.0

    # show on screen the values of rows and column for debugging other functions of the code
    if debug:
        image = omega   # Type image to be displayed
        if image is omega:
            print ('Combined ramped images:  ')
            print ('   AFTER zeroing negative pixels')
        else:
            print ('   BEFORE zeroing negative pixels')
        print ('max_image = ', image.max())
        print ('min_image = ', image.min())
        for j in range(np.shape(image)[0]):
            print (j, image[j, :])#, alpha[j, :], beta[j, :])
    

    print('(readimage): Image processed!')
        
    # Return the extracted image
    return omega

