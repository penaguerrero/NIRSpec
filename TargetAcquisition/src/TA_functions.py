from __future__ import print_function, division
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
import copy
import os
import collections
import string

# Tommy's code
import jwst_targloc as jtl

# other code
import coords_transform as ct
import least_squares_iterate as lsi

# Header
__author__ = "Maria A. Pena-Guerrero"
__version__ = "1.0"

"""

This script has Target Acquisition functions that are auxiliary to the TA functions in the jwst_targloc.py
script.

* Functions are ordered alphabetically.

** Functions specific to test if the averaging of centroids in pixel space, sky, or individually returns the best
results are located at the end of this file.

"""

# FUNCTIONS AUXILIARY TO JWST CENTROID ALGORITHM

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

    Args:
        img        -- Image
        bg_method  -- Either None value or string: "fixed" or "frac"
        bg_value   -- Float, fixed value to subtract from each pixel (this has to
                      be set if bg_method = "fixed")
        bg_frac    -- Float, fractional value to subtract from image (this has to
                      be set if bg_method = "frac")
    
    Returns:
        img_bgcorr -- The group of 3 background subtracted images

    Example usage:
        >> img_bgcorr = bg_correction(master_img, bg_method='frac', bg_value=0.4)
    """
    # Make sure to return the image as is if None is selected
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
        # Find the pixel value (bg) that represents that fraction of the population
        img_original = copy.deepcopy(img)
        sorted_img = np.sort(np.ravel(img))   # flatten the image and sort it
        xsize = np.shape(img)[1]
        ysize = np.shape(img)[0]
        idx_bg = np.floor(bg_frac * xsize * ysize)
        # If at the edge, correct
        if idx_bg == np.shape(sorted_img)[0]:
            idx_bg -= 1
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
    
    Args:
        cb_centroid_list           -- List, centroid window based centroid determined by TA algorithm in
                                      terms of 32 by 32 pixels for centroid window sizes 3, 5, and 7
        true_center                -- List, actual (true) position of star in terms of full detector
    
    Returns:
        cb_centroid_list_fulldetector  -- List of centroid locations determined with the TA algorithm in
                                          terms of full detector. List is for positions determined with
                                          3, 5, and 7 centroid window sizes.
        loleftcoords                   -- List, Coordinates of the lower left corner of the 32x32 pixel box
        true_center32x32               -- List, true center given in coordinates of 32x32 pix
        differences_true_TA            -- List, difference of true-observed positions
    """
        
    # Get the lower left corner coordinates in terms of full detector. We subtract 16.0 because indexing
    # from centroid function starts with 1
    corrected_x = true_center[0]
    corrected_y = true_center[1]
    loleft_x = np.floor(corrected_x) - 16.0
    loleft_y = np.floor(corrected_y) - 16.0
    loleftcoords = [loleft_x, loleft_y]

    # get center in terms of 32x32 cutout
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
    if len(corr_cb_centroid_list) != 1:   # make sure this function works even for one centroid window
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
    """
    This function obtains the differences of the input arrays with the reference or benchmark data.
    Args:
        case        -- string, for example 'Scene2_rapid_real_bgFrac'
        bench_stars -- numpy array of the star numbers being used
        benchV2     -- numpy array of the benchmark V2s
        benchV3     -- numpy array of the benchmark V3s
        stars       -- list of the star numbers being studied
        V2in        -- numpy array of measured V2s
        V3in        -- numpy array of measured V3s
        arcsecs     -- True or False (if False results are in degrees)

    Returns:
        4 lists: diffV2, diffV3, bench_V2_list, bench_V3_list
        diffV2 = benchmark V2 - measured V2
        diffV3 = benchmark V3 - measured V3
        bench_V2_list = benchmark V2 converted to arcsecs (if arcsecs=True)
        bench_V3_list = benchmark V3 converted to arcsecs (if arcsecs=True)
    """

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


def convert2MSAcenter(xin, yin, xtin, ytin):
    """
    This function is a python translation of Tonys IDL equivalent function. It converts the
    measured coordinates of each star into the frame relative to the center of the MSA.
    Args:
        xin: numpy array of the measured V2s
        yin: numpy array of the measured V2s
        xtin: numpy array of true V2s
        ytin: numpy array of true V3s

    Returns:
        4 numpy arrays: x, y, xt, yt - Corrected arrays
    """
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


def display_centroids(detector, st, case, psf, corr_true_center_centroid,
                      corr_cb_centroid_list, show_disp, vlims=None, savefile=False, 
                      fig_name=None, redos=False, display_master_img=False):
    """
    This function displays de the centroids for the 32x32 pixel cutout images, showing
    the true position as wel as the measured centroids.
    Args:
        detector                  -- integer, either 491 or 492
        st                        -- integer, star number
        case                      -- string, for example 'Scene2_rapid_real_bgFrac'
        psf                       --  numpy array of shape (3, 32, 32) -- cutout of 3 ramped images
        corr_true_center_centroid -- list of x and y true pixel positions
        corr_cb_centroid_list     -- list of 3 lists, x and y pixel positions for centroiding window 3, 5, and 7
        show_disp                 -- True or False, show the 32x32 image with true center and measured positions
        vlims                     -- tuple, example: (10.0, 50.0)
        savefile                  -- True or False, save or not the image as a .jpg
        fig_name                  -- string, name for the figure
        redos                     -- True or False, use or not the directories with _redo
        display_master_img        -- True or False, show the initial image (before background subtraction)

    Returns:
        Nothing.
    """
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
    ax.plot(corr_cb_centroid_list[0][0], corr_cb_centroid_list[0][1], marker='*', ms=20, mec='black', mfc='blue', ls='', label='CentroidWin=3')
    if len(corr_cb_centroid_list) != 1:
        ax.plot(corr_cb_centroid_list[1][0], corr_cb_centroid_list[1][1], marker='*', ms=17, mec='black', mfc='green', ls='', label='CentroidWin=5')
        ax.plot(corr_cb_centroid_list[2][0], corr_cb_centroid_list[2][1], marker='*', ms=15, mec='black', mfc='red', ls='', label='CentroidWin=7')
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
    
    
def display_ns_psf(image, vlim=(), fsize=(8, 8), interp='nearest', title='',
                   cmap='gray', extent=None, savefile=None, cb=False):
    """
    Custom display a PSF generated with WEBBPSF or similar tool.
    A quick tool for displaying NIRSpec images in native size 
    (2048x2048) with additional options for display.
    Args:
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
    Returns:
        Nothing.
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

    # See plots when not in Notebook environment
    plt.show()

    if savefile is not None:
        fig.savefig(savefile)


def do_Piers_correction(detector, cb_centroid_list):
    """ This function performs the average correction found by Pierre in each of
        the x and y measured centroids.
    Args:
        detector          -- integer, either 491 or 492
        cb_centroid_list  -- list of 3 lists (measured centroids for centroid window sizes 3, 5, and 7)

    Returns:
        List of 3 lists corresponding to corrected centroids for centroid window sizes 3, 5, and 7
    """
    xy3, xy5, xy7 = cb_centroid_list
    xy3corr = Pier_correction(detector, xy3)
    xy5corr = Pier_correction(detector, xy5)
    xy7corr = Pier_correction(detector, xy7)
    corr_cb_centroid_list = [xy3corr, xy5corr, xy7corr]
    return corr_cb_centroid_list 
    
    
def find_centroid(fits_file, bg_corr_info, recursive_centroids_info, display_centroids_info, x_centroids, y_centroids,
                  fits_names, output_file_path, centroids_info):
    """ This function reads the image, finds the centroid, and displays the result.
    It returns the centroid values.

    Args:
        fits_file                -- name of the fits file being studied
        bg_corr_info             -- list of the information concerning background subtraction
        recursive_centroids_info -- list of information for running centroid algorithm
        display_centroids_info   -- list of information to show and display the centroids
        x_centroids              -- list of centroids for centroid window sizes of 3, 5, and 7 for x position
        y_centroids              -- list of centroids for centroid window sizes of 3, 5, and 7 for y position
        fits_names               -- list to the append the studied files (so that it ends up being the same
                                    length as the list of the measured centroids -- in case of fractional background)
        output_file_path         -- path for the output file
        centroids_info           -- list of the information concerning the true centroids, the output in full
                                    detector coordinates, and the on-screen measured centroids
    Returns:
        x_centroids = list of 3 lists corresponding to pixel x-positions for centroid window sizes 3, 5, and 7
        y_centroids = list of 3 lists corresponding to pixel y-positions for centroid window sizes 3, 5, and 7
    """

    # unfold information
    backgnd_subtraction_method, background_method, bg_value, bg_frac, debug = bg_corr_info
    xwidth_list, ywidth_list, centroid_win_size, max_iter, threshold, determine_moments, display_master_img, vlim = recursive_centroids_info
    true_center, output_full_detector, show_centroids = centroids_info
    case, show_disp, save_centroid_disp = display_centroids_info
    x_centroids3, y_centroids3 = x_centroids[0], y_centroids[0]
    x_centroids5, y_centroids5 = x_centroids[1], y_centroids[1]
    x_centroids7, y_centroids7 = x_centroids[2], y_centroids[2]

    # get detector and name of base name of each fits file
    ff = os.path.basename(fits_file)
    ff1 = string.split(ff, sep="_")
    detector = ff1[2]
    #fits_trial = ff1[1]
    #fits_base = ff1[0]

    # Read FITS image
    #img = fits.open(fits_file)
    #img.info()
    #raw_input()
    #hdr = fits.getheader(fits_file, 0)
    #print("** HEADER:", hdr)
    master_img = fits.getdata(fits_file, 0)
    print ('Master image shape: ', np.shape(master_img))
    # Obtain the combined FITS image that combines all frames into one image AND
    # check if all image is zeros, take the image that still has a max value
    psf = readimage(master_img, backgnd_subtraction_method, bg_method=background_method,
                        bg_value=bg_value, bg_frac=bg_frac, debug=debug)
    cb_centroid_list_in32x32pix = run_recursive_centroids(psf, bg_frac, xwidth_list, ywidth_list,
                                               centroid_win_size, max_iter, threshold,
                                               determine_moments, debug)
    cb_centroid_list, loleftcoords, true_center32x32, differences_true_TA = centroid2fulldetector(cb_centroid_list_in32x32pix,
                                                                                        true_center)
    if not output_full_detector:
        cb_centroid_list = cb_centroid_list_in32x32pix
    if show_centroids:
        print ('***** Measured centroids:')
        print ('      cb_centroid_list = ', cb_centroid_list)
        #print ('           True center = ', true_center)

    x_centroids3.append(cb_centroid_list[0][0])
    y_centroids3.append(cb_centroid_list[0][1])
    if len(xwidth_list) != 1:
        x_centroids5.append(cb_centroid_list[1][0])
        y_centroids5.append(cb_centroid_list[1][1])
        x_centroids7.append(cb_centroid_list[2][0])
        y_centroids7.append(cb_centroid_list[2][1])
    x_centroids = [x_centroids3, x_centroids5, x_centroids7]
    y_centroids = [y_centroids3, y_centroids5, y_centroids7]

    # Show the display with the measured and true positions
    ff = string.replace(ff, ".fits", "")
    fits_names.append(ff)
    fig_name = os.path.join(output_file_path, ff+".jpg")
    # Display the combined FITS image that combines all frames into one image
    m_img = display_master_img
    if display_master_img:
        m_img = readimage(master_img, backgnd_subtraction_method=None, bg_method=None,
                          bg_value=None, bg_frac=None, debug=False)
    display_centroids(detector, ff, case, psf, true_center32x32, cb_centroid_list_in32x32pix,
                         show_disp, vlim, savefile=save_centroid_disp, fig_name=fig_name, display_master_img=m_img)
    return x_centroids, y_centroids


def find_std(arr):
    """
    This function determines the standard deviation of the given array.
    Args:
        arr = numpy array for which the standard deviation and means are to be determined

    Returns:
        std = standard deviation of the given array
        mean = mean value of the given array

    Usage:
         import TA_functions as taf
         std, mean = taf.find_std(y_positions)
    """
    N = float(len(arr))
    mean = sum(arr) / N
    diff2meansq_list = []
    for a in arr:
        diff = a - mean
        diffsq = diff * diff
        diff2meansq_list.append(diffsq)
    std = (1.0 / N * sum(diff2meansq_list)) ** 0.5
    #print ('sigma = ', std, '    mean = ', mean)
    return std, mean


def get_mindiff(d1, d2, d3):
    """
    This function determines the minimum difference from centroid window sizes 3, 5, and 7,
    and counts the number of repetitions.
    Args:
        d1: list of differences of measured values with respect to true for centroid window size 3
        d2: list of differences of measured values with respect to true for centroid window size 5
        d3: list of differences of measured values with respect to true for centroid window size 7

    Returns:
        min_diff = the minimum difference centroid window
        counter = a dictionary that has the centroid window sizes and their repetitions in order of most
                    repeated to least, example: {{7: 13}, {5: 8}, {3: 2}}
    """
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
    Possible paths to Scenes 1 and 2 directories in directory PFforMaria:
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
    Args:
        path4starfiles -- string, path to get to the files
        scene          -- integer, either 1 or 2
        shutters       -- string, shutter velocity: 'rapid' or 'slow'
        noise          -- string, noise level: 'nonoise' or 'real'
        redo           -- True or False, go to (or not) to the directories that have a _redo at the end
    Returns:
        dir2test_list = A list of strings with the paths to position files 1 and 2
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
        position1 += "_redo"
        position2 += "_redo"
    dir2test_list = [position1, position2]
    return dir2test_list    
    
    
def Pier_correction(detector, XandYarr):
    """ This function corrects the measured centroids for the average values.
    Args:
        Pier_corr         -- Perform average correction suggested by Pier: True or False
        
    Returns:
        cb_centroid_list  -- List, values corrected for Pier's values
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


def print_file_lines(output_file, save_text_file, xwidth_list, ff, background2use,
                     i, x_centroids, y_centroids):
    """ This function prints the info on screen AND in a text file. It expects that the output file
    already exists (it appends information to it). Columns are fits file name, background value used
    and method (None, fixed, or fractional), x and y coordinates for centroid window sizes of 3, 5, and 7
    pixels. Eight columns in total.
    Args:
        output_file    -- name of the output file = string
        save_text_file -- do you want to save the file = True or False
        xwidth_list    -- list of x width size = list can have 1 to 3 elements
        ff             -- name of the fits file the info corresponds to = string
        background2use -- value of the background to use for the background method = float
        i              -- the index of line to append to the file = integer
        x_centroids    -- list of the x measured centroids for centroid window sizes 3, 5, and 7
        y_centroids    -- list of the y measured centroids for centroid window sizes 3, 5, and 7
    Returns:
        Nothing.
    """
    x_centroids3, y_centroids3 = x_centroids[0], y_centroids[0]
    x_centroids5, y_centroids5 = x_centroids[1], y_centroids[1]
    x_centroids7, y_centroids7 = x_centroids[2], y_centroids[2]
    if len(xwidth_list)==1:
        line1 = "{:<40} {:>4} {:>16} {:>16}".format(
                                   ff, background2use, x_centroids3[i], y_centroids3[i])
    else:
        line1 = "{:<40} {:>4} {:>16} {:>14} {:>18} {:>14} {:>18} {:>14}".format(
                                                               ff, background2use,
                                                               x_centroids3[i], y_centroids3[i],
                                                               x_centroids5[i], y_centroids5[i],
                                                               x_centroids7[i], y_centroids7[i])
    print (line1)
    if save_text_file:
        f = open(output_file, "a")
        f.write(line1+"\n")


def Nsigma_rejection(N, x, y, max_iterations=10):
    """ This function will reject any residuals that are not within N*sigma in EITHER coordinate. 
    Args:
         - x and y must be the numpy arrays of the differences with respect to true values: True-Measured
         - N is the factor (integer or float) by which sigma will be multiplied
         - max_iterations is the maximum integer allowed iterations
    Returns:
         - sigma_x = the standard deviation of the new array x
         - mean_x  = the mean of the new array x
         - sigma_y = the standard deviation of the new array y
         - mean_y  = the mean of the new array y
         - x_new   = the new array x (with rejections)
         - y_new   = the new array y (with rejections)
         - niter   = the number of iterations to reach a convergence (no more rejections)
    Usage:
         import TA_functions as taf
         sigma_x, mean_x, sigma_y, mean_y, x_new, y_new, niter = taf.Nsigma_rejection(N, x, y, max_iterations=10)
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
    """
    This function reads the fits table that contains the flux and converts to magnitude for the
    simulated stars.
    Args:
        list_file_name    -- string, fits file of the list file to be used
        detector          -- integer, either 491 or 492
        background_method -- None or string, either 'fractional' or 'fixed'

    Returns:
        5 numpy arrays and a string: star_number, xpos, ypos, factor, mag, bg_method
        star_number = number of all stars of the selected detector
        xpos = x pixel position of star_number
        ypos = y pixel position of star_number
        factor = simulated variance of a star of magnitude 23 (factor=1)
        mag = conversion of factor into magnitudes
        bg_method = string, method used for background subtraction 'None', 'fractional', 'fixed'
    """
    listfiledata = fits.getdata(list_file_name)
    star_number, xpos, ypos, orient, factor = np.array([]), np.array([]), np.array([]), np.array([]), np.array([])  
    for row in listfiledata:
        star_number = np.append(star_number, row[0])
        xpos = np.append(xpos, row[1]) 
        ypos = np.append(ypos, row[2])
        orient = np.append(orient, row[3])
        factor = np.append(factor, row[4])
    # convert the flux into magnitude (factor=1.0 is equivalent to magnitude=23.0,
    #  and factor=100.0 is equivalent to magnitude=18.0)
    mag = -2.5*np.log10(factor) + 23.0
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
    """
    This function reads the fits table that contains the true full detector positions of all simulated stars.
    Args:
        positions_file_name -- string, name of the fits file to be read
        detector            -- integer, 491 or 492

    Returns:
        5 numpy arrays: star_number, xpos, ypos, trueV2, trueV3
        star_number = number of all stars of the selected detector
        xpos = x pixel position of star_number
        ypos = y pixel position of star_number
        trueV2 = V2s of those star_number
        trueV3 = V3s of those star_number
    """

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
    return star_number, xpos, ypos, trueV2, trueV3


def run_recursive_centroids(psf, background, xwidth_list, ywidth_list, checkbox_size, max_iter,
                            threshold, determine_moments, debug):   
    """
    This function determines the centroid given that the background is already subtracted.
    Args:
        psf: numpy array of shape (3, 32, 32) -- cutout of 3 ramped images
        background: float, fixed or fractional value used for background subtraction
        xwidth_list: list of integer numbers -- the size of centroiding window in the x-direction
        ywidth_list: list of integer numbers -- the size of centroiding window in the y-direction
        checkbox_size: integer, the size of the checkbox to use -- usually 3
        max_iter: integer, maximum iteration number for centroiding function
        threshold: float, convergence threshold of accepted difference between checkbox centroid and coarse location
        determine_moments: True or False
        debug: True or False

    Returns:
        cb_centroid_list = list of 3 lists with the x and y pixel positions for centroid window sizes 3, 5, and 7.
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


def read_star_param_files(test_case, path4starfiles=None, paths_list=None):
    """
    This function reads the corresponding star parameters file and returns the data for P1 and P2.
    Args:
        test_case     -- string, for example 'Scene2_rapid_real_bgFrac'
        path4starfile -- string, path to find the fits file that is in the postage stamps directory
        list_file1    -- string, name of the file to get magnitudes for position 1

    Returns:
        2 lists: benchmark_data, magnitudes
        benchmark_data = list of benchmark positions 1 and 2
        magnitudes = list of magnitudes (same for both positions)
    """

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
            
    # *** THE star_parameters.txt FILES in the postage stamps directory HAVE THE SAME DATA FOR BOTH DETECTORS.

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
    """
    This function gets the 'true' positions from the fits files Pierre sent.
    Args:
        path2listfile: string, path to find the fits file
        list_file1: string, name of the file to get magnitudes for position 1
        positions_file1: string, name of file to get 1st positions
        list_file2: string, name of the file to get magnitudes for position 2
        positions_file2: string, name of file to get 2nd positions
        test_case: string, for example 'Scene2_rapid_real_bgFrac'
        detector: integer, either 491 or 492

    Returns:
        2 lists: benchmark_data, magP1
        benchmark_data = list of benchmark positions 1 and 2
        magP1 = list of magnitudes (same for both positions)
    """
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
    Args:
        master_img                 -- 3-framed image (as per NIRSpec output images)
        backgnd_subtraction_method -- 1 = Do background subtraction on final image (after subtracting 3-2 and 2-1),
                                          before converting negative values into zeros
                                      2 = Do background subtraction on 3-2 and 2-1 individually
    Returns:
        omega -- A combined FITS image that combines all frames into one image.
    """
    
    # Read in input file, and generate the alpha and beta images
    # (for subtraction)
    alpha = master_img[1] - master_img[0]
    beta = master_img[2] - master_img[1]

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
            print (j, image[j, :])
    

    print('(readimage): Image processed!')
        
    # Return the extracted image
    return omega


#######################################################################################################################


"""
*Functions specific to test if the averaging of centroids in pixel space, sky, or individually returns the best results:
     TEST1 - Average positions P1 and P2, transform to V2-V3 space, and compare to average
             reference positions (V2-V3 space)
     TEST2 - Transform individual positions P1 and P2 to V2-V3 space, average V2-V3 space
             positions, and compare to average reference positions.
     TEST3 - Transform P1 and P2 individually to V2-V3 space and compare star by star and
             position by position.
"""


def TEST1(detector, transf_direction, stars, case, bench_starP1, avg_benchV23, P1P2data,
          filter_input, tilt, diffs_in_arcsecs, debug):
    """
    TEST 1: (a) Avg P1 and P2
            (b) transform to V2-V3
            (c) compare to avg reference positions (V2-V3 space)
    Args:
        detector: integer, either 491 or 492
        transf_direction: string, 'forward' or 'backward'
        stars: list, star sample to be studied
        case: string, for example 'Scene2_rapid_real_bgFrac'
        bench_starP1: list, benchmark star positions
        bench_Vs: list of 2 lists, benchmark V2s and V3s
        P1P2data: list of 6 lists with x and y positions 1 and 2 for centroid window sizes 3, 5, and 7
        filter_input: string, exapmple 'F140X'
        tilt: True or False
        diffs_in_arcsecs: True or False (if False result is degrees)
        debug: True or False (if True functions are very verbose)

    Returns:
        3 lists: T1_transformations, T1_diffs, T1_benchVs_list
        T1_transformations = list of 6 lists with V2s and V3s for centroid window sizes 3, 5, and 7
        T1_diffs = list of 6 lists with differences of V2 and V3 measured with respect to benchmark
        T1_benchVs_list = list of 2 lists of averaged V2s and V3s
    """

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
    # Step (c) - comparison
    T1_diffV2_3, T1_diffV3_3, T1bench_V2_list, T1bench_V3_list = compare2ref(case, bench_starP1, avg_benchV2, avg_benchV3, stars, T1_V2_3, T1_V3_3, arcsecs=diffs_in_arcsecs)
    T1_diffV2_5, T1_diffV3_5, _, _ = compare2ref(case, bench_starP1, avg_benchV2, avg_benchV3, stars, T1_V2_5, T1_V3_5, arcsecs=diffs_in_arcsecs)
    T1_diffV2_7, T1_diffV3_7, _, _ = compare2ref(case, bench_starP1, avg_benchV2, avg_benchV3, stars, T1_V2_7, T1_V3_7, arcsecs=diffs_in_arcsecs)
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

def TEST2(detector, transf_direction, stars, case, bench_starP1, avg_benchV23, P1P2data,
          filter_input, tilt, diffs_in_arcsecs, debug):
    """
    TEST 2: (a) Transform individual P1 and P2 to V2-V3
            (b) avg V2-V3 space positions
            (c) compare to avg reference positions

    Args:
        detector: integer, either 491 or 492
        transf_direction: string, 'forward' or 'backward'
        stars: list, star sample to be studied
        case: string, for example 'Scene2_rapid_real_bgFrac'
        bench_starP1: list, benchmark star positions
        bench_Vs: list of 2 lists, benchmark V2s and V3s
        P1P2data: list of 6 lists with x and y positions 1 and 2 for centroid window sizes 3, 5, and 7
        filter_input: string, exapmple 'F140X'
        tilt: True or False
        diffs_in_arcsecs: True or False (if False result is degrees)
        debug: True or False (if True functions are very verbose)

    Returns:
        3 lists: T2_transformations, T2_diffs, T2_benchVs_list
        T2_transformations = list of 6 lists with averaged V2s and V3s for centroid window sizes 3, 5, and 7
        T2_diffs = list of 6 lists with differences of V2 and V3 measured with respect to benchmark
        T2_benchVs_list = list of 2 lists of averaged V2s and V3s
    """

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
    T2_diffV2_3, T2_diffV3_3, T2bench_V2_list, T2bench_V3_list = compare2ref(case, bench_starP1, avg_benchV2, avg_benchV3, stars, T2_V2_3, T2_V3_3, arcsecs=diffs_in_arcsecs)
    T2_diffV2_5, T2_diffV3_5, _, _ = compare2ref(case, bench_starP1, avg_benchV2, avg_benchV3, stars, T2_V2_5, T2_V3_5, arcsecs=diffs_in_arcsecs)
    T2_diffV2_7, T2_diffV3_7, _, _ = compare2ref(case, bench_starP1, avg_benchV2, avg_benchV3, stars, T2_V2_7, T2_V3_7, arcsecs=diffs_in_arcsecs)
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


def TEST3(detector, transf_direction, stars, case, bench_starP1, bench_Vs, P1P2data,
          filter_input, tilt, diffs_in_arcsecs, debug):
    """
    TEST 3: (a) Transform P1 and P2 individually to V2-V3
            (b) compare star by star and position by position

    Args:
        detector: integer, either 491 or 492
        transf_direction: string, 'forward' or 'backward'
        stars: list, star sample to be studied
        case: string, for example 'Scene2_rapid_real_bgFrac'
        bench_starP1: list, benchmark star positions
        bench_Vs: list of 2 lists, benchmark V2s and V3s
        P1P2data: list of 6 lists with x and y positions 1 and 2 for centroid window sizes 3, 5, and 7
        filter_input: string, exapmple 'F140X'
        tilt: True or False
        diffs_in_arcsecs: True or False (if False result is degrees)
        debug: True or False (if True functions are very verbose)

    Returns:
        3 lists: T3_transformations, T3_diffs, T3_benchVs_list
        T3_transformations = list of 2 lists with V2s and V3s for positions 1 and 2
        T3_diffs = list of 2 lists with differences of V2 and V3 measured with respect to benchmark for
                    positions 1 and 2
        T3_benchVs_list = list of 4 lists of V2s and V3s for positions 1 and 2
    """
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
    T3_diffV2_13, T3_diffV3_13, T3bench_V2_listP1, T3bench_V3_listP1 = compare2ref(case, bench_starP1, bench_V2P1, bench_V3P1, stars, T3_V2_13, T3_V3_13, arcsecs=diffs_in_arcsecs)
    T3_diffV2_23, T3_diffV3_23, T3bench_V2_listP2, T3bench_V3_listP2 = compare2ref(case, bench_starP1, bench_V2P2, bench_V3P2, stars, T3_V2_23, T3_V3_23, arcsecs=diffs_in_arcsecs)
    T3_diffV2_15, T3_diffV3_15, _, _ = compare2ref(case, bench_starP1, bench_V2P1, bench_V3P1, stars, T3_V2_15, T3_V3_15, arcsecs=diffs_in_arcsecs)
    T3_diffV2_25, T3_diffV3_25, _, _ = compare2ref(case, bench_starP1, bench_V2P2, bench_V3P2, stars, T3_V2_25, T3_V3_25, arcsecs=diffs_in_arcsecs)
    T3_diffV2_17, T3_diffV3_17, _, _ = compare2ref(case, bench_starP1, bench_V2P1, bench_V3P1, stars, T3_V2_17, T3_V3_17, arcsecs=diffs_in_arcsecs)
    T3_diffV2_27, T3_diffV3_27, _, _ = compare2ref(case, bench_starP1, bench_V2P2, bench_V3P2, stars, T3_V2_27, T3_V3_27, arcsecs=diffs_in_arcsecs)
    if debug:
        print ("TEST 3: ")
        print ("transformations: detector P1 and P2 (x, y),  sky P1 and P2 (V2, V3),  true P1 and P2 (V2, V3)")
        print (" Centroid window 3: ", x13[0], y13[0], x23[0], y23[0], T3_V2_13[0], T3_V3_13[0], T3_V2_23[0], T3_V3_23[0], bench_V2P1[0], bench_V3P1[0], bench_V2P2[0], bench_V3P2[0])
        print (" Centroid window 5: ", x15[0], y15[0], x25[0], y25[0], T3_V2_13[0], T3_V3_13[0], T3_V2_23[0], T3_V3_23[0], bench_V2P1[0], bench_V3P1[0], bench_V2P2[0], bench_V3P2[0])
        print (" Centroid window 7: ", x17[0], y17[0], x27[0], y27[0], T3_V2_13[0], T3_V3_13[0], T3_V2_23[0], T3_V3_23[0], bench_V2P1[0], bench_V3P1[0], bench_V2P2[0], bench_V3P2[0])
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


def runTest_and_append_results(test2run, data4test, Vs, diffs, benchVs, filter_input, tilt, diffs_in_arcsecs, debug):
    """
    This function runs the test for the specified detector and sliced arrays, and appends it to the results.
    Args:
        test2run: string, either 'T1', 'T2', or 'T3'
        data4test: list of 7 elements (detector, transfromation direction, case, stars sample,
                   positions 1 and 2data, list of benchmark values for position 1 or 2, and list of
                   benchmark values for V2 and V3)
        Vs: list transformed V2 and V3s for centroid window sizes 3, 5, and 7
        diffs: list of 6 lists of differences with respect to benchmark V2s and V3s
        benchVs: list of 2 lists, benchmark V2s and V3s
        filter_input: string
        tilt: True or false
        diffs_in_arcsecs: True or False
        debug: True or False

    Returns:
        4 lists: P1P2data, Vs, diffs, benchVs
    """
    detector, transf_direction, case, stars, P1P2data, bench_starP1, benchV23 = data4test
    bench_V2P1, bench_V3P1, bench_V2P2, bench_V3P2 = benchV23
    avg_benchV2 = (bench_V2P1 + bench_V2P2)/2.0
    avg_benchV3 = (bench_V3P1 + bench_V3P2)/2.0
    avg_benchV23 = [avg_benchV2, avg_benchV3]
    T_V2_3, T_V3_3, T_V2_5, T_V3_5, T_V2_7, T_V3_7 = Vs
    T_diffV2_3, T_diffV3_3, T_diffV2_5, T_diffV3_5, T_diffV2_7, T_diffV3_7 = diffs
    Tbench_V2_list, Tbench_V3_list = benchVs
    if test2run == "T1":
        transformations, diffs, benchVs_list = TEST1(detector, transf_direction, stars, case, bench_starP1, avg_benchV23, P1P2data,
                                                    filter_input, tilt, diffs_in_arcsecs, debug)
    if test2run == "T2":
        transformations, diffs, benchVs_list = TEST2(detector, transf_direction, stars, case, bench_starP1, avg_benchV23, P1P2data,
                                                    filter_input, tilt, diffs_in_arcsecs, debug)
    if test2run == "T3":
        transformations, diffs, benchVs_list = TEST3(detector, transf_direction, stars, case, bench_starP1, benchV23, P1P2data,
                                                    filter_input, tilt, diffs_in_arcsecs, debug)
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


def runTEST(test2run, detectors, transf_direction, case, stars, P1P2data, bench_starP1, trueVsP1, trueVsP2,
          filter_input, tilt, diffs_in_arcsecs, debug):
    """
    This function runs the test for both detectors and returns the results for the star sample.
    Pier_corr is set to False in case it was corrected before.
    Args:
        test2run: string, either 'T1', 'T2', or 'T3'
        detectors: integer, 491 or 492
        transf_direction: string, 'forward' or 'backward'
        case: string, for example 'Scene2_rapid_real_bgFrac'
        stars: list, the star sample
        P1P2data: list of positions 1 and 2 for x and y for centroid window size 3, 5, and 7
        bench_starP1: list of benchmark positions for either position 1 or 2
        trueVsP1: list of true V2s and V3s of position 1
        trueVsP2: list of true V2s and V3s of position 2
        filter_input: string
        tilt: True or false
        diffs_in_arcsecs: True or False
        debug: True or False

    Returns:
        resultsTEST = [P1P2data, Vs, diffs, benchVs]
        List of 4 lists, list of positions 1 and 2 for x and y for centroid window size 3, 5, and 7, V2s and V3,
        differences of measured with respect to benchmark V2s and V3s, and list of benchmark V2s and V3s.
    """
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
    P1P2data, Vs, diffs, benchVs = runTest_and_append_results(test2run, data4test, Vs, diffs, benchVs,
                                                                filter_input, tilt, diffs_in_arcsecs, debug)
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
        P1P2data, Vs, diffs, benchVs = runTest_and_append_results(test2run, data4test, Vs, diffs, benchVs,
                                                                  filter_input, tilt, diffs_in_arcsecs, debug)
    resultsTEST = [P1P2data, Vs, diffs, benchVs]
    return resultsTEST


def combine2arrays(arr1, arr2, combined_arr):
    """ This function combines 2 (of positions 1 and 2) into a single array. This is to have both positions
    for the calculation of least squares and sigma rejection in Test 3.
    Args:
        arr1: Either x or y measurement of position 1
        arr2: Either x or y measurement of position 2
        combined_arr: The already defined numpy array that the measurements are to be appended to.

    Returns:
         combined_arr = A combined numpy array of one dimension.
    """
    for item in arr1:
        combined_arr = np.append(combined_arr, item)
    for item in arr2:
        combined_arr = np.append(combined_arr, item)
    return combined_arr


def get_rejected_stars(stars_sample, rejected_elements_idx):
    """
    This function finds which are the star numbers of the rejected stars.
    Args:
        stars_sample: list of stars in the sample
        rejected_elements_idx: the index list of the rejected stars
    Returns:
        rejected_elements = The list of rejected stars
    """
    rejected_elements = []
    if len(rejected_elements_idx) != len(stars_sample):
        for i in rejected_elements_idx:
            rejected_elements.append(stars_sample[i])
    return rejected_elements


def get_stats(case, T_transformations, T_diffs, T_benchVs_list, Nsigma, max_iterations):
    """
    This function obtains the standard deviations through regular statistics as well as through
    a sigma clipping algorithm and an iterative least square algorithm. It also obtains the minimum
    differences from centroid window sizes 3, 5, and 7, and returns the counter for each.
    Args:
        case: string, for example 'Scene2_rapid_real_bgFrac'
        T_transformations: list of sky transformation for centroid window sizes 3, 5, and 7
        T_diffs: list of differences of measured with respect to benchmark sky values
        T_benchVs_list: list of benchmark V2 and V3s
        Nsigma: float or integer, the number of sigmas to reject
        max_iterations: integer, maximum iterations for the Nsigma routine

    Returns:
        results_stats: List with standard deviations and means, dictionary, minimum differences from centroid window
                       sizes 3, 5, and 7, and their repetitions, benchmark values, results from least squares
                       and Nsigma rejection routines, and list of rejected elements from each of those 2 routines.
        results_stats = [st_devsAndMeans, diff_counter, bench_values, sigmas_deltas, sigma_reject,
                        rejected_elementsLS, rejected_elementsNsig]
    """
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
    Tstdev_V2_3, Tmean_V2_3 = find_std(T_diffV2_3)
    Tstdev_V2_5, Tmean_V2_5 = find_std(T_diffV2_5)
    Tstdev_V2_7, Tmean_V2_7 = find_std(T_diffV2_7)
    Tstdev_V3_3, Tmean_V3_3 = find_std(T_diffV3_3)
    Tstdev_V3_5, Tmean_V3_5 = find_std(T_diffV3_5)
    Tstdev_V3_7, Tmean_V3_7 = find_std(T_diffV3_7)
    Tbench_V2, Tbench_V3 = np.array(Tbench_V2_list), np.array(Tbench_V3_list)
    # get the minimum of the differences
    T_min_diff, T_counter = get_mindiff(T_diffV2_3, T_diffV2_5, T_diffV2_7)
    # calculate least squares but first convert to MSA center
    T_V2_3, T_V3_3, Tbench_V2, Tbench_V3 = convert2MSAcenter(T_V2_3, T_V3_3, Tbench_V2, Tbench_V3)
    T_V2_5, T_V3_5, _, _ = convert2MSAcenter(T_V2_5, T_V3_5, Tbench_V2, Tbench_V3)
    T_V2_7, T_V3_7, _, _ = convert2MSAcenter(T_V2_7, T_V3_7, Tbench_V2, Tbench_V3)
    # to express in arcsecs multiply by 3600.0
    TLSdeltas_3, TLSsigmas_3, TLSlines2print_3, rejected_elements_3 = lsi.ls_fit_iter(max_iterations, T_V2_3*3600.0, T_V3_3*3600.0, Tbench_V2*3600.0, Tbench_V3*3600.0)
    TLSdeltas_5, TLSsigmas_5, TLSlines2print_5, rejected_elements_5 = lsi.ls_fit_iter(max_iterations, T_V2_5*3600.0, T_V3_5*3600.0, Tbench_V2*3600.0, Tbench_V3*3600.0)
    TLSdeltas_7, TLSsigmas_7, TLSlines2print_7, rejected_elements_7 = lsi.ls_fit_iter(max_iterations, T_V2_7*3600.0, T_V3_7*3600.0, Tbench_V2*3600.0, Tbench_V3*3600.0)
    # Do N-sigma rejection
    TsigmaV2_3, TmeanV2_3, TsigmaV3_3, TmeanV3_3, TnewV2_3, TnewV3_3, Tniter_3, Tlines2print_3, rej_elements_3 = Nsigma_rejection(Nsigma, T_diffV2_3, T_diffV3_3, max_iterations)
    TsigmaV2_5, TmeanV2_5, TsigmaV3_5, TmeanV3_5, TnewV2_5, TnewV3_5, Tniter_5, Tlines2print_5, rej_elements_5 = Nsigma_rejection(Nsigma, T_diffV2_5, T_diffV3_5, max_iterations)
    TsigmaV2_7, TmeanV2_7, TsigmaV3_7, TmeanV3_7, TnewV2_7, TnewV3_7, Tniter_7, Tlines2print_7, rej_elements_7 = Nsigma_rejection(Nsigma, T_diffV2_7, T_diffV3_7, max_iterations)
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
    results_stats = [st_devsAndMeans, diff_counter, bench_values, sigmas_deltas, sigma_reject,
                     rejected_elementsLS, rejected_elementsNsig]
    return results_stats


def printTESTresults(stars_sample, case, test2perform, diffs_in_arcsecs, Tstdev_Vs, Tmean_Vs,
                     T_diff_counter, save_text_file, TLSlines2print, Tlines2print, Tbench_Vs_list,
                     T_Vs, T_diffVs, rejected_elementsLS, rejected_eleNsig,
                     background_method, background2use, path4results):
    """
    This function writes to text the results from the Test performed.
    Args:
        stars_sample: list of the stars to be studied
        case: string, for example 'Scene2_rapid_real_bgFrac'
        test2perform: string, either 'T1', 'T2', or 'T3'
        diffs_in_arcsecs: True or False (if False results are degrees)
        Tstdev_Vs: list of test standard deviations for V2 and V3 for centroid window sizes 3, 5, and 7
        Tmean_Vs: list of mean V2 and V3 for centroid window sizes 3, 5, and 7
        T_diff_counter: dictionary, minimum differences from centroid window sizes 3, 5, and 7, and their repetitions
        save_text_file: True or False
        TLSlines2print: list of lines to print from least square routine results
        Tlines2print: list of lines to print on screen and as column headers in the text file
        Tbench_Vs_list: list of benchmark V2s and V3s for the star sample
        T_Vs: list of measured V2s and V3s
        T_diffVs: list of differences of measured V2, V3s with respect to benchmark
        rejected_elementsLS: list of rejected star numbers from the least squares routine
        rejected_eleNsig: list of rejected star numbers from the Nsigma routine
        background_method: None or string ('fix' or 'frac')
        background2use: value to use for fixed or fractional background methods (0.0 is used if background_method=None)
        path4results: string with the path for the resulting text file

    Returns:
        Nothing. The text file is written.
    """
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
    line3bisAb = "   centroid window 3: {} ".format(rejected_elements_3)
    line3bisAc = "   centroid window 5: {} ".format(rejected_elements_5)
    line3bisAd = "   centroid window 7: {} ".format(rejected_elements_7)
    line3bisAe = "   centroid window 3: {} ".format(Nsig_rej_elements_3)
    line3bisAf = "   centroid window 5: {} ".format(Nsig_rej_elements_5)
    line3bisAg = "   centroid window 7: {} ".format(Nsig_rej_elements_7)
    # Print number of repetitions to find best centroid window
    line3bisB = "\n *** Repetitions Diffs: {}".format(T_counter)
    line4 = "{:<5} {:<20} {:<40} {:<40} {:<38} {:<28} {:<7}".format(
                    "Star", "BG_value", "Pos_centroid_win_3", "Pos_centroid_win_5", "Poscentroid_win_7",
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
        to.write('       Centroid window 3:  \n')
        for line2print in TLSlines2print_3:
            to.write(line2print+"\n")
        to.write('       Centroid window 5:  \n')
        for line2print in TLSlines2print_5:
            to.write(line2print+"\n")
        to.write('       Centroid window 7:  \n')
        for line2print in TLSlines2print_7:
            to.write(line2print+'\n')
        # print standard deviations and means after n-sigma rejection
        to.write('\n * From N-sigma rejection routine:  \n')
        to.write(line3bisAa+"\n")
        to.write(line3bisAe+"\n")
        to.write(line3bisAf+"\n")
        to.write(line3bisAg+"\n")
        to.write(" Centroid window 3:  \n")
        for line2print in Tlines2print_3:
            to.write(line2print+"\n")
        to.write(" Centroid window 5:  \n")
        for line2print in Tlines2print_5:
            to.write(line2print+"\n")
        to.write(" Centroid window 7:  \n")
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
        j += 1
        if st == stars_sample[-1]:
            j = 0
    if save_text_file:
        to.close()
        print (" * Results saved in file: ", txt_out)


def writePixPos(save_text_file, show_centroids, output_file, lines4screenandfile,
                stars_sample, background2use, data2write):
    """
    This function writes a text file with the measured centroids (in pixel space).
    Args:
        save_text_file: True or False
        show_centroids: True or False
        output_file: name of the output file
        lines4screenandfile: list of lines to be shown on screen and in the file at the top (column headers)
        stars_sample: list of stars to be studied
        background2use: a float value with the background value to be used as fixed or fractional
        data2write: list of 7 lists of the same length (x and y positions, true center, lower left coordinates,
                    star magnitudes, centroid window size with minimum difference with respect to true
                    center -- either 3, 5, or 7)

    Returns:
        Nothing. The text file is written.
    """
    line0, line2a, line2b = lines4screenandfile
    x_pixpos, y_pixpos, corr_true_center_centroid, loleftcoords, mag, min_diff_pixposX, min_diff_pixposY = data2write
    x3, x5, x7 = x_pixpos
    y3, y5, y7 = y_pixpos
    counterX = collections.Counter(min_diff_pixposX)
    counterY = collections.Counter(min_diff_pixposY)
    line1a = "Counter for X positions: {}".format(counterX)
    line1b = "Counter for Y positions: {}".format(counterY)
    if show_centroids:
        print()
        print(line0)
        print(line1a)
        print(line1b)
        print(line2a)
        print(line2b)
    if save_text_file:
        f = open(output_file, "a")
        f.write(line1a+"\n")
        f.write(line1b+"\n")
        f.write(line2a+"\n")
        f.write(line2b+"\n")
        f.close()
    for i, st in enumerate(stars_sample):
        line3 = "{:<5} {:<10} {:<14} {:<16} {:<14} {:<16} {:<14} {:<18} {:<16} {:<16} {:<10} {:<10} {:<11.2f} {:<2} {:<10}".format(
                                                    st, background2use,
                                                    x3[i], y3[i], x5[i], y5[i], x7[i], y7[i],
                                                    corr_true_center_centroid[i][0], corr_true_center_centroid[i][1],
                                                    loleftcoords[i][0], loleftcoords[i][1],
                                                    mag[i],
                                                    min_diff_pixposX[i], min_diff_pixposY[i])
        if save_text_file:
            f = open(output_file, "a")
            f.write(line3+"\n")
            f.close()
        if show_centroids:
            print(line3)


def remove_bad_stars(stars_sample):
    """ This function reads the text files of bad stars, compares the sample data, removes
    the bad stars, and returns the sample without bad stars.
    Args:
        stars_sample = list of stars to be studied
    Returns:
        stars_sample = list of stars with the 'bad stars' removed.
    """
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


