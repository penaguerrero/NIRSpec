from __future__ import print_function, division
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt

# Tommy's code
import tautils as tu
import jwst_targloc as jtl

# Header
__author__ = "Maria A. Pena-Guerrero"
__version__ = "1.0"

def run_recursive_centroids(psf, background, xwidth_list, ywidth_list, checkbox_size, max_iter, 
                            threshold, determine_moments, debug, display_master_img, vlim=()):   
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


def transform2fulldetector(detector, centroid_in_full_detector, cb_centroid_list, ESA_center, true_center, perform_avgcorr=False):
    """
    Transform centroid coordinates into full detector coordinates.
    
    Keyword arguments:
    detector                   -- Either 491 or 492
    centroid_in_full_detector  -- Resuling coordinates in terms of full detector, True or False
    cb_centroid_list           -- Centroid determined by target acquisition (TA) algorithm in terms of 32 by 32 pixels
    ESA_center                 -- Centroid determined with the ESA version of TA algorithm in terms of full detector
    true_center                -- Actual (true) position of star in terms of full detector  
    perform_avgcorr            -- Add the average correction given by Pierre
    
    Output(s):
    cb_centroid_list_fulldetector  -- List of centroid locations determined with the TA algorithm in 
                                      terms of full detector. List is for positions determined with
                                      3, 5, and 7 checkbox sizes. 
    
    """
    # Corrections for offsets in positions (see section 2.5 of Technical Notes in Documentation directory)
    corrected_x = true_center[0]
    corrected_y = true_center[1]
    if perform_avgcorr:
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
    loleftcoords = [loleft_x, loleft_y]

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
    return corr_true_center_centroid, corr_cb_centroid_list, loleftcoords, differences_true_TA


def write2file(data2write, lines4screenandfile):
    line0, line0a, line0b = lines4screenandfile
    save_text_file, output_file, st, bg, corr_cb_centroid_list, corr_true_center_centroid, loleftcoords, factor, differences_true_TA = data2write
    line1 = "{:<5} {:<10} {:<14} {:<16} {:<14} {:<16} {:<14} {:<16} {:<12} {:<14} {:<10} {:<14} {:<10.2f} {:<20} {:<22} {:<20} {:<22} {:<20} {:<22}\n".format(
                                                    st, bg, 
                                                    corr_cb_centroid_list[0][0], corr_cb_centroid_list[0][1],
                                                    corr_cb_centroid_list[1][0], corr_cb_centroid_list[1][1],
                                                    corr_cb_centroid_list[2][0], corr_cb_centroid_list[2][1],
                                                    corr_true_center_centroid[0], corr_true_center_centroid[1],
                                                    loleftcoords[0], loleftcoords[1],
                                                    factor,
                                                    differences_true_TA[0][0][0], differences_true_TA[0][0][1],
                                                    differences_true_TA[0][1][0], differences_true_TA[0][1][1],
                                                    differences_true_TA[0][2][0], differences_true_TA[0][2][1])
    if save_text_file:
        f = open(output_file, "a")
        f.write(line1)
        f.close()
    print(line0)
    print(line0a)
    print(line0b)
    print(line1) 


def read_listfile(list_file_name, detector, background_method):    
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
    #mag = 2.5*np.log10(factor) + 18.0
    # Get the correct slices according to detector
    if detector == 491:   # slice from star 101 to 200
        star_number, xpos, ypos, factor, mag = star_number[100:], xpos[100:], ypos[100:], factor[100:], mag[100:]
    elif detector == 492:   # slice from star 1 to 100
        star_number, xpos, ypos, factor, mag = star_number[:100], xpos[:100], ypos[:100:], factor[:100], mag[:100]
    bg_method = background_method
    if background_method is None:   # convert the None value to string
        bg_method = 'None'
    return star_number, xpos, ypos, factor, mag, bg_method
    

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


def find_std(arr):
    """ This function determines the standard deviation of the given array. """
    N = float(len(arr))
    mean = np.sum(arr) / N
    diff2meansq_list = []
    for a in arr:
        diff = a - mean
        diffsq = diff * diff
        diff2meansq_list.append(diffsq)
    std = ( 1.0/(N-1.0) * sum(diff2meansq_list) )**(0.5)
    #print ('sigma = ', std, '    mean = ', mean)
    return std, mean


def get_frac_stdevs(frac_data):
    sig3, mean3 = [], []
    sig5, mean5 = [], []
    sig7, mean7 = [], []
    for f in frac_data:
        s3, m3 = find_std(f[1])
        s5, m5 = find_std(f[3])
        s7, m7 = find_std(f[5])
        sig3.append(s3)
        sig5.append(s5)
        sig7.append(s7)
        mean3.append(m3)
        mean5.append(m5)
        mean7.append(m7)
    return sig3, mean3, sig5, mean5, sig7, mean7


def display_centroids(st, case, psf, corr_true_center_centroid, corr_cb_centroid_list, show_disp, vlims=None, savefile=None):    
    fig_title = "star_"+str(st)+"_"+case
    # Display both centroids for comparison.
    _, ax = plt.subplots(figsize=(8, 8))
    ax.set_title(fig_title)
    ax.autoscale(enable=False, axis='both')
    ax.imshow(psf, cmap='gray', interpolation='nearest')
    ax.set_ylim(0.0, np.shape(psf)[0])
    ax.set_xlim(0.0, np.shape(psf)[1])
    ax.plot(corr_cb_centroid_list[0][0], corr_cb_centroid_list[0][1], marker='*', ms=20, mec='black', mfc='blue', ls='', label='Checkbox=3')
    ax.plot(corr_cb_centroid_list[1][0], corr_cb_centroid_list[1][1], marker='*', ms=17, mec='black', mfc='green', ls='', label='Checkbox=5')
    ax.plot(corr_cb_centroid_list[2][0], corr_cb_centroid_list[2][1], marker='*', ms=15, mec='black', mfc='red', ls='', label='Checkbox=7')
    ax.plot(corr_true_center_centroid[0], corr_true_center_centroid[1], marker='o', ms=8, mec='black', mfc='yellow', ls='', label='True Centroid')
    # Shrink current axis by 10%
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.9, box.height])
    ax.legend(loc='upper right', bbox_to_anchor=(1.26, 1.0), prop={"size":"small"})   # put legend out of the plot box   
    # Add plot with different sensitivity limits
    if vlims is None:
        vlims = (1, 10)
    fig, ax = plt.subplots(figsize=(8, 8))
    ax.set_title(fig_title)
    ax.autoscale(enable=False, axis='both')
    ax.imshow(psf, cmap='gray', interpolation='nearest')
    ax.set_ylim(0.0, np.shape(psf)[0])
    ax.set_xlim(0.0, np.shape(psf)[1])
    ax.plot(corr_cb_centroid_list[0][0], corr_cb_centroid_list[0][1], marker='*', ms=20, mec='black', mfc='blue', ls='', label='Checkbox=3')
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
    if savefile is not None:
        path4fig = "PFforMaria/centroid_figs/"
        if "scene1" in fig_title:
            if "slow" in fig_title:
                if "real" in fig_title:
                    in_dir = "Scene1_slow_real/"
                else:
                    in_dir = "Scene1_slow_nonoise/"
            elif "rapid" in fig_title:
                if "real" in fig_title:
                    in_dir = "Scene1_rapid_real/"
                else:
                    in_dir = "Scene1_rapid_nonoise/"
        if "scene2" in fig_title:
            if "slow" in fig_title:
                if "real" in fig_title:
                    in_dir = "Scene2_slow_real/"
                else:
                    in_dir = "Scene2_slow_nonoise/"
            elif "rapid" in fig_title:
                if "real" in fig_title:
                    in_dir = "Scene2_rapid_real/"
                else:
                    in_dir = "Scene2_rapid_nonoise/"
        fig_name = path4fig+in_dir+fig_title+".jpg"
        fig.savefig(fig_name)
        print ("Figure ", fig_name, " was saved!")
    
    
def get_cutouts(ref_star_position, detector_x, detector_y):
    """
    This function does the cutouts for the reference star position.
    Inputs:
        - ref_star_position = [true_x_center, true_y_center]
        - detector_x = array of x detector values
        - detector_y = array of y detector values
    """
    xc, yc = ref_star_position[0], ref_star_position[1]
    # Get the box coordinates in terms of full detector (subtract 16.0 because indexing
    # from centroid function starts with 1) 
    lo_x = np.floor(xc) - 16.0
    lo_y = np.floor(yc) - 16.0
    up_x = lo_x + 32.0
    up_y = lo_y + 32.0
    cutout_x = detector_x[(detector_x >= lo_x) & (detector_x <= up_x)]
    cutout_y = detector_y[(detector_y >= lo_y) & (detector_y <= up_y)]
    print (len(cutout_x), len(cutout_y))
    
    

# Print diagnostic load message
print("(testing_functions): testing functions script Version {} loaded!".format(__version__))

