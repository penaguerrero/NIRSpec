from __future__ import print_function, division
import numpy as np

# Header
__author__ = "Maria A. Pena-Guerrero"
__version__ = "1.0"

# HISTORY
#    1. Jul 2015 - Vr. 1.0: Initial Python translation of IDL code.

"""
This script is a Python translation from the ta_lsfitdbl.pro by Tony Keyes, which 
in turn is a sample code from the fitting algorithm with roll correction 
(Jakobsen, 2006) draft coding for spatial offset sigmas added.
     -- needs verification and does not contain calculation for roll sigma

Keyword arguments:
    niter       -- number of maximum iterations
    x_input     -- input measured x-centroid converted to V2 
    y_input     -- input measured y-centroid converted to V3 
        xcentroids  -- input measured x-centroid converted to V2 
        ycentroids  -- input measured y-centroid converted to V3 
    xtrue       -- numpy array of true V2 positions
    ytrue       -- numpy array of true V3 positions

Output(s):
    deltas = [delta_x, delta_y, delta_theta]
    sigmas = [sigma_x, sigma_y, sigma_theta]

Example usage:
    import least_squares_iterate as lsi
    deltas, sigmas = lsi.ls_fit_iter(niter, x_input, y_input, xtrue, ytrue)


*** Testing suite of the script at the bottom

"""

def ls_fit_iter(niter, xt, yt, x, y):
    """
    Inputs:  niter  = number of max iterations
             xt     = input measured x-centroid converted to V2
             yt     = input measured y-centroid converted to V3
             x      = numpy array of true V2 positions
             y      = numpy array of true V3 positions
    """  
    # do up to niter iterations of sigma-clipping (first time through is 
    # initial calculation, then up to niter iterations)
    n = len(x)
    for nit in range(niter):
        # Initialize the sums
        sum_tot = 0.0
        sum_x = 0.0
        sum_y = 0.0
        sum_xt = 0.0
        sum_yt = 0.0
        sum_xt2 = 0.0
        sum_yt2 = 0.0
        sum_xyt = 0.0
        sum_xty = 0.0
        
        for i in range(n):
            sum_tot = sum_tot + 1.0
            sum_x = sum_x + x[i]
            sum_y = sum_y + y[i]
            sum_xt = sum_xt + xt[i]
            sum_yt = sum_yt + yt[i]
            sum_xt2 = sum_xt2 + xt[i]*xt[i]
            sum_yt2 = sum_yt2 + yt[i]*yt[i]
            sum_xyt = sum_xyt + x[i]*yt[i]
            sum_xty = sum_xty + xt[i]*y[i]
        
        det = sum_tot*sum_tot*(sum_xt2 + sum_yt2) - sum_tot*(sum_xt*sum_xt + sum_yt*sum_yt)
        
        delta_x = (sum_tot*(sum_xt2 + sum_yt2) - sum_xt*sum_xt) * (sum_x - sum_xt)
        delta_x = delta_x - (sum_y-sum_yt)*sum_xt*sum_yt - sum_tot*sum_yt*(sum_xyt-sum_xty)
        delta_x = delta_x/det    # offset in x-coordinate
        
        delta_y = (sum_tot*(sum_xt2+sum_yt2) - sum_yt*sum_yt) * (sum_y-sum_yt)
        delta_y = delta_y - (sum_x-sum_xt)*sum_xt*sum_yt - sum_tot*sum_xt*(sum_xty-sum_xyt)
        delta_y = delta_y/det    # offset in y-coordinate
        
        delta_theta = sum_tot*((sum_xt*sum_y-sum_xt*sum_yt) + sum_tot*(sum_xyt-sum_xty))
        delta_theta = delta_theta/det    # roll angle correction  (what units??)
        
        # outputs:  delta_x, delta_y, delta_theta, sigma_x, sigma_y, sigma_theta
        print ('(least_squares_iterate):  iteration no.: ', nit)
        print ('(least_squares_iterate):  delta_x = {}   delta_y = {}   delta_theta = {}'.format(delta_x, delta_y, delta_theta))
        deltas = [delta_x, delta_y, delta_theta]
        
        # verify this coding for sigma_xtrue and sigma_ytrue
        sum_delta_x2 = 0.0
        sum_delta_y2 = 0.0
        sum_delta_theta2 = 0.0
        for i in range(n):
            sum_delta_x2 = sum_delta_x2 + (-xt[i] + x[i] - delta_x) * (-xt[i] + x[i] - delta_x) 
            sum_delta_y2 = sum_delta_y2 + (-yt[i] + y[i] - delta_y) * (-yt[i] + y[i] - delta_y) 
        
        sigma_x = np.sqrt(sum_delta_x2/n)   # sigma for xtrue-offset  -- is this right?
        sigma_y = np.sqrt(sum_delta_y2/n)   # sigma for ytrue-offset  -- is this right?
        
        # for now set sigma_thera to bogus value  (we don't presently know how to calculate it)
        sigma_theta = -999.0
        
        print ('(least_squares_iterate):  sigma_x = {}   sigma_y = {}   sigma_theta = {}'.format(sigma_x, sigma_y, sigma_theta))
        sigmas = [sigma_x, sigma_y, sigma_theta]
        
        # calculate new best position and residuals for each coordinate (neglect any roll correction for now)
        xnewpos = x - delta_x
        xdiff = xnewpos - xt
        ynewpos = y - delta_y
        ydiff = ynewpos - yt
        
        # Rejection sequence:
        # reject any residuals that are not within 3*sigma in EITHER coordinate
        # (this is slightly more restrictive than doing this for the vector sum, 
        # but easier to implement right now)
        thres_x = 3.0*sigma_x
        thres_y = 3.0*sigma_y
        var_clip = xdiff[(np.where((np.abs(xdiff)<=thres_x) & (np.abs(ydiff)<=thres_y)))]
        xcentroids_new = xt[(np.where((np.abs(xdiff)<=thres_x) & (np.abs(ydiff)<=thres_y)))]
        ycentroids_new = yt[(np.where((np.abs(xdiff)<=thres_x) & (np.abs(ydiff)<=thres_y)))]
        x_new = x[(np.where((np.abs(xdiff)<=thres_x) & (np.abs(ydiff)<=thres_y)))]
        y_new = y[(np.where((np.abs(xdiff)<=thres_x) & (np.abs(ydiff)<=thres_y)))]
        
        # This commented following section would do rejection for 3.0*sigma AND a specified
        # physical distance threshold zerop7, a 0.7 arcsec threshold
        """
        zerop7 = 0.7
        var_clip = xdiff[(np.where((np.abs(xdiff)<=thres_x) & (np.abs(ydiff)<=thres_y) & (np.abs(xdiff)<=zerop7) & (np.abs(ydiff)<=zerop7)))]
        xcentroids_new = xt[(np.where((np.abs(xdiff)<=thres_x) & (np.abs(ydiff)<=thres_y) & (np.abs(xdiff)<=zerop7) & (np.abs(ydiff)<=zerop7)))]
        ycentroids_new = yt[(np.where((np.abs(xdiff)<=thres_x) & (np.abs(ydiff)<=thres_y) & (np.abs(xdiff)<=zerop7) & (np.abs(ydiff)<=zerop7)))]
        x_new = x[(np.where((np.abs(xdiff)<=thres_x) & (np.abs(ydiff)<=thres_y) & (np.abs(xdiff)<=zerop7) & (np.abs(ydiff)<=zerop7)))]
        y_new = y[(np.where((np.abs(xdiff)<=thres_x) & (np.abs(ydiff)<=thres_y) & (np.abs(xdiff)<=zerop7) & (np.abs(ydiff)<=zerop7)))]
        """
        if len(xcentroids_new) == len(xt):
            break   # exit the loop since no additional rejections on this iteration
        else:
            xt = xcentroids_new
            yt = ycentroids_new
            x = x_new
            y = y_new
    
    return deltas, sigmas
    
    # Still do not know how to do delta_theta sigma  -- is this calculation needed?


# Print diagnostic load message
print("(least_squares_iterate): Least squares iteration algorithm Version {} loaded!".format(__version__))


testing = False
if testing: 
    # Set test values  for arrays
    n = 10            # Number of positions
    xtrue = np.array(range(10))     # true x-coordinate of each reference star: from 0 to 9
    ytrue = np.array(range(1, 11))  # true y-coordinate of each reference star: from 1 to 10
    xinput = xtrue + 0.02     # measured centroid x-coordinate of each reference star
    yinput = ytrue + 0.01     # measured centroid y-coordinate of each reference star
    deltas, sigmas = ls_fit_iter(n, xinput, yinput, xtrue, ytrue)
    """
    With these parameters output should be:
        (least_squares_iterate): Least squares iteration algorithm Version 1.0 loaded!
        (least_squares_iterate):  iteration no.:  0
        (least_squares_iterate):  delta_x = -0.02   delta_y = -0.01   delta_theta = 7.57912251477e-16
        (least_squares_iterate):  sigma_x = 4.564982887e-15   sigma_y = 3.69348008392e-15   sigma_theta = -999.0
    """
    

