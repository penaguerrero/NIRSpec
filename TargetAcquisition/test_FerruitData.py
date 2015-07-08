# Import necessary modules
from __future__ import division     
from __future__ import print_function
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt

# Tommy's code
import tautils as tu
import jwst_targloc as jtl


"""

This script tests selected synthetic data from Pierre Ferruit in order to compare it to the centroids
calulated by him and by IDL code written by Tony Keyes.

"""

# Benchmark values:
star101_pierre = [16.865, 16.857]
star111_pierre = [16.262, 16.697] 


""" Running Tommy's code on these 2 Stars """
# Read input image: 3-ram integration master image
star = 'PierreFerruitData/TA_cutouts/postageout_star_     101 quad_       3 quad_star        1.fits'
#img = fits.open(star)
#img.info()
#hdr = img[0].header
#for item in hdr:
#    print (item)
#master_img = img[0].data
master_img = fits.getdata(star, 0)
print ('Master image shape: ', np.shape(master_img))
#print (master_img)
#tu.display_ns_psf(master_img[2,:,:])   # to show the individual images

# Obtain and display the combined FITS image that combines all frames into one image.
psf = tu.readimage(master_img)
#tu.display_ns_psf(psf)#, vlim=(0.001, 0.01))

# Test checkbox piece
cb_cen, cb_hw = jtl.checkbox_2D(psf, 3, debug=True)
print('Got coarse location. \n')
#raw_input()

# Checkbox center, in base 1
print('Checkbox Output:')
print('Checkbox center: [{}, {}]'.format(cb_cen[0], cb_cen[1]))
print('Checkbox halfwidths: xhw: {}, yhw: {}'.format(cb_hw[0], cb_hw[1]))
print()

# Now calculate the centroid based on the checkbox region calculated above
cb_centroid, cb_sum = jtl.centroid_2D(psf, cb_cen, cb_hw, max_iter=50, threshold=1e-5, debug=True)
print('Final sum: ', cb_sum)
print('cb_centroid: ', cb_centroid)
print()

'''
# Plot the rows/columns encompasing the centroid to verify psf shape
x_radial_1 = psf[:, 15]
x_radial_2 = psf[:, 16]
y_radial_1 = psf[15, :]
y_radial_2 = psf[16, :]
fig, ax = plt.subplots(figsize=(8, 8))
ax.set_xlim(0, 31)
ax.set_xlabel('Pixel Location')
ax.set_ylabel('Counts (of total 1)')
ax.plot(xrange(32), x_radial_1, ls='-', lw=2., label='Column 15')
ax.plot(xrange(32), x_radial_2, ls='-', lw=2., label='Column 16')
ax.plot(xrange(32), y_radial_1, ls='-', lw=2., label='Row 15')
ax.plot(xrange(32), y_radial_2, ls='-', lw=2., color='purple', label='Row 16')
ax.legend(numpoints=1, loc=1, prop={"size":"medium"})
plt.show()
'''