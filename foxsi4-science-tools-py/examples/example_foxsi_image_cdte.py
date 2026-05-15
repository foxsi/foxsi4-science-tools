"""
This code provides a simple example of making an image from FOXSI-4 
Level 3 CdTe (HXR) data.

To make an image, create a 2D histogram in X and Y spatial dimensions.
Don't forget to restrict the event list to the time and energy ranges
you would like.

If you do not restrict the time range, you might end up including a 
little bit of pointing drift at the very start and end of the observation,
or during slews.
"""

import numpy as np
from copy import deepcopy
from astropy.table import Table
from astropy.time import Time
from astropy.io import fits
import matplotlib.pyplot as plt

# Choose file and directory
# File should be a Level 3 FOXSI-4 CdTe data file. (Example: 'cdte1_level3_version1.fits')
# Example here chooses a CdTe1 file. Also define a string with this information.
# Don't forget to change the directory path to match YOUR path to the data.
dir = "/Users/glesener/Dropbox/data/foxsi/20240417/cdte_lev3/"
file = 'cdte1_level3_v2.fits'
name = 'CdTe1'

# Read in event data from the file
with fits.open( dir + file ) as hdul:
    hdus = deepcopy(hdul)

evt = Table(hdus[1].data)
# EVT now holds the event list for the selected detector.
print(len(evt),'events read from file for '+name)

# Restrict to events with nominal quality flag
evt = evt[ evt["flag_quality"]==1 ]
print(len(evt),'perfect events')

# Select a time range for the image
# Note that we might change the time format of the event list to a different format.
# This code will need to change if so.
t1 = Time('2024-04-17 22:14:45')
t2 = Time('2024-04-17 22:17:45')
evt = evt[ (Time(evt['utc']) >= t1) & (Time(evt['utc']) <= t2) ]
print(len(evt),'events between',t1,'and',t2)

# Select an energy range for the image
e_range = [4.,15.]
evt = evt[ (evt['doi_corrected_energy'] >= e_range[0]) & (evt['doi_corrected_energy'] <= e_range[1]) ]
print(len(evt),'events between',e_range,'keV')

# Make an image via a 2D histogram of event locations
# Choose the x and y ranges and the pixel (bin) size.
# Be careful about the bin size. If it's too small, you can get aliasing in the image.
# A bin size of 8 arcsec tends to work well, but this isn't a strict limit.
x_range = [-550,-250]	# arcsec
y_range = [-200,100]	# arcsec
bin_size = 8			# image pixel size in arcsec
n_bins = [int((x_range[1]-x_range[0])/bin_size),int((y_range[1]-y_range[0])/bin_size)]
img, xedges, yedges = np.histogram2d( evt['solar_x_unaligned'], evt['solar_y_unaligned'], 
                                      bins=n_bins, range=[x_range,y_range]  )

# plot the image
plt.pcolormesh(xedges, yedges, img.T, shading='auto')
plt.colorbar(label='Counts')
plt.xlabel('x [arcsec]')
plt.ylabel('y [arcsec]')
plt.title('FOXSI-4 '+name+' Image')
plt.show()
