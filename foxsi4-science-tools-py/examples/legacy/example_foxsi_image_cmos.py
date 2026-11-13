"""
This code provides a simple example of making an image from FOXSI-4 
CMOS data, using the May 2025 data provided by NAOJ.

Please note that this example will require updating as the data product
advances.

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

import foxsi4_science_tools_py
import foxsi4_science_tools_py.util.solar_coords_cmos as sol

# Choose file and directory
# File should be a CMOS event list. (Example: 'foxsi4_telescope-0_CMOS-1_PHOTON_EVENT_V25MAY21.fits')
# Example here chooses a CMOS1 file. Also define a string with this information.
# Don't forget to change the directory path to match YOUR path to the data.
dir = "/Users/glesener/Dropbox/data/foxsi/20240417/lev0_data_may2025/"
file = 'foxsi4_telescope-1_CMOS-2_PHOTON_EVENT_V25MAY21.fits'
name = 'CMOS2'

# Read in event data from the file
with fits.open( dir + file ) as hdul:
    hdus = deepcopy(hdul)

evt = Table()
for hdu in hdus[1:]:
    evt.add_column( hdu.data, name=hdu.name )

# EVT now holds the event list for the selected detector.
print(len(evt),'events read from file for '+name)

# Apply some processing to get CMOS info into solar coordinates and UTC.

# Transfer pixel numbers to solar coordinates with center of detector at solar center.
evt = sol.compute_solar_coords( evt )

# Convert linetime to UTC
if name == 'CMOS1':
    evt = sol.compute_utc( evt, 0 )		# Second argument is telescope 0 or 1 (Tel0=CMOS1)
elif name == 'CMOS2':
    evt = sol.compute_utc( evt, 1 )		# Second argument is telescope 0 or 1 (Tel0=CMOS1)
else:
    print("NAME needs to be CMOS1 or CMOS2")

# Apply SPARCS pointing info - you'll need to specify where your pointing file is.
# This gets you coarsely in the right place but does not account for 
# uncalibrated SPARCS-EXP offsets or the observed drift.
sparcs_file = '~/Dropbox/foxsi/foxsi4/science-analysis/analysis-code/foxsi4-science-tools/foxsi4-science-tools-py/foxsi4_science_tools_py/util/pointing_v2.csv'
evt = sol.apply_sparcs_pointing( evt, sparcs_file=sparcs_file )
evt.info()

# Select a time range for the image
# Note that we might change the time format of the event list to a different format.
# This code will need to change if so.
# This snippet tries to handle either Time or float formats.
t1 = Time('2024-04-17 22:14:45')
t2 = Time('2024-04-17 22:17:45')
if isinstance( evt['utc'], Time ):
    evt = evt[ (evt['utc'].unix >= t1.unix) & (evt['utc'].unix <= t2.unix) ]
else:
    evt = evt[ (evt['utc'] >= t1.unix) & (evt['utc'] <= t2.unix) ]
print(len(evt),'events between',t1,'and',t2)

"""
# Select an energy range for the image
# Since we are waiting on CMOS energy calibration, there's not a good way to do this yet.
e_range = [4.,15.]
evt = evt[ (evt['doi_corrected_energy'] >= e_range[0]) & (evt['doi_corrected_energy'] <= e_range[1]) ]
print(len(evt),'events between',e_range,'keV')
"""

# Make an image via a 2D histogram of event locations
# Choose the x and y ranges and the pixel (bin) size.
# Be careful about the bin size. If it's too small, you can get aliasing in the image.
# A bin size of 8 arcsec tends to work well, but this isn't a strict limit.
x_range = [-800, 0]	# arcsec
y_range = [-600,200]	# arcsec
bin_size = 9			# image pixel size in arcsec
n_bins = [int((x_range[1]-x_range[0])/bin_size),int((y_range[1]-y_range[0])/bin_size)]
img, xedges, yedges = np.histogram2d( evt['solar_x_unaligned'], evt['solar_y_unaligned'], bins=n_bins, 
                                      range=[x_range,y_range]  )

# plot the image
plt.pcolormesh(xedges, yedges, img.T, shading='auto')
plt.colorbar(label='Counts')
plt.xlabel('x [arcsec]')
plt.ylabel('y [arcsec]')
plt.title('FOXSI-4 '+name+' Image')
plt.show()
