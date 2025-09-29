;
; Make images for FOXSI-4 data
; This is IDL CODE
;

;; Make sure to change the directory to wherever your data are.
dir = '/Users/glesener/Dropbox/data/foxsi/20240417/cdte_lev3/'
filename = 'cdte1_level3_version1.fits'
file = file_search( dir+filename )
name = 'CdTe1'

evt = mrdfits( file, 1, hdr )
; EVT now holds the event list for the selected detector.
print, n_elements(evt), ' events read from file for '+name

; Restrict to events with nominal quality flag
evt = evt[ where(evt.flag_quality eq 1) ]
print, n_elements(evt), ' perfect events'

; Select a time range for the image
; Note that we might change the time format of the event list to a different format.
; This code will need to change if so.
t1 = '2024-04-17 22:14:45'
t2 = '2024-04-17 22:17:45'
t1 = anytim(t1)
t2 = anytim(t2)
evt = evt[ where( anytim(evt.utc,fiducial='sys') ge t1 and anytim(evt.utc,fiducial='sys') le t2 ) ]
print, n_elements(evt), ' events between ', anytim(t1,/yo), ' and ', anytim(t2,/yo)

; Select an energy range for the image
e_range = [4.,15.]
evt = evt[ where( evt.doi_corrected_energy ge e_range[0] and evt.doi_corrected_energy le e_range[1] ) ]
print, n_elements(evt), ' events in keV range ', e_range

; Make an image via a 2D histogram of event locations
; Choose the x and y ranges and the pixel (bin) size.
; Be careful about the bin size. If it's too small, you can get aliasing in the image.
; A bin size of 8 arcsec tends to work well, but this isn't a strict limit.
x_range = [-550,-250]	; arcsec
y_range = [-200,100]	; arcsec
bin_size = 8			; image pixel size in arcsec
n_bins = [fix((x_range[1]-x_range[0])/bin_size),fix((y_range[1]-y_range[0])/bin_size)]

img = hist_2d( evt.solar_x_unaligned, evt.solar_y_unaligned, $
			   min1=x_range[0], max1=x_range[1], min2=y_range[0], max2=y_range[1], $
			   bin1=bin_size, bin2=bin_size )

map = make_map( img, xcen=mean(x_range), ycen=mean(y_range), dx=bin_size, dy=bin_size, $
	time=anytim( mean(anytim(evt.utc,fid='sys')),/yo), dur=t2-t1 )
plot_map, map
