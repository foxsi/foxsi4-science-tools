"""
Series of functions to handle the conversion from CMOS pixel to solar coordinates.

"""

from astropy import units as u
from astropy.io import ascii
from astropy.time import Time, TimeDelta
from astropy.io import fits
from astropy.table import Table
import numpy as np
import keyword
from scipy.interpolate import interp1d
from scipy.spatial.transform import Rotation as R

#from cdte_tools_py.util.loops import progress_printer
#import cdte_tools_py.visualization.image as image

UNPHYSICAL_SOLAR_COORD = np.nan
# Note: have not implemented this anywhere; revisit if needed.

def pixel_micrometers(pix, dim ):
    """ This function returns the pixel position on the detector in um.
    Origin is at the center of the detector, defined as the corner of pixels (1023,1024) 
    in the X direction and (959,960) in Y. In other words, the origin is not the center 
    of any pixel because the number of pixels is even in both directions.

    Parameters
    ----------
    pix : `int` or `numpy.array`
        The pixel(s) for which the user wants the position.
        Integer values refer to the center of the pixel.
        Non-integer values are allowed.
    dim : `int`
    	0 or 1, defines whether the coordinate is an X coord (0) or Y coord (1).

    Returns
    -------
    `astropy.units.quantity.Quantity` (`numpy.array`-like) : 
        The pixel position(s) in um where the origin is the centre 
        of the detector.
    """

    # Check if dim is a sensible value (0 or 1)
    if dim != 0 and dim != 1:
        print("Parameter DIM needs to be 0 or 1.")
        return np.nan
    # Note: here are a bunch of parameters that should go in a config file.
    cmos_pitch = 11 #11 um per pixel
    n_pix_X = 2048	# number of pixels in X direction
    n_pix_Y = 1920	# number of pixels in Y direction
    pix0_X = 0		# actually, values only occupy pixels 1-2046 but this is probably just excluding edge pixels.
    pix0_Y = 64		# pixels occupy values 65-1982
    
    # calculate distance from detector center along respective dimension, in pixel units
    if dim==0:
        pix_from_center = pix-pix0_X - (n_pix_X/2-0.5)
    if dim==1:
        pix_from_center = pix-pix0_Y - (n_pix_Y/2-0.5)

    return (pix_from_center*cmos_pitch) << u.um

def pixel_arcsec(pix, dim):
    """ Similar to `pixel_micrometers` function but units in arcsec.    
    
    Parameters
    ----------
    pix : `int` or `numpy.array`
        The pixel(s) for which the user wants the position.
        Integer values refer to the center of the pixel.
        Non-integer values are allowed.
    dim : `int`
    	0 or 1, defines whether the coordinate is an X coord (0) or Y coord (1).

    Returns
    -------
    `astropy.units.quantity.Quantity` (`numpy.array`-like) : 
        The pixel position(s) in arcseconds where the origin is the 
        center of the detector.
    """
    # Get the position in um from the other function:
    pos_um = pixel_micrometers( pix, dim )
    # Calculate arcsec per micron
    plate_scale = np.degrees(np.arctan(1./1000 / 2000.))*3600
    
    # Convert position to arcsec
    return (pos_um.value*plate_scale)<<u.arcsec

def read_cmos_event_list( file ):
    """ Reads in a FITS file containing a CMOS event list and produces an
    astropy table with the event data.
    
    Example:
    evt1 = read_cmos_event_list( "lev0_data_may2025/foxsi4_telescope-0_CMOS-1_PHOTON_EVENT_V25MAY21.fits" )
    evt2 = read_cmos_event_list( "lev0_data_may2025/foxsi4_telescope-1_CMOS-2_PHOTON_EVENT_V25MAY21.fits" )
    
    """
    # Open the FITS file
    hdulist = fits.open( file )
    
    evt = Table()
    for hdu in hdulist[1:]:
        evt.add_column( hdu.data, name=hdu.name )
        
    hdulist.close()
    return evt

def compute_solar_coords( evt ):
    """ Computes a position on the Sun for each photon in an event list.
    
    It does not incorporate SPARCS information, so it puts the center of 
    the detector at solar center.

    For now, assuming no rotation is needed; REVISIT this later.
    
    No filtering is done for event quality. We don't judge here.

    Parameters
    ----------
    evt : `astropy.table.table.Table`
        The CMOS event list.

    Returns
    -------
    `astropy.table.table.Table` :
        The original event list but now with the added fields: 
        - "solar_x_unaligned": 
        - "solar_y_unaligned": same as "solar_x_unaligned"" but for y-axis
    """
    evt_copy = evt.copy()
    
    ## convert detector position to arcsec on the Sun.
    _solar_x_unaligned = pixel_arcsec( evt_copy['PIXEL_X'], 0 )
    _solar_y_unaligned = pixel_arcsec( evt_copy['PIXEL_Y'], 1 )
    
    ## Add new columns to the EVT data structure.
    ## Note: y axis points toward Payload 0.
    evt_copy.add_column(_solar_x_unaligned, name="solar_x_unaligned") # X coord in arcsec
    evt_copy.add_column(_solar_y_unaligned, name="solar_y_unaligned") # Y coord in arcsec

    return evt_copy

def compute_utc( evt, telescope_no, cmos=None ):

    """ Convert linetime (CMOS time counter) to UTC (coordinated universal time)
    
    This code is based on the IDL code convert_linetime_utc by NAOJ, specifically
    ver FOXSI4_CMOS_codes_20250410T07
    
    CALLING SEQUENCE
    ----------------
    
        compute_utc( evt, telescope_no, [CMOS = CMOS] )
    
    INPUTS
    ------
    
        linetime - (long or long array) Linetime value(s) output from CMOS camera
        telescope_no - (int) Telescope No. : 0 (Marshall optics) or 1 (Nagoya optics)
        
    KEYWORDS
    --------
    
        CMOS   - [Optional] (int) If set as "CMOS = 1", use CMOS1 conversion
                                  If set as "CMOS = 2", use CMOS2 conversion
               - If this keyword is used, telescope_no can be omitted.
       
    Returns
    -------
    `astropy.table.table.Table` :
        The event list with a new field for UTC time ("utc") in float format
        
    """
    evt_copy = evt.copy()
    
    # Note: here are more parameters that should go in a config file.

    #reference absolute time determined by the timing of LISS pitch change
    time0 = Time('2024-04-17 22:18:19.902', format="iso", scale="utc")  # time of LISS change
    linetime0_CMOS1  = 53312257              # reference linetime derived from QL data of CMOS1
    offset_for_CMOS2_linetime = 37898995     # difference in linetime between CMOS1 and CMOS2 derived from QL data
    clock_hz = 25.e6              # CMOS camera operating clock [sec^-1]
    T_clk_pix = 1 / clock_hz    # [sec / clock]
    T_line = 513 * T_clk_pix    # [sec / linetime]
    
    # First, convert linetime to seconds (don't worry about reference time yet).
    # This is the same for both telescopes.
    ####linetime_sec = evt_copy["LINETIME"] * T_line  ; [sec]
    
    if keyword.iskeyword( cmos ):
        if cmos==1:
            telescope_no = 0
        if cmos==0:
            telescope_no = 1
    
    if telescope_no == 0:
        new_time = evt_copy["LINETIME"]
    elif telescope_no == 1:
    	new_time = evt_copy["LINETIME"] + offset_for_CMOS2_linetime
    else:
        print("Error: telescope_no shall be 0 or 1, or Specify CMOS=1 or CMOS=2")
        return -1
        
    utc = time0 + TimeDelta( T_line * (new_time.data - linetime0_CMOS1), format="sec" )

    evt_copy.add_column(utc, name="utc")

    return evt_copy

def apply_sparcs_pointing(evt, sparcs_file=None):
    """ Reads in the SPARCS pointing information and offsets the EVT 
    solar coordinates to approximately the right place on the Sun.

    The known SPARCS-EXP offset is included.

    Not included are unknown SPARCS-EXP additional offsets, imperfections 
    in detector placement, strange solar drift, etc.

    To REALLY get the right flare position you need to coalign with AIA, 
    and that is not included in this function.
    
    Parameters
    ----------
    evt : `astropy.table.table.Table`
        The event list with fields "solar_x_unaligned" and
        "solar_y_unaligned".

    sparcs_file : `str`
        The name of the SPARCS pointing file. Passed to `ascii.read`.
        Nominally: `os.path.join(pathlib.Path(__file__).parent, 
                                 "solar_coords_data", 
                                 f"pointing_{_spa_fver}.csv")`
    
    Returns
    -------
    `astropy.table.table.Table` :
        The event list with "solar_x_unaligned" and "solar_y_unaligned"
        adjusted for SPARCS pointing.
        
    Important note! This function is the same for CdTe and CMOS, so in principle it 
    doesn't need to exist twice.
        
    """
    evt_copy = evt.copy()

    ## ** Addition in future PR **
    ## # if a sparcs file is given, use it, if not, use the default
    ## _spa_fver = cdte_tools_py.contextCdTeInfo.get("paths")["path2"]["spa_fver"]
    ## sparcs_file = os.path.join(pathlib.Path(__file__).parent, 
    ##                                 "solar_coords_data", 
    ##                                 f"pointing_{_spa_fver}.csv") if sparcs_file is None else sparcs_file
    
    # Read in SPARCS pointing data file.
    sparcs = ascii.read(sparcs_file)
    sparcs_times = Time(sparcs["Command initiation UTC"])
    # SPARCS pitch = arcsec from solar center WESTward
    sparcs_x = np.array(sparcs["Pitch (\'\')"],dtype=float)
    # SPARCS yaw = arcsec from solar center SOUTHward
    sparcs_y = np.array(-sparcs["Yaw (\'\')"], dtype=float)
    
    # Apply known SPARCS - EXP offset.
    # This is hard-coded because it shouldn't change; it is the offset that we 
    # actually used at PFRR when coordinating with SPARCS personnel.
    offset = (149.21, -156.15) # This is SPARCS pos minus FOXSI pos
    sparcs_x -= offset[0]
    sparcs_y -= offset[1]
    
    ## Compare times to see which pointing interval each photon is in.
    evt_times_utc = Time(evt_copy["utc"], format="unix_tai", scale="utc")
    for i in range(0,len(sparcs_times)-1):
        ind = np.nonzero((evt_times_utc >= sparcs_times[i]) & (evt_times_utc < sparcs_times[i+1]))
        evt_copy["solar_x_unaligned"][ind] += sparcs_x[i]
        evt_copy["solar_y_unaligned"][ind] += sparcs_y[i]
    
    return evt_copy