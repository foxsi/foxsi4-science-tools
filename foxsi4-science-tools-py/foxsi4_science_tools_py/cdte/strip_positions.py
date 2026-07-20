from astropy import units as u
import numpy as np
from scipy.interpolate import interp1d

def strip_ref_micrometers(chan, allow_float=False):
    """ This function returns the strip reference point in um.

    Origin is at the center of the detector.
    The "strip reference point" is the 1D location of the center of the 
    electrode.

    For many strips (but not all), it is the midpoint between the strip 
    bin edges.

    Eventually this could be extended to handle any point, not just the 
    reference point.

    Parameters
    ----------
    chan : `int` or `numpy.array`
        The channel(s) for which the user wants the electrode center.

    allow_float : `bool`
        Set to True if a non-integer `chan` value is given (e.g., if the
        `chan` values represent energy averaged positions). This will 
        interpolate the float onto the physical size array.
        Default: False

    Returns
    -------
    `astropy.units.quantity.Quantity` (`numpy.array`-like) : 
        The electrode position(s) in um where the origin is the centre 
        of the detector.
    """
    C = np.arange(28)*100
    B = np.arange(20)*80 
    A = np.arange(16)*60  
    
    # For convenience, start with a 64-element array that will hold the position of just 
    # the positive numbers (i.e. second half of electrode array). Then, we will construct 
    # the actual array by mirroring that to the negative numbers and appending together.
    ## **STRIPS_PER_SIDE will replace 128**
    n_strips = 128
    n_half = int(n_strips/2)
    array = np.zeros(n_half, dtype=float)
    
    # Loop through the electrode locations, starting from center.
    for i in range(0,16):
        array[i] = A[i]+30 # Fine pitch strips
    for i in range(0,20):
        array[i+16] = array[16-1]+B[i]+80 # Medium pitch strips
    for i in range(0,28):
        array[i+16+20] = array[16+20-1]+C[i]+100 # Coarse pitch strips
    
    electrodes = np.concatenate((-array[::-1], array))

    if allow_float:
        i = interp1d(np.arange(len(electrodes)), 
                     electrodes, 
                     kind="linear",
                     bounds_error=True)
        return i(chan) << u.um
    
    return electrodes[chan] << u.um

def strip_ref_arcsec(chan, allow_float=False):
    """ Similar to `strip_ref_micrometers` function but units in arcsec.    
    
    Parameters
    ----------
    chan : `int` or `numpy.array`
        The channel(s) for which the user wants the electrode center.

    allow_float : `bool`
        Set to True if a non-integer `chan` value is given (e.g., if the
        `chan` values represent energy averaged positions). This will 
        interpolate the float onto the physical size array.
        Default: False

    Returns
    -------
    `astropy.units.quantity.Quantity` (`numpy.array`-like) : 
        The electrode position(s) in arcseconds where the origin is the 
        center of the detector.
    """
    # Get the position in um from the other function:
    pos_um = strip_ref_micrometers(chan, allow_float=allow_float)
    # Calculate arcsec per micron
    plate_scale = np.degrees(np.arctan(1./1000 / 2000.))*3600
    
    # Convert position to arcsec
    return (pos_um.value*plate_scale)<<u.arcsec

def strip_edges_micrometers():
    """ 
    Function to define the physical sizes of the different pitch widths 
    and what to do as they transition. 

    First edge starts at 0 um.
    
    Pitch widths are 100, 80, 60 um and the spaced between are 90 
    and 70 um (100/2+80/2 and 80/2+60/2, respectively).

    **Intended to replace `image.strip_edges_micrometers` in the pending
    "great visualization purge".

    Returns
    -------
    `astropy.units.quantity.Quantity` (`numpy.array`-like) : 
        The strip boundary edges in um where the origin is the center of 
        the detector.
    """
    ## ** Addition in future [PR#46] **
    ## STRIPS_PER_SIDE = cdte_tools_py.contextCdTeInfo["general"]["strips_per_side"] will replace "128"
    electrodes = strip_ref_micrometers(np.arange(128))
    edges_without_edges = np.mean(np.concatenate((electrodes[:-1][None,:],
                                                  electrodes[1:][None,:]), 
                                                 axis=0),
                                  axis=0)
    lower_edge = 2*edges_without_edges[0]-edges_without_edges[1]
    upper_edge = 2*edges_without_edges[-1]-edges_without_edges[-2]
    edges = np.insert(edges_without_edges, [0, len(edges_without_edges)], [lower_edge, upper_edge])
    
    return edges