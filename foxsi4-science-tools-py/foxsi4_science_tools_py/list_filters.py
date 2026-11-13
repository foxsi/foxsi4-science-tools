import astropy.units as u
import numpy as np
from shapely import points

__all__ = ["by_general_region", "by_signal", "by_time"]

def by_general_region(event_data, shape):
    """Filter and return event list filtered for a general region.

    Parameters
    ----------
    event_data: FITS data structure
        This should be an hdu.data structure from a NuSTAR FITS file.

    shape : `~shapely.geometry.polygon.Polygon`
    """
    # enforce everything in arcseconds before operations
    sol_x_arc = event_data["solar_x"] << u.arcsec
    sol_y_arc = event_data["solar_y"] << u.arcsec
    
    map_points = points(np.vstack((sol_x_arc.value, sol_y_arc.value)).T)
    in_shape = np.array([(mp.within(shape) or mp.touches(shape)) for mp in map_points])
    return event_data[in_shape]

def by_signal(event_data, signal_low, signal_high):
    """Apply signal filtering to the data.

    Parameters
    ----------
    event_data: FITS data class
        The event list.

    signal_low: float
        Low-side signal bound for the map you want to produce in 
        relevant units. Inclusive.

    signal_high: float
        High-side signal bound for the map you want to produce in 
        relevant units. Exclusive.
    """
    signal = event_data["signal"] 

    signal_low <<= signal.unit
    signal_high <<= signal.unit
    
    sig_filter = (signal_low <= signal) & (signal < signal_high)
    goodinds = (sig_filter).nonzero()[0]

    return event_data[goodinds]

def by_time(event_data, start, end):
    """Only include counts within a given time range.

    Parameters
    ----------
    event_data: 
        The event list.

    start : `str`
        Input two times in the form 'yyyy/mm/dd, HH:MM:SS'
        (e.g. '2019/03/06, 16:45:30'). First time is inclusive and 
        second time is exclusive.
    """

    tstart = Time(start, format="isot", scale="utc")
    tend = Time(end, format="isot", scale="utc")
    times = event_data["utc"]

    time_filter = (tstart <= times) & (times < tend)
    goodinds = (time_filter).nonzero()[0]

    return event_data[goodinds]