"""Core functions for worjing with FOXSI-4/5 data."""

from copy import deepcopy

import astropy.units as u
import numpy as np
import sunpy.map
from astropy.coordinates import SkyCoord
from astropy.time import Time

from foxsi4_science_tools_py.io.detectors import load_cdte, load_cmos
from foxsi4_science_tools_py.util.clean import cdte_usual_suspects, cmos_usual_suspects
from foxsi4_science_tools_py.util.maps import make_sunpy_map, assign_plot_settings, correlation_map_shift


def cdte_events(cdte_file, name):
    """Given a file, load in the CdTe data."""
    cdte_data = load_cdte(cdte_file)[1]
    return cdte_usual_suspects(cdte_data, name=name)

def cmos_events(cmos_file, name):
    """Given a file, load in the CdTe data."""
    cmos_data = load_cmos(cmos_file)[1]
    return cmos_usual_suspects(cmos_data, name=name)

def full_disk_map(event_data, **kwargs):
    """Return a Sunpy map of the FOXSI data."""
    bottom_left = [-1200, -1200] << u.Unit("arcsec")
    top_right = [1200, 1200] << u.Unit("arcsec")
    foxsi_map = make_sunpy_map(event_data, **kwargs)
    return get_submap(foxsi_map, bottom_left, top_right)

def disk_map(event_data, bottom_left, top_right, **kwargs):
    """Return a Sunpy map of the FOXSI data."""
    foxsi_map = make_sunpy_map(event_data, **kwargs)
    return get_submap(foxsi_map, bottom_left, top_right)

def get_submap(sunpy_map_obj, bottom_left: u.Quantity, top_right: u.Quantity):
    """Create and return a submap of a Sunpy map object.

    Parameters
    ----------
    sunpy_map_obj : `~sunpy.map.Map`
        The map object.

    bottom_left : `~astropy.units.Quantity`
        The bottom left coordinates of the submap. Unit aware (e.g.,
        arc-seconds)

    top_right : `~astropy.units.Quantity`
        The top right coordinates of the submap. Unit aware (e.g.,
        arc-seconds)
    """
    # fix really large plot, instead of going from -3600 to 3600 in x and y
    bl = SkyCoord(*bottom_left, frame=sunpy_map_obj.coordinate_frame)
    tr = SkyCoord(*top_right, frame=sunpy_map_obj.coordinate_frame)
    _submap = sunpy_map_obj.submap(
        bl, top_right=tr
    )  # submaps act differently to normal maps
    return assign_plot_settings(sunpy.map.Map(_submap.data, _submap.meta))

def event_data_xy_shift(
    event_data, dx: u.Quantity | None = None, dy: u.Quantity | None = None
):
    """Apply a shift in the positions of each event.

    Parameters
    ----------
    event_data :
        The event list.

    dx : `u.Quantity`|`None`
        The change in the x-position of the counts. Needs to be
        convertible to arcseconds or will be assumed to be arcseconds
        if no units exist.
        Default: None

    dy : `u.Quantity`|`None`
        The change in the y-position of the counts. Needs to be
        convertible to arcseconds or will be assumed to be arcseconds
        if no units exist.
        Default: None
    """
    _event_data = deepcopy(event_data)
    if dx is not None:
        _event_data.solar_x += (dx << u.arcsec)
    if dy is not None:
        _event_data.solar_y += (dy << u.arcsec)
    return _event_data


def coalign_event_lists(reference_list, new_list, map_kwargs=None, correlation_kwargs=None):
    """Returns altered new list so that it is coaligned with the reference."""
    map_kwargs = {} if map_kwargs is None else map_kwargs
    _map_kwargs = {"bottom_left":[-1200, -1200]<<u.arcsec, 
                   "top_right":[-1200, -1200]<<u.arcsec, 
                   "number_of_bins":300} | map_kwargs
    correlation_kwargs = {} if correlation_kwargs is None else correlation_kwargs
    reference_map = disk_map(reference_list, **_map_kwargs)
    new_map = disk_map(new_list, **map_kwargs)
    shift_x, shift_y = correlation_map_shift(reference_map, new_map, **correlation_kwargs)
    return event_data_xy_shift(new_list, dx=shift_x, dy=shift_y), (shift_x, shift_y)

def count_time_profile_array(
    event_data=None,
    time_bins=None,
    time_binning=None,
    start_time=None,
    end_time=None,
):
    """Get counts and time bins for plotting a NuSTAR time profile.

    Parameters
    ----------
    event_data :
        The event list. If `None` then `self.cleaned_evt_data` will
        be used.
        Default: None

    time_bins :
        Takes priority over `time_binning`, `start_time`, and
        `stop_time`.
        Default: None

    time_binning :
        If `None` then a default value is used.
        Default: None

    start_time :
        If `None` then a default value is used.
        Default: None

    end_time :
        If `None` then a default value is used.
        Default: None

    Returns
    -------
    """
    if time_bins is None:
        return count_time_profile_array_from_start_stop_binning(
            event_data=event_data,
            time_binning=time_binning,
            start_time=start_time,
            end_time=end_time,
        )
    return count_time_profile_array_from_binning_array(
        event_data=event_data, time_bins=time_bins
    )

def count_time_profile_array_from_start_stop_binning(
    event_data=None, time_binning=None, start_time=None, end_time=None
):
    """Get counts and time bins for plotting a NuSTAR time profile.

    User given start, stop, and/or time bin size.

    Parameters
    ----------
    event_data :
        The event list. If `None` then `self.cleaned_evt_data` will
        be used.
        Default: None

    time_binning :

        Default: None

    start_time :

        Default: None

    end_time :

        Default: None

    Returns
    -------
    """
    time_binning = 10 << u.second if time_binning is None else time_binning << u.second
    start_time = np.min(event_data["utc"]) if start_time is None else start_time
    end_time = np.max(event_data["utc"]) if end_time is None else end_time

    _st = 0 << u.second
    _et = (end_time-start_time).sec << u.second

    time_bins = np.arange(
        _st.value, _et.value + time_binning.value, time_binning.value
    )
    count_times = (event_data["utc"]-start_time).sec << u.second
    counts, time_bins = np.histogram(count_times.value, time_bins)
    return start_time+(time_bins<<u.second), counts << u.ct

def count_time_profile_array_from_binning_array(
    event_data=None, time_bins=None
):
    """Get counts and time bins for plotting a NuSTAR time profile.

    User given time bins.

    Parameters
    ----------
    event_data :
        The event list. If `None` then `self.cleaned_evt_data` will
        be used.
        Default: None

    time_bins :

        Default: None

    Returns
    -------
    """
    start_time = time_bins[0]
    time_bins_s = np.array((time_bins-start_time).sec) << u.second
    count_times = np.array((event_data["utc"]-start_time).sec) << u.second
    counts, time_bins = np.histogram(count_times.value, time_bins_s.value)
    return start_time+(time_bins<<u.second), counts << u.ct