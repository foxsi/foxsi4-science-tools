""" Module with utility functions for maps."""

import astropy.units as u
import numpy as np
from scipy import signal
import sunpy.map

from foxsi4_science_tools_py import IMAGE_COLORMAP, IMAGE_SCALE

def make_sunpy_map(event_list, xkey=None, ykey=None, tkey=None, number_of_bins=None):
    number_of_bins = 300 if number_of_bins is None else number_of_bins
    xkey = "solar_x" if xkey is None else xkey
    ykey = "solar_y" if ykey is None else ykey
    tkey = "utc" if tkey is None else tkey

    solx = event_list[xkey].value
    soly = event_list[ykey].value
    
    x_max = 1500
    x_min = -1500
    y_max = 1500
    y_min = -1500
    # Define the center pixel and spatial scale
    x_scale = (x_max-x_min)/number_of_bins
    y_scale = (y_max-y_min)/number_of_bins

    mid_obs_time = np.mean(event_list[tkey])

    # Define structural metadata
    dict_header = {
        "OBJECT": "Sun",
        "CTYPE1": "HPLN-TAN",  # Helioprojective longitude (typical for SunPy)
        "CTYPE2": "HPLT-TAN",  # Helioprojective latitude
        "CUNIT1": "arcsec",
        "CUNIT2": "arcsec",
        "CRPIX1": -(x_min/x_scale),
        "CRPIX2": -(y_min/y_scale),
        "CDELT1": x_scale,  # Arcseconds per pixel
        "CDELT2": y_scale,
        "DATE-OBS": str(event_list[tkey][0]), # Add a timestamp (Required for many solar coordinate transformations)
        "PIXLUNIT": "counts",
        "EXPTIME": float((event_list[tkey][-1] - event_list[tkey][0]).sec),
        "HGLT_OBS": sunpy.coordinates.sun.B0(mid_obs_time),  # get_sun_B0(mid_obs_time),
        "HGLN_OBS": 0,
        "RSUN_OBS": sunpy.coordinates.sun.angular_radius(
            mid_obs_time
        ).value,  
        "RSUN_REF": sunpy.sun.constants.radius.value,  
        "DSUN_OBS": sunpy.coordinates.sun.earth_distance(mid_obs_time).value * u.astrophys.au,
    }
    header = sunpy.util.MetaDict(dict_header)
    
    H, _, _ = np.histogram2d(soly, solx, bins=number_of_bins, range=[[y_min, y_max], [x_min, x_max]])
    return assign_plot_settings(sunpy.map.Map(H, header))

def assign_plot_settings(sunpy_map):
    """Function to set plot settings that are multiple lines."""
    sunpy_map.plot_settings["norm"] = IMAGE_SCALE
    sunpy_map.plot_settings["cmap"] = IMAGE_COLORMAP
    return sunpy_map

def draw_grid(sunpy_map_obj, axes):
    """Draw a grid representing the Sun."""
    sunpy_map_obj.draw_limb(
        axes=axes, color="black", linewidth=1, linestyle="dashed", zorder=0
    )
    sunpy_map_obj.draw_grid(
        axes=axes,
        color="grey",
        alpha=0.5,
        linewidth=1,
        linestyle="dashed",
        zorder=0,
        annotate=False,
    )

def corrolation_map(reference_map, new_map, **kwargs):
    defaults = {"boundary":"symm", "mode":"same"} | kwargs
    # check image scales
    if not all(np.isclose([reference_map.meta["CDELT1"]], [new_map.meta["CDELT1"]])):
        print("Image scales are different. No correlation being done.")
        return
    reference_array, new_array = reference_map.data, new_map.data
    corr = signal.correlate2d(reference_array, new_array, **defaults)
    return assign_plot_settings(sunpy.map.Map(corr, reference_map.meta))

def correlation_map_shift(reference_map, new_map, **kwargs):
    cor_map = corrolation_map(reference_map, new_map, **kwargs)
    pixel_pos = np.argwhere(cor_map.data == cor_map.data.max()) * u.pixel
    hpc_max = cor_map.wcs.pixel_to_world(pixel_pos[:, 1], pixel_pos[:, 0])
    return hpc_max.Tx-reference_map.center.Tx, hpc_max.Ty-reference_map.center.Ty
