"""Functions to perform common filtering/conversionson/extract CdTe data. 

Filters
-------

Filters return an Astropy table where triggers have been filter out.
    - E.g., filter out triggers that are not single hit detections.

All filter and conversion functions should state what CdTe data for which
they are valid.

The filters can be chained together. For example, to filter triggers that 
have energies on both detector sides and only came from the Telemetry 
data source:

>>> from astropy.table import Table
>>> from cdte_tools_py.calibration.path1 import get_path1
>>> from cdte_tools_py.io.fits_tools import load_fits

>>> filename = f"path/to/file/cdte1_merged_v0_0_3.fits"
>>> existing_cdte_hdus = load_fits(filename)

>>> path1_cdte_hdus = get_path1(existing_cdte_hdus)

>>> cdte_data = Table(path1_cdte_hdus[1].data)
>>> filtered_cdte_data = both_sides_hit_filter(telemetry_filter(cdte_data))

Converters
----------

Converters return an Astropy table with field formats edited.
    - E.g., when only looking at single Pt-side hits, replace the 
      "pt_merged_energy_list" with one where the dummy value has been 
      removed making it converting it from an array of shape (N, 64) to
      a shape of (N, 1) where N is the number of triggers.

Extractors
----------

Extractors will return data from the Astropy table in an unusual form.
    - E.g., arrays need to be rectangular with every entry having the 
      same number of values; however, this might not always be wanted.
      A user might want to see the strip/energy clumps in a list format
      where each trigger entry is a non-uniform list of grouped strips
      and energies.
"""

import numpy as np

UNPHYSICAL_ENERGY = -1.0

# Filters
# -------
def telemetry_filter(cdte_data):
    """Leave only the CdTe triggers that came from the Telemetry source.

    Filter is valid for any merged data product where the table has a
    field called "source" with potential values "telemetry", "common",
    and "de".

    Parameters
    ----------
    cdte_data : `astropy.table.Table`
        The Astropy table of the CdTe data with energy information.

    Returns
    -------
    `astropy.table.Table` : 
        The CdTe data with only the triggers that appeared in the Telemetry 
        data source.

    Examples
    --------
    >>> from astropy.table import Table
    >>> from cdte_tools_py.calibration.path1 import get_path1
    >>> from cdte_tools_py.io.fits_tools import load_fits

    >>> filename = f"path/to/file/cdte1_merged_v0_0_3.fits"
    >>> existing_cdte_hdus = load_fits(filename)

    >>> path1_cdte_hdus = get_path1(existing_cdte_hdus)

    >>> cdte_data = Table(path1_cdte_hdus[1].data)
    >>> filtered_cdte_data = telemetry_filter(cdte_data)
    """
    bool_tel_data = (cdte_data["source"]=="telemetry") | (cdte_data["source"]=="common")
    return cdte_data[bool_tel_data]

def exclusive_telemetry_filter(cdte_data):
    """Leave only the CdTe triggers that came **exclusively** from the Telemetry source.

    Filter is valid for any merged data product where the table has a
    field called "source" with potential values "telemetry", "common",
    and "de".

    Parameters
    ----------
    cdte_data : `astropy.table.Table`
        The Astropy table of the CdTe data with energy information.

    Returns
    -------
    `astropy.table.Table` : 
        The CdTe data with only the triggers that only appeared in the 
        Telemetry data source.

    Examples
    --------
    >>> from astropy.table import Table
    >>> from cdte_tools_py.calibration.path1 import get_path1
    >>> from cdte_tools_py.io.fits_tools import load_fits

    >>> filename = f"path/to/file/cdte1_merged_v0_0_3.fits"
    >>> existing_cdte_hdus = load_fits(filename)

    >>> path1_cdte_hdus = get_path1(existing_cdte_hdus)

    >>> cdte_data = Table(path1_cdte_hdus[1].data)
    >>> filtered_cdte_data = exclusive_telemetry_filter(cdte_data)
    """
    bool_extel_data = (cdte_data["source"]=="telemetry")
    return cdte_data[bool_extel_data]

def de_filter(cdte_data):
    """Leave only the CdTe triggers that came from the DE source.

    Filter is valid for any merged data product where the table has a
    field called "source" with potential values "telemetry", "common",
    and "de".

    Parameters
    ----------
    cdte_data : `astropy.table.Table`
        The Astropy table of the CdTe data with energy information.

    Returns
    -------
    `astropy.table.Table` : 
        The CdTe data with only the triggers that appeared in the DE 
        data source.

    Examples
    --------
    >>> from astropy.table import Table
    >>> from cdte_tools_py.calibration.path1 import get_path1
    >>> from cdte_tools_py.io.fits_tools import load_fits

    >>> filename = f"path/to/file/cdte1_merged_v0_0_3.fits"
    >>> existing_cdte_hdus = load_fits(filename)

    >>> path1_cdte_hdus = get_path1(existing_cdte_hdus)

    >>> cdte_data = Table(path1_cdte_hdus[1].data)
    >>> filtered_cdte_data = de_filter(cdte_data)
    """
    bool_de_data = (cdte_data["source"]=="de") | (cdte_data["source"]=="common")
    return cdte_data[bool_de_data]

def exclusive_de_filter(cdte_data):
    """Leave only the CdTe triggers that came **exclusively** from the DE source.

    Filter is valid for any merged data product where the table has a
    field called "source" with potential values "telemetry", "common",
    and "de".

    Parameters
    ----------
    cdte_data : `astropy.table.Table`
        The Astropy table of the CdTe data with energy information.

    Returns
    -------
    `astropy.table.Table` : 
        The CdTe data with only the triggers that only appeared in the DE 
        data source.

    Examples
    --------
    >>> from astropy.table import Table
    >>> from cdte_tools_py.calibration.path1 import get_path1
    >>> from cdte_tools_py.io.fits_tools import load_fits

    >>> filename = f"path/to/file/cdte1_merged_v0_0_3.fits"
    >>> existing_cdte_hdus = load_fits(filename)

    >>> path1_cdte_hdus = get_path1(existing_cdte_hdus)

    >>> cdte_data = Table(path1_cdte_hdus[1].data)
    >>> filtered_cdte_data = exclusive_de_filter(cdte_data)
    """
    bool_exde_data = (cdte_data["source"]=="de")
    return cdte_data[bool_exde_data]

def common_filter(cdte_data):
    """Leave only the CdTe triggers from both the DE and Telemetry source.

    Filter is valid for any merged data product where the table has a
    field called "source" with potential values "telemetry", "common",
    and "de".

    Parameters
    ----------
    cdte_data : `astropy.table.Table`
        The Astropy table of the CdTe data with energy information.

    Returns
    -------
    `astropy.table.Table` : 
        The CdTe data with only the triggers that appeared in both the 
        Telemetry and DE  data source.

    Examples
    --------
    >>> from astropy.table import Table
    >>> from cdte_tools_py.calibration.path1 import get_path1
    >>> from cdte_tools_py.io.fits_tools import load_fits

    >>> filename = f"path/to/file/cdte1_merged_v0_0_3.fits"
    >>> existing_cdte_hdus = load_fits(filename)

    >>> path1_cdte_hdus = get_path1(existing_cdte_hdus)

    >>> cdte_data = Table(path1_cdte_hdus[1].data)
    >>> filtered_cdte_data = common_filter(cdte_data)
    """
    bool_exde_data = (cdte_data["source"]=="common")
    return cdte_data[bool_exde_data]

def both_sides_hit_filter(cdte_data):
    """Leave only the CdTe triggers that contain hits on both detector sides. 

    Filter is valid for CdTe data that has already had its ADC values 
    converted to energy (e.g., `path2` data from
    cdte_tools_py.calibration.path2.get_path2).
    
    This filter mimics the filter performed for `path2` in the C++ code
    in `path2_cal_and_merge_eachside.cpp`.

    Parameters
    ----------
    cdte_data : `astropy.table.Table`
        The Astropy table of the CdTe data with energy information.

    Returns
    -------
    `astropy.table.Table` : 
        The CdTe data with only the triggers that have strips above the 
        energy thresholds used on both the Pt- and Al-side of the detector.
        - Filters out pseudo triggers.
        - Filters out averaged strips in the guard ring of the detector.
            - Edge strips of values <4 and >122.

    Examples
    --------
    >>> from astropy.table import Table
    >>> from cdte_tools_py.calibration.path1 import get_path1
    >>> from cdte_tools_py.calibration.path2 import get_path2
    >>> from cdte_tools_py.io.fits_tools import load_fits

    >>> filename = f"path/to/file/cdte1_merged_v0_0_3.fits"
    >>> existing_cdte_hdus = load_fits(filename)

    >>> path1_cdte_hdus = get_path1(existing_cdte_hdus)

    >>> path2_cdte_hdus = get_path2(path1_cdte_hdus)

    >>> cdte_data = Table(path2_cdte_hdus[1].data)
    >>> filtered_cdte_data = both_sides_hit_filter(cdte_data)
    """
    
    cdte_data = _single_hit_no_guard_ring_filter(cdte_data)

    valid_trig = (cdte_data["pt_merged_nhit"]>=1)\
                 &(cdte_data["al_merged_nhit"]>=1)\
                 &(cdte_data["flag_pseudo"]!=1)
    
    return cdte_data[valid_trig]

def no_pseudo_filter(cdte_data):
    """Leave only real (non-pseudo) triggers. 
    
    Filter is valid for CdTe data that comes out of the parser, so the 
    earliest human readable data. Only requires the field "flag_pseudo".

    Parameters
    ----------
    cdte_data : `astropy.table.Table`
        The Astropy table of the CdTe data with energy information.

    Returns
    -------
    `astropy.table.Table` : 
        The CdTe data with only real triggers.

    Examples
    --------
    >>> from astropy.table import Table
    >>> from cdte_tools_py.calibration.path1 import get_path1
    >>> from cdte_tools_py.io.fits_tools import load_fits

    >>> filename = f"path/to/file/cdte1_merged_v0_0_3.fits"
    >>> existing_cdte_hdus = load_fits(filename)

    >>> path1_cdte_hdus = get_path1(existing_cdte_hdus)

    >>> cdte_data = Table(path1_cdte_hdus[1].data)
    >>> filtered_cdte_data = no_pseudo_filter(cdte_data)
    """
    valid_trig = (cdte_data["flag_pseudo"]==0)
    return cdte_data[valid_trig]

def pseudo_filter(cdte_data):
    """Leave only pseudo (injected) triggers. 
    
    Filter is valid for CdTe data that comes out of the parser, so the 
    earliest human readable data. Only requires the field "flag_pseudo".

    Parameters
    ----------
    cdte_data : `astropy.table.Table`
        The Astropy table of the CdTe data with energy information.

    Returns
    -------
    `astropy.table.Table` : 
        The CdTe data with only pseudo triggers.

    Examples
    --------
    >>> from astropy.table import Table
    >>> from cdte_tools_py.calibration.path1 import get_path1
    >>> from cdte_tools_py.io.fits_tools import load_fits

    >>> filename = f"path/to/file/cdte1_merged_v0_0_3.fits"
    >>> existing_cdte_hdus = load_fits(filename)

    >>> path1_cdte_hdus = get_path1(existing_cdte_hdus)

    >>> cdte_data = Table(path1_cdte_hdus[1].data)
    >>> filtered_cdte_data = pseudo_filter(cdte_data)
    """
    valid_trig = (cdte_data["flag_pseudo"]==1)
    return cdte_data[valid_trig]

def _single_hit_no_guard_ring_filter(cdte_data):
    """Filter out CdTe triggers that average strip positions are in the guard ring. 

    Will only check the first merged strip position in the 
    "XX_merged_position_list" lists (XX means either "pt" or "al").

    Filter is valid for CdTe data that has already had its ADC values 
    converted to energy (e.g., `path2` data from
    cdte_tools_py.calibration.path2.get_path2).
    
    Helps filter mimics the filter performed for `path2` in the C++ code
    in `path2_cal_and_merge_eachside.cpp`.

    Parameters
    ----------
    cdte_data : `astropy.table.Table`
        The Astropy table of the CdTe data with energy information.

    Returns
    -------
    `astropy.table.Table` : 
        Filters out averaged strips in the guard ring of the detector.
        - Edge strips of values <4 and >122.

    Examples
    --------
    >>> from astropy.table import Table
    >>> from cdte_tools_py.calibration.path1 import get_path1
    >>> from cdte_tools_py.calibration.path2 import get_path2
    >>> from cdte_tools_py.io.fits_tools import load_fits

    >>> filename = f"path/to/file/cdte1_merged_v0_0_3.fits"
    >>> existing_cdte_hdus = load_fits(filename)

    >>> path1_cdte_hdus = get_path1(existing_cdte_hdus)

    >>> path2_cdte_hdus = get_path2(path1_cdte_hdus)

    >>> cdte_data = Table(path2_cdte_hdus[1].data)
    >>> filtered_cdte_data = _single_hit_no_guard_ring_filter(cdte_data)
    """
    bool_pos = (cdte_data["pt_merged_position_list"][:,0]>3)\
                &(cdte_data["pt_merged_position_list"][:,0]<123)\
                &(cdte_data["al_merged_position_list"][:,0]>3)\
                &(cdte_data["al_merged_position_list"][:,0]<123)
    return cdte_data[bool_pos]

def single_pt_filter(cdte_data):
    """Leave only the CdTe triggers with only a single strip hit on Pt-side. 

    This also requires that the Al-side has at least one hit strip.

    Filter is valid for CdTe data that has already had its ADC values 
    converted to energy (e.g., `path2` data from
    cdte_tools_py.calibration.path2.get_path2).
    
    Helps filter mimics the filter performed for `path2` in the C++ code
    in `path2_cal_and_merge_eachside.cpp`.

    Parameters
    ----------
    cdte_data : `astropy.table.Table`
        The Astropy table of the CdTe data with energy information.

    Returns
    -------
    `astropy.table.Table` : 
        Filters for "pt_nhit"==1 and "al_nhit">0.

    Examples
    --------
    >>> from astropy.table import Table
    >>> from cdte_tools_py.calibration.path1 import get_path1
    >>> from cdte_tools_py.calibration.path2 import get_path2
    >>> from cdte_tools_py.io.fits_tools import load_fits

    >>> filename = f"path/to/file/cdte1_merged_v0_0_3.fits"
    >>> existing_cdte_hdus = load_fits(filename)

    >>> path1_cdte_hdus = get_path1(existing_cdte_hdus)

    >>> path2_cdte_hdus = get_path2(path1_cdte_hdus)

    >>> cdte_data = Table(path2_cdte_hdus[1].data)
    >>> filtered_cdte_data = single_pt_filter(cdte_data)
    """
    bool_pt_hit = (cdte_data["pt_nhit"]==1)&(cdte_data["al_nhit"]>0)
    return cdte_data[bool_pt_hit]

def single_al_filter(cdte_data):
    """Leave only the CdTe triggers with only a single strip hit on Al-side. 

    This also requires that the Pt-side has at least one hit strip.

    Filter is valid for CdTe data that has already had its ADC values 
    converted to energy (e.g., `path2` data from
    cdte_tools_py.calibration.path2.get_path2).
    
    Helps filter mimics the filter performed for `path2` in the C++ code
    in `path2_cal_and_merge_eachside.cpp`.

    Parameters
    ----------
    cdte_data : `astropy.table.Table`
        The Astropy table of the CdTe data with energy information.

    Returns
    -------
    `astropy.table.Table` : 
        Filters for "al_nhit"==1 and "pt_nhit">0.

    Examples
    --------
    >>> from astropy.table import Table
    >>> from cdte_tools_py.calibration.path1 import get_path1
    >>> from cdte_tools_py.calibration.path2 import get_path2
    >>> from cdte_tools_py.io.fits_tools import load_fits

    >>> filename = f"path/to/file/cdte1_merged_v0_0_3.fits"
    >>> existing_cdte_hdus = load_fits(filename)

    >>> path1_cdte_hdus = get_path1(existing_cdte_hdus)

    >>> path2_cdte_hdus = get_path2(path1_cdte_hdus)

    >>> cdte_data = Table(path2_cdte_hdus[1].data)
    >>> filtered_cdte_data = single_al_filter(cdte_data)
    """
    bool_al_hit = (cdte_data["pt_nhit"]>0)&(cdte_data["al_nhit"]==1)
    return cdte_data[bool_al_hit]

def ti_range_filter(cdte_data, lower_ti, upper_ti):
    """Leave only the CdTe triggers between the lower and upper Ti. 

    Filter is valid for CdTe data that comes out of the parser, so the 
    earliest human readable data. Only requires the field "ti":
    - lower_ti<=ti<upper_ti

    Note: For later data products, a suitable Ti range was selected to 
    ensure the data includes the flight and not too much time, e.g., on 
    the rail. If one of the original raw files (without this suitable 
    time range filtering) is used then be aware that Ti rolls over in 
    value periodically.

    Parameters
    ----------
    cdte_data : `astropy.table.Table`
        The Astropy table of the CdTe data with energy information.

    lower_ti, upper_ti : `int` or `float`, `int` or `float`
        The lower bound of the Ti time range (inclusive) and the upper
        bound Ti of the time range (exclusive), respectively. 
            The field `lower_ti`/`upper_ti` is inclusive/exclusive in 
        order to allow multiple ranges being selected to be stacked 
        without unwanted overlap. E.g., `ti_range_filter(data, 1, 2)` and
        `ti_range_filter(data, 2, 3)` represent two time range selections 
        where no double counting of triggers can occur.

    Returns
    -------
    `astropy.table.Table` : 
        Filters for triggers lower_ti<=ti<upper_ti.

    Examples
    --------
    >>> from astropy.table import Table
    >>> from cdte_tools_py.calibration.path1 import get_path1
    >>> from cdte_tools_py.calibration.path2 import get_path2
    >>> from cdte_tools_py.io.fits_tools import load_fits

    >>> filename = f"path/to/file/cdte1_merged_v0_0_3.fits"
    >>> existing_cdte_hdus = load_fits(filename)

    >>> path1_cdte_hdus = get_path1(existing_cdte_hdus)

    >>> cdte_data = Table(path1_cdte_hdus[1].data)
    >>> filtered_cdte_data = ti_range_filter(cdte_data, 
                                             path1_cdte_hdus[0].header["TI1"], 
                                             path1_cdte_hdus[0].header["TI7"])
    """
    bool_ti = (lower_ti<=cdte_data["ti"])&(cdte_data["ti"]<upper_ti)
    return cdte_data[bool_ti]

def energy_range_filter(cdte_data, lower_energy, upper_energy, _pre_gap_loss=False):
    """Leave only the CdTe triggers between the lower and upper Energy. 

    Filter is valid for CdTe data post-level 2 (after calibration step
    path 2 has been applied to the data). 
    
    If `_pre_gap_loss` is False, requires fields "pt_gap_loss_corrected_merged_energy_list" 
    and "al_gap_loss_corrected_merged_energy_list":
    - lower_energy<=pt_gap_loss_corrected_merged_energy_list[:,0]<upper_energy and
      lower_energy<=al_gap_loss_corrected_merged_energy_list[:,0]<upper_energy

    If `_pre_gap_loss` is True, requires fields "pt_merged_energy_list" 
    and "al_merged_energy_list":
    - lower_energy<=pt_merged_energy_list[:,0]<upper_energy and
      lower_energy<=al_merged_energy_list[:,0]<upper_energy

    Note: This filter stipulates that both the Pt-side and Al-side energy
    for a trigger must be within the energy range. The filter also assumes 
    that the triggers being filtered only consist of one strip clump on 
    both detector sides.

    Parameters
    ----------
    cdte_data : `astropy.table.Table`
        The Astropy table of the CdTe data with energy information.

    lower_energy, upper_energy : `int` or `float`, `int` or `float`
        In keV, the lower bound of the energy range (inclusive) and the 
        upper bound of the energy range (exclusive), respectively. 
            The field `lower_energy`/`upper_energy` is inclusive/exclusive 
        in order to allow multiple ranges being selected to be stacked 
        without unwanted overlap. E.g., 
        `energy_range_pre_gap_loss_filter(data, 1, 2)` and
        `energy_range_pre_gap_loss_filter(data, 2, 3)` represent two 
        energy range selections where no double counting of triggers can 
        occur.

    _pre_gap_loss : bool
        If True, then applies the energy filter to the energy fields 
        created before the gap-loss crrected energies were calculated.
        Applies the filter to fields "pt_merged_energy_list" and 
        "al_merged_energy_list". If False (default), will apply the 
        filter to the fields "pt_gap_loss_corrected_merged_energy_list"
        and "al_gap_loss_corrected_merged_energy_list".
        Default: False

    Returns
    -------
    `astropy.table.Table` : 
        Filters for triggers lower_energy<=(Pt & Al energy)<upper_energy.

    Examples
    --------
    (1) Filter level 2 gap-loss corrected energies to 5--6.9 keV
    >>> from cdte_tools_py.io.fits_tools import load_fits

    >>> cdte_level2 = load_fits("cdteN_level2_versionX.fits")

    >>> cdte_data = Table(cdte_level2)
    >>> filtered_cdte_data = energy_range_pre_gap_loss_filter(cdte_data, 
                                                              5, 
                                                              6.9)

    (2) Filter level 2 pre-gap-loss corrected energies to 5--6.9 keV
    >>> from cdte_tools_py.io.fits_tools import load_fits

    >>> cdte_level2 = load_fits("cdteN_level2_versionX.fits")

    >>> cdte_data = Table(cdte_level2)
    >>> filtered_cdte_data = energy_range_pre_gap_loss_filter(cdte_data, 
                                                              5, 
                                                              6.9,
                                                              _pre_gap_loss=True)

    (3) Filter Path 2 data to energies between 5 and 6.9 keV
    >>> from astropy.table import Table
    >>> from cdte_tools_py.calibration.path1 import get_path1
    >>> from cdte_tools_py.calibration.path2 import get_path2
    >>> from cdte_tools_py.io.fits_tools import load_fits

    >>> filename = f"path/to/file/cdte1_merged_v0_0_3.fits"
    >>> existing_cdte_hdus = load_fits(filename)

    >>> path1_cdte_hdus = get_path1(existing_cdte_hdus)
    >>> path2_cdte_hdus = get_path2(path1_cdte_hdus)

    >>> cdte_data = Table(path2_cdte_hdus[1].data)
    >>> filtered_cdte_data = energy_range_pre_gap_loss_filter(cdte_data, 
                                                              5, 
                                                              6.9,
                                                              _pre_gap_loss=True)
    """
    pt_energy_field = "pt_gap_loss_corrected_merged_energy_list"
    al_energy_field = "al_gap_loss_corrected_merged_energy_list"

    if _pre_gap_loss:
        pt_energy_field = "pt_merged_energy_list"
        al_energy_field = "al_merged_energy_list"

    bool_e = (lower_energy<=cdte_data[pt_energy_field][:,0])\
                &(cdte_data[pt_energy_field][:,0]<upper_energy)\
                    &(lower_energy<=cdte_data[al_energy_field][:,0])\
                        &(cdte_data[al_energy_field][:,0]<upper_energy)
    return cdte_data[bool_e]


# Converters
# ----------

# Extractors
# ----------
def strip_energy_clumps_extractor(cdte_data):
    """Return lists of the strip and energy clumpe per trigger.
    
    This is without dummy values.

    Helps filter mimics the filter output for `path2` in the C++ code
    in `path2_cal_and_merge_eachside.cpp`.

    Parameters
    ----------
    cdte_data : `astropy.table.Table`
        The Astropy table of the CdTe data with energy information.

    Returns
    -------
    `dict[list[list[numpy.array]]]` : 
        Two lists where the entries each contain a list of Numpy arrays 
        which are the strip/energy clumps and correspond to each trigger.
        Lists contain the strip and energy clumps for both detector sides 
        (keys: "pt_strip_clumps", "pt_energy_clumps", "pt_strip_clumps", 
        "pt_energy_clumps").

    Examples
    --------
    >>> from astropy.table import Table
    >>> from cdte_tools_py.calibration.path1 import get_path1
    >>> from cdte_tools_py.calibration.path2 import get_path2
    >>> from cdte_tools_py.io.fits_tools import load_fits

    >>> filename = f"path/to/file/cdte1_merged_v0_0_3.fits"
    >>> existing_cdte_hdus = load_fits(filename)

    >>> path1_cdte_hdus = get_path1(existing_cdte_hdus)

    >>> path2_cdte_hdus = get_path2(path1_cdte_hdus)

    >>> cdte_data = Table(path2_cdte_hdus[1].data)
    >>> clump_cdte_data = strip_energy_clumps_extractor(cdte_data)

    >>> print(clump_cdte_data)
    {'pt_strip_clumps': [[array([45, 46, 47, 48, 49, 50, 51, 52]),
                          array([112, 113, 114, 115, 116, 117, 118, 119, 120, 121])],
                        [array([], dtype=int64)],
                        [array([69, 70, 71])], 
                        ...],
     'pt_energy_clumps': [[array([2.16712927, 1.97407416, 6.86221546, 5.13268624, 6.30540948,
                                  4.36840715, 1.72018145, 3.6211327 ]),
                           array([2.34532438, 2.50425944, 2.22570308, 2.10514853, 3.97429507,
                                  2.27513143, 4.17305479, 5.11113677, 3.64069373, 4.328372  ])],
                         [array([], dtype=float64)],
                         [array([4.06837925, 2.3098551 , 5.07257457])],
                         ...],
     'al_strip_clumps': [[array([], dtype=int64)],
                         [array([ 8,  9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24,
                                 25, 26, 27, 28])],
                         [array([109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121,
                                 122])],
                         ...],
     'al_energy_clumps': [[array([], dtype=float64)],
                          [array([5.89762788, 5.11891298, 6.06494816, 4.69267884, 5.82365456,
                                  4.90829596, 5.10940535, 4.48578127, 4.42433919, 4.48548716,
                                  3.75547534, 4.45320757, 4.51964432, 3.92108359, 2.49585796,
                                  2.71299407, 1.72115739, 1.68886558, 2.45447938, 2.66269441,
                                  1.66666141])],
                          [array([2.68924521, 2.07152897, 2.73367524, 2.69262692, 2.85876764,
                                  2.36293439, 2.97037131, 3.44922074, 4.25779235, 3.18408042,
                                  2.9769731 , 3.64827478, 3.44199856, 3.60150054])],
                         ...]}
    """
    all_pt_strip_clumps, all_pt_energy_clumps = [], []
    all_al_strip_clumps, all_al_energy_clumps = [], []
    # loop through all triggers
    for tigger_pt_energies, tigger_al_energies in zip(cdte_data["pt_valid_energy_array"], 
                                                      cdte_data["al_valid_energy_array"]):
        # get the non-dummy values
        valid_pt_energy_inds = np.nonzero(tigger_pt_energies!=UNPHYSICAL_ENERGY)[0]
        valid_al_energy_inds = np.nonzero(tigger_al_energies!=UNPHYSICAL_ENERGY)[0]

        # split the consecutive indices into wither own lists 
        # (these are the strip clumps for the trigger)
        stepsize = 1
        trigger_pt_strip_clumps = np.split(valid_pt_energy_inds, 
                                        np.nonzero(np.diff(valid_pt_energy_inds)!=stepsize)[0]+1)
        trigger_al_strip_clumps = np.split(valid_al_energy_inds, 
                                        np.nonzero(np.diff(valid_al_energy_inds)!=stepsize)[0]+1)

        # loop through the strip clumps and get the energy clumps
        trigger_pt_energy_clumps = [tigger_pt_energies[inds] for inds in trigger_pt_strip_clumps]
        trigger_al_energy_clumps = [tigger_al_energies[inds] for inds in trigger_al_strip_clumps]

        # save the strip/energy clumps for the trigger and move onto the next
        all_pt_strip_clumps.append(trigger_pt_strip_clumps)
        all_pt_energy_clumps.append(trigger_pt_energy_clumps)
        all_al_strip_clumps.append(trigger_al_strip_clumps)
        all_al_energy_clumps.append(trigger_al_energy_clumps)
        
    return {"pt_strip_clumps":all_pt_strip_clumps, 
            "pt_energy_clumps":all_pt_energy_clumps, 
            "al_strip_clumps":all_al_strip_clumps, 
            "al_energy_clumps":all_al_energy_clumps}