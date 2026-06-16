"""
Module containing methods to load and save fits files
"""
from copy import deepcopy
import logging
import os

from astropy.io import fits
from astropy.table import Table, vstack
from astropy.utils.diff import report_diff_values

from foxsi4_science_tools_py.util.loops import progress_printer

def load_fits(filename):
    """
    Read a .fits file

    Parameters
    ----------
    filename : `str`, `file-like` or `pathlib.Path`
        A .fits file 

    Returns
    -------
    'list' :
        HDUList or PrimaryHDU of file.
    """

    with fits.open(filename) as hdul:
        hdus = deepcopy(hdul)

    return hdus

def save_fits(hdul, filename, **kwargs):
    """
    Write data to a .fits file

    Parameters
    ----------
    hdul : 'list', HDUList or PrimaryHDU
        A list of all hdus to be saved into the .fits file

    filename : `str`, `file-like` or `pathlib.Path`
        A .fits file
    
    **kwargs : stuff to be passed to fits.writeto(), such as overwrite

    Returns
    -------
    None
    """

    try:
        if not isinstance(hdul, fits.hdu.hdulist.HDUList):
            hdul = fits.HDUList(hdul)

        hdul.writeto(filename, **kwargs)

        logging.info(f"Succesfully saved data to {filename}")
    
    except Exception as err:
        logging.warning(f"Error while saving to fits: {err}")
        raise

def create_fits(data, metadata):
    """Create a FITS structure for the telemetry and DE merged data. 
    
    Parameters
    ----------
    data : `numpy.array`
        The data to be saved. Here should be an event list. Passed to
        `fits.BinTableHDU`.
    
    metadata : `dict`
        A dictionary compatible with FITS standard. Passed to 
        `fits.PrimaryHDU`.
    
    Returns
    -------
    : `list`, HDUList or PrimaryHDU
        A list of all HDUs to be saved into the .fits file

    Examples
    --------
    >>> from cdte_tools_py.merging.de_telemetry import merge_event_lists
    >>> from cdte_tools_py.merging.merging_metadata import merging_metadata
    >>> cdte_no = 1
    >>> telemetry_file = "/path/to/telemetry/file"
    >>> de_file = "/path/to/de/file"
    >>> data, meta = merge_event_lists(cdte_no, 
                                       telemetry_file, 
                                       de_file)
    >>> fits_struct = create_fits(data=data, 
                                  metadata=merging_metadata(cdte_no, 
                                                            telemetry_file, 
                                                            de_file, 
                                                            **meta))
    """
    header = fits.Header(metadata)

    primary_hdu = fits.PrimaryHDU(header=header)
    table_hdu = fits.BinTableHDU(data=data)

    hdu_list = fits.HDUList([primary_hdu, table_hdu])
    
    return hdu_list

def merge_hduls(hdul1, hdul2):
    """Function to merge two HDULs. 

    The resulting HDUL file will contain the header from `hdul1` and the
    merged event list table from both `hdul1` and `hdul2`.
    
    Parameters
    ----------
    hdul1, hdul2 : `str`
        The HDULs to be merged.
    
    Returns
    -------
    : `astropy.fits.HDUList`

    Examples
    --------
    """
    return fits.HDUList([hdul1[0], 
                         fits.BinTableHDU(data=vstack([Table(hdul1[1].data), 
                                                       Table(hdul2[1].data)
                                                       ])
                                          )
                         ])

def merge_fits(fits1, fits2, new_fits_name, clean_up=False, safely_clean_up=False, **kwargs):
    """Function to merge two FITS files. 

    The resulting FITS file will contain the header from `fits1` and the
    merged event list table from both `fits1` and `fits2`.
    
    Parameters
    ----------
    fits1, fits2 : `str`
        The FITS files to be merged.

    new_fits_name : `str`
        The name for the new FITS file. Must be different from `fits1` 
        and `fits2`.
    
    clean_up : `bool`
        If True, delete `fits1` and `fits2` when finished. For any issues, 
        this clean up input only relies on an error occurring during 
        loading/saving and checking if a file by the name `new_fits_name` 
        exists for the files to be deleted.
            Since this is destructive, please consider using 
        `safely_clean_up`; however, it will be more time and memory 
        intensive.
    
    safely_clean_up : `bool`
        If True, delete `fits1` and `fits2` while performing the same 
        checks as `clean_up=True` but also load in the created file and 
        check the contents.
            This option will be more time and memory intensive than using
        `clean_up=True`.

    **kwargs : Passed to `cdte_tools_py.io.fits_tools.save_fits`.
    
    Returns
    -------
    : `None`

    Examples
    --------
    """
    if (fits1==new_fits_name) or (fits2==new_fits_name):
        raise ValueError(f"Either fits1={fits1} or fits2={fits2} is the same as new_fits_name={new_fits_name}.")

    hdul1, hdul2 = load_fits(fits1), load_fits(fits2)

    new_hdul = merge_hduls(hdul1, hdul2)

    save_fits(new_hdul, new_fits_name, **kwargs)

    # make sure a file with the new name exists
    if not os.path.isfile(new_fits_name):
        raise AssertionError(f"File {new_fits_name} was not created in merging process.")
    
    if safely_clean_up:
        # check the new file saved properly (True if no diff, else False).
        new_contents = load_fits(new_fits_name)
        if not report_diff_values(Table(new_hdul[1].data), Table(new_contents[1].data)):
            raise AssertionError(f"File {new_fits_name} contents not as expected.")
        del new_contents
        clean_up = True
    
    # delete the old files if needed but check the new file at least exists
    if clean_up:
        os.remove(fits1)
        os.remove(fits2)

    # avoid waiting for the system to clean up
    del hdul1, hdul2, new_hdul

def merge_multifits(fits_list, new_fits_name, progress=False, **kwargs):
    """Function to merge two FITS files. 

    The resulting FITS file will contain the header from the first FITS
    file and the merged event list table from all FITS files.

    **Note:** Consider setting `clean_up` or `safe_clean_up` to `True` 
    if you want to delete the FITS files in `fits_list` as they're being 
    merged.

    **Note:** Files will be merged 2 at a time to avoid opening loads of 
    files at once. The new intermediate merged file will then be merged 
    wiith the next in `fits_list`. Intermediate merged files will be 
    deleted after they've been merged with the next file in `fits_list`.
    
    Parameters
    ----------
    fits_list : `list[str]`
        The list of FITS files to be merged.

    new_fits_name : `str`
        The name for the new final FITS file.
    
    progress : `bool`
        Print a statement indicating the progress while looping the file 
        list.
        Default: False

    **kwargs : Passed to `cdte_tools_py.io.fits_tools.merge_fits`.
        Consider setting `clean_up` or `safe_clean_up` to `True`.
    
    Returns
    -------
    : `None`

    Examples
    --------
    """
    fits_list_len = len(fits_list)
    merge_fits(fits_list[0], fits_list[1], f"{new_fits_name}{1}", **kwargs)
    for i in range(2, fits_list_len):
        if progress:
            progress_printer("Merging FITS files:", i+1, fits_list_len)

        prev_file = f"{new_fits_name}{i-1}"
        merge_fits(prev_file, fits_list[i], f"{new_fits_name}{i}", **kwargs)

        # make sure to clean up intermediate merged files
        if os.path.isfile(prev_file):
            os.remove(prev_file)

    os.rename(f"{new_fits_name}{i}", new_fits_name)