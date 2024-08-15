"""File to host code that downloads SDO data."""

import astropy.units as u
from sunpy.net import Fido, attrs as a
import os

def sdo_download(start_time=None, end_time=None, directory="./", wave_values=None, get_hmi=False):
    """Allows the download of SDO/AIA and /HMI data.
    
    Parameters
    ----------
    start_time : `str`, `None` 
        The start time for the data. In the following format; e.g., 
        "2024-04-17T21:57:00" which is the value taken if input is 
        `None`.
        Default: None
    
    end_time : `str`, `None`
        The end time for the data. In the following format; e.g., 
        "2024-04-17T22:30:00" which is the value taken if input is 
        `None`.
        Default: None
    
    directory : `str`
        The directory where the SDO files will be saved. Should end in 
        a "/". Default is the current directory.
        Default: "./" 
    
    wave_values : `list[int]`, None
        A list of the wavelengths (in angstroms) to be downloaded. For 
        example, [94, 131, 171, 211, 1600, 1700] which is the value 
        taken if input is `None`.
        Default: None 
    
    get_hmi : `Bool`
        A boolean flag where, if True, HMI LOS data will also be 
        downloaded.
        Default: False
    
    Returns
    -------
    : `list[str]`
        A list object containing the names of the files that failed to
        download.

    Examples
    --------
    ## Default behaviour for the FOXSI M1 flare time
    >>> sdo_download()

    ## Extend the time interval, download a subset of wavelengths, 
    ## define a specific directory, and get the SDO/HMI data as well
    >>> sdo_download(start_time="2024-04-17T21:30:00", 
                     end_time="2024-04-17T22:45:00", 
                     directory="path/to/save/directory/", 
                     wave_values=[94,171], 
                     get_hmi=True)

    ## Only download SDO/HMI LOS data
    >>> sdo_download(wave_values=[], get_hmi=True)
    """
    # e.g., start_time, end_time = '2018-09-10T12:30:00', '2018-09-10T13:15:00'
    # wave_values in angstroms

    start_time = "2024-04-17T21:57:00" if start_time is None else start_time
    end_time = "2024-04-17T22:30:00" if end_time is None else end_time
    wave_values = [94, 131, 171, 193, 211, 304, 335, 1600, 1700] if wave_values is None else wave_values

    time = a.Time(start_time, end_time)
    aia = a.Instrument('AIA')
    waves = wave_values*u.angstrom 
    wave_dirs = [f"{wv}angstrom/" for wv in wave_values]

    needed_files = []

    # for aia
    for wave, wd in zip(waves, wave_dirs):
        print("Doing "+directory+wd)
        os.makedirs(directory+wd, exist_ok=True) # make the directory if it isn't there

        search = (time, aia, a.Wavelength(wave))
        fetch = {"path":directory+wd+"{file}"}

        filepaths = fido_download(search, fetch)

        # well try again with failed files
        needed_files = handle_retries(filepaths, needed_files=needed_files)

    # for hmi
    if get_hmi:
        print("Doing "+directory+"hmi")
        os.makedirs(directory+"hmi", exist_ok=True) # make the directory if it isn't there

        search = (time, a.Instrument('hmi'), a.Physobs('LOS_magnetic_field'))
        fetch = {"path":directory+"hmi/{file}"}

        filepaths = fido_download(search, fetch)

        # well try again with failed files
        needed_files = handle_retries(filepaths, needed_files=needed_files)

    print("Files Needed are:")
    print(needed_files)

    return needed_files

def fido_download(search, fetch=None):
    """Handle the common download elements for `sunpy.Fido`.
    
    Parameters
    ----------
    search : `tuple`
        The inputs given to the `Fido.search` query.
    
    fetch : `dict`
        Kyeword arguments to be passed to the `Fido.fetch` method. If 
        `None` then an empty `dict` object is passed.
        Default: None
    
    Returns
    -------
    : `sunpy.UnifiedResponse` object
        A list object containing the names of the files that failed to
        download.

    Examples
    --------
    ## Download some SDO/AIA data
    >>> from sunpy.net import attrs as a
    >>> search = (a.Time("2024-04-17T21:30:00", 
                         "2024-04-17T22:30:00"), 
                  a.Instrument('AIA'), 
                  a.Wavelength(wave))
    >>> fetch = {"path":"path/to/save/{file}"}
    >>> filepaths = fido_download(search, fetch)

    ## Download some SDO/HMI data
    >>> from sunpy.net import attrs as a
    >>> search = (a.Time("2024-04-17T21:30:00", 
                         "2024-04-17T22:30:00"), 
                  a.Instrument('hmi'), 
                  a.Physobs('LOS_magnetic_field'))
    >>> fetch = {"path":"path/to/save/{file}"}
    >>> filepaths = fido_download(search, fetch)
    """
    unifresp = Fido.search(*search)

    fetch = {} if fetch is None else fetch
    filepaths = Fido.fetch(unifresp, **fetch)

    return filepaths

def handle_retries(filepaths, tries=5, needed_files=None):
    """If a file(s) isn't downloaded straight away, retry the download.
    
    Parameters
    ----------
    filepaths : `sunpy.UnifiedResponse` object
        An object used to store the failed downloaded files.
    
    tries : `int`
        The number of times to retry downlading the given files.
        Default: 5
    
    needed_files : `list[str]`
        A list of files already needed to be returned if failed after 
        retrying `tries` times.For example, [] which is the value taken 
        if input is `None`. If this is a list object then the failed 
        files will be appended to this list.
        Default: None
    
    Returns
    -------
    : `list[str]`
        A list object containing the file names of the failed downloads
        after retrying.

    Examples
    --------
    ## Retry failed downloads and return a list still un-downloaded
    >>> needed_files = handle_retries(filepaths)

    ## If a list of un-downloaded alreadt exists then pass to append
    >>> needed_files = handle_retries(filepaths, 
                                      needed_files=needed_files)
    """

    needed_files = [] if needed_files is None else needed_files

    for _ in range(tries):
        print(f"Trying again for {filepaths.errors}")
        filepaths = Fido.fetch(filepaths)
        if filepaths.errors == []:
            break

    if filepaths.errors != []:
        needed_files.append(filepaths.errors)

    return needed_files

if __name__=="__main__":
    sdo_download(directory="./data/")
