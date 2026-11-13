"""
Module containing methods to load and save fits files
"""

import numpy as np
from astropy.io import fits
from astropy.table import Table, vstack
from astropy.time import Time

from foxsi4_science_tools_py.io.fits_tools import load_fits

def load_cdte(f):
    hdul = load_fits(f)
    return [hdul[0].header, Table(hdul[1].data)]

def load_cmos(f):
    hdul = load_fits(f)
    return [hdul[0].header, Table(hdul[1].data)]

def load_cmos_old(f):
    # pre-energy calibration
    return [cmos_header(f), cmos_construct(f)]

def cmos_iter_event_chunks(path, chunk_rows):
    with fits.open(path, memmap=True) as hdul:
        events = hdul["EVENTS"].data
        for start in range(0, len(events), chunk_rows):
            yield events[start : min(len(events), start + chunk_rows)]

def cmos_header(path):
    with fits.open(path, memmap=True) as hdul:
        h = hdul[0].header
    return h

def cmos_construct(path, chunk_rows=500_000):
    main_table = Table()
    for rows in cmos_iter_event_chunks(path, chunk_rows):
        # do grade selection
        datatype = np.asarray(rows["DATATYPE"], dtype=np.int16)
        group_ev = np.asarray(rows["GROUP_EV"], dtype=np.int32)
        edge_flg = np.asarray(rows["EDGE_FLG"], dtype=np.uint8)
        selection = (datatype == 1) & (group_ev == 1) & (edge_flg == 0)

        tmp_table = Table()
        # time
        tmp_table["utc"] = Time(np.asarray([utc.replace("+00:00", "") for utc in rows["UTC"]]))[selection]
        # space
        tmp_table["solar_x"] = np.asarray(rows["SOLAR_X"], dtype=np.float64)[selection]
        tmp_table["solar_y"] = np.asarray(rows["SOLAR_Y"], dtype=np.float64)[selection]
        # energy
        tmp_table["dn"] = np.asarray(rows["CORR_CNT"], dtype=np.float64)[selection]

        main_table = vstack([main_table, tmp_table])
    return main_table
