""" Module with utility functions to perform typical cleaning of event list."""

from dataclasses import dataclass, field, fields

import astropy.units as u
from astropy.time import Time
import numpy as np

from foxsi4_science_tools_py.cdte.strip_positions import get_hit_pitch

@dataclass
class UsualSuspects:
    name: str
    utc: Time
    solar_x: u.Quantity 
    solar_y: u.Quantity
    signal: u.Quantity 
    extra_count_info: dict = field(default_factory=dict)
    misc_info: dict = field(default_factory=dict)

    def __getitem__(self, key):
        if isinstance(key, int):
            # Access by integer index (e.g., pt[0])
            field_name = fields(self)[key].name
            return getattr(self, field_name)
        elif isinstance(key, str):
            # Access by field name string (e.g., pt['x'])
            return getattr(self, key)
        elif isinstance(key, np.ndarray):
            # this must be indexing
            utc = self.utc[key]
            solar_x = self.solar_x[key]
            solar_y = self.solar_y[key]
            signal = self.signal[key]
            _extra_count_info = {}
            if len(self.extra_count_info)>0:
                for k,v in self.extra_count_info.items():
                    _extra_count_info[k] = v[key]
            return UsualSuspects(name=self.name, utc=utc, solar_x=solar_x, solar_y=solar_y, signal=signal, extra_count_info=_extra_count_info, misc_info=self.misc_info)
        raise TypeError(f"Invalid subscript type: {type(key)}")

def clean_cdte(event_list):
    return event_list[event_list["flag_quality"]==1]

def _cdte_usual_suspects(event_list, name="cdte", signal_source="Pt-side"):
    time = Time(event_list["utc"])
    solx = event_list["solar_x_unaligned"] << u.arcsec
    soly = event_list["solar_y_unaligned"] << u.arcsec
    if signal_source=="DOI":
        signal = event_list["doi_corrected_energy"]
    elif signal_source=="Pt-side":
        signal = event_list["pt_gap_loss_corrected_merged_energy_list"][:, 0]

    # extra CdTe information needed for analysis
    _pt_strips = event_list["pt_merged_position_list"][:,0]
    _pitch = np.array([get_hit_pitch(strip) for strip in _pt_strips])
    _livetime = event_list["livetime"]*1e-8 << u.second
    extra_count_info = {"pitch":_pitch, "livetime":_livetime}

    misc_info = {"signal_source":signal_source}

    return UsualSuspects(name=name, utc=time, solar_x=solx, solar_y=soly, signal=signal<<u.keV, extra_count_info=extra_count_info, misc_info=misc_info)

def cdte_usual_suspects(event_list, name="cdte", signal_source="Pt-side"):
    return _cdte_usual_suspects(clean_cdte(event_list), name=name, signal_source=signal_source)

def clean_cmos(event_list):
    return event_list

def _cmos_usual_suspects(event_list, name="cmos"):
    # pre-energy calibration
    time = Time("2024-04-17T00:00:00") + (event_list["TIME"]<<u.second)
    solx = event_list["SOLAR_X"] << u.arcsec
    soly = event_list["SOLAR_Y"] << u.arcsec
    signal = event_list["ENERGY"] << u.keV
    return UsualSuspects(name=name, utc=time, solar_x=solx, solar_y=soly, signal=signal)

def _cmos_usual_suspects_old(event_list, name="cmos"):
    # pre-energy calibration
    time = Time(event_list["utc"])
    solx = event_list["solar_x"] << u.arcsec
    soly = event_list["solar_y"] << u.arcsec
    signal = event_list["dn"] << u.dn
    return UsualSuspects(name=name, utc=time, solar_x=solx, solar_y=soly, signal=signal)

def cmos_usual_suspects(event_list, name="cmos"):
    return _cmos_usual_suspects(clean_cmos(event_list), name=name)