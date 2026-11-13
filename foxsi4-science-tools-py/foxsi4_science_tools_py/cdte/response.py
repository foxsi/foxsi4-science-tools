"""Helper unctions to get a CdTe detector's response."""

import re
import warnings

import numpy as np
import response_tools.responses as responses

def pitch_fractions(event_object):
    """Get the fraction of count detected by each pitch region."""
    pitches = event_object.extra_count_info["pitch"]
    total_counts = len(pitches)
    _0 = pitches[np.nonzero(pitches==0)]
    _1 = pitches[np.nonzero(pitches==1)]
    _2 = pitches[np.nonzero(pitches==2)]
    pitch0_fraction = np.sum(_0[0])/total_counts if len(_0>0) else 0
    pitch1_fraction = np.sum(_1[0])/total_counts if len(_1>0) else 0
    pitch2_fraction = np.sum(_2[0])/total_counts if len(_2>0) else 0
    _gr = pitches[np.nonzero(pitches==-1)]
    _np = pitches[np.nonzero(pitches==-2)]
    guard_ring = np.sum(_gr[0])/total_counts if len(_gr>0) else 0
    nonphysical = np.sum(_np[0])/total_counts if len(_np>0) else 0
    return {"pitch0_fraction":pitch0_fraction,
            "pitch1_fraction":pitch1_fraction,
            "pitch2_fraction":pitch2_fraction,
            "guard_ring":guard_ring,
            "nonphysical":nonphysical,
            }

def download_response_files():
    import response_tools.io.fetch_response_data as fetch
    fetch.foxsi4_download_required(verbose=True)

def construct_cdte_response(event_object, off_axis_angle=None, time_range=None, det_to_tel_mapping=None):
    # first, we need the CdTe detector number to do any of this
    
    det_id_candidates = re.findall(r'cdte\d+', event_object.name.lower())
    if (len(det_id_candidates)!=1):
        warnings.warn("Unable to understand `event_object.name` ({event_object.name}) attribute.")
        return
    cdte_id_num = int(det_id_candidates[0][-1])
    if cdte_id_num not in list(det_to_tel_mapping.keys()):
        warnings.warn(f"The `event_object.name` ({event_object.name}) attribute\'s number ID ({cdte_id_num}) is not in det_to_tel_mapping ({det_to_tel_mapping}).")
        return
    
    telescope_number = det_to_tel_mapping[cdte_id_num]
    try:
        rmf = rmf_construct(event_object.misc_info["signal_source"], 
                            telescope_rmf_func(telescope_number), 
                            pitch_fractions(event_object))
    except FileNotFoundError:
        download_response_files()
        rmf = rmf_construct(event_object.misc_info["signal_source"], 
                            telescope_rmf_func(telescope_number), 
                            pitch_fractions(event_object))

    arf = arf_construct(telescope_arf_func(telescope_number), rmf, off_axis_angle=off_axis_angle, time_range=time_range)

    return responses.foxsi4_telescope_spectral_response(arf, rmf)

def arf_construct(arf_func, rmf_obj, off_axis_angle=None, time_range=None):
    rmf_mid_energies = (rmf_obj.input_energy_edges[:-1]+rmf_obj.input_energy_edges[1:])/2
    return arf_func(mid_energies=rmf_mid_energies,
                            off_axis_angle=off_axis_angle,
                            time_range=time_range)

def telescope_rmf_func(telescope_number):
    """Return the correct telescope RMF function."""
    rmf_func = {2:responses.foxsi4_telescope2_rmf, 
                3:responses.foxsi4_telescope3_rmf, 
                4:responses.foxsi4_telescope4_rmf, 
                5:responses.foxsi4_telescope5_rmf,
                }
    return rmf_func[telescope_number]

def _rmf_breakdown(signal_source, rmf_func):
    """Get the three RMFs for a CdTes three pitches."""
    if signal_source=="Pt-side":
        _signal_source = "pt"
    elif signal_source=="DOI":
        _signal_source = "merged"
    else:
        warnings.warn("The `signal_source` input should be \'Pt-side\' or \'DOI\', not \'{signal_source}\'.")
        return
    pos_rmf0 = rmf_func(region=0, _side=_signal_source)
    pos_rmf1 = rmf_func(region=1, _side=_signal_source)
    pos_rmf2 = rmf_func(region=2, _side=_signal_source)
    return {"pitch0_rmf":pos_rmf0,
            "pitch1_rmf":pos_rmf1,
            "pitch2_rmf":pos_rmf2,
            }

def _rmf_construct(rmf_objs, pitch_fracs):
    """Given the count pitch fractions, construct the CdTe RMF array."""
    pos_rmf0_contribution = rmf_objs["pitch0_rmf"].response * pitch_fracs["pitch0_fraction"]
    pos_rmf1_contribution = rmf_objs["pitch1_rmf"].response * pitch_fracs["pitch1_fraction"]
    pos_rmf2_contribution = rmf_objs["pitch2_rmf"].response * pitch_fracs["pitch2_fraction"]
    return pos_rmf0_contribution + pos_rmf1_contribution + pos_rmf2_contribution

def rmf_construct(signal_source, rmf_func, pitch_fracs):
    """Get the combined RMF object made for the combinations of pitches."""
    rmf_objs = _rmf_breakdown(signal_source, rmf_func)
    rmf_array = _rmf_construct(rmf_objs, pitch_fracs)
    return responses.Response2DOutput(filename="No-File", 
                                      function_path="No-Path",
                                      input_energy_edges=rmf_objs["pitch0_rmf"].input_energy_edges,
                                      output_energy_edges=rmf_objs["pitch0_rmf"].output_energy_edges,
                                      response=rmf_array,
                                      response_type="Combined-RFM",
                                      telescope=f"{rmf_objs["pitch0_rmf"].telescope}",
                                      elements=(rmf_objs["pitch0_rmf"],
                                              rmf_objs["pitch1_rmf"],
                                              rmf_objs["pitch2_rmf"],
                                              ),
                                     )

def telescope_arf_func(telescope_number):
    """Return the correct telescope ARF function."""
    arf_func = {2:responses.foxsi4_telescope2_flight_arf, 
                3:responses.foxsi4_telescope3_flight_arf, 
                4:responses.foxsi4_telescope4_flight_arf, 
                5:responses.foxsi4_telescope5_flight_arf}
    return arf_func[telescope_number]
