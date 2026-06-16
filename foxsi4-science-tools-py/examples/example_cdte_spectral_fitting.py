"""
Looking to fit CdTe data
========================

Kris' example on using the FOXSI code suite and Sunkit-spex to perform
spectral fitting on some FOXSI-4 CdTe data.

The file has a lot of lines but far fewer than 100 lines are actually 
doing anything, the rest are all plotting code or comments.

Directly below, I:
1. Import all the packages I want to use.
2. Define some global variables so they are easily changed:
   - CDTE_NUM is the CdTe number you want to look at (1-4)
   - CDTE_FILE is the path and file you want loaded in
      - I have it so that CDTE_NUM is used it it as well
   - ENERGIES, the type of energies you want to fit 
      - "Pt-side" (recommended) or "DOI"
   - SAVE_DIR=string if you want figures to be saved, else None
   - VARY_FE (bool), if you want to vary the Fe abundance during fitting
   - VARY_GAIN (bool), if you want to vary the gain during fitting
3. Define a custom thermal model for the fitting that allows me to vary 
   only the Fe abundance if I want (Feldman 1992 coronal abundances used 
   by default).

Below, in the `if __name__=="__main__:" section, I go through:
1. filtering the CdTe event list in time, then selecting a region, 
   before inspecting the count spectrum.
2. I then take the selections made to the data and use them to produce 
   the responses I need (time range, off-axis angle, strips used, etc.).
3. Pass the filtered data and the selected response to Sunkit-spex where
   I fit an isothermal model (with varying Fe abundance).

You will need the CdTe level files you want to load in.

Packages you will need:
- Astropy, matplotlib, numpy
- cdte_tools_py [https://github.com/foxsi/cdte-tools]
- response_tools [https://github.com/foxsi/response-tools]
- sunkit_spex [https://github.com/sunpy/sunkit-spex]
"""
#%% 
# Imports

# general packages
from astropy.table import Table
from astropy.time import Time
import astropy.units as u
from astropy.visualization import time_support
from matplotlib.dates import DateFormatter
from matplotlib.colors import LogNorm
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import numpy as np

# FOXSI code suite
import cdte_tools_py
import cdte_tools_py.calibration.filter_convert_extract as fce
from cdte_tools_py.calibration.path3 import get_hit_region
from cdte_tools_py.io.fits_tools import load_fits
from cdte_tools_py.solar_coords.solar_coords import strip_ref_arcsec
import response_tools.responses as responses

# imports for the spectral fitting
import warnings
from numpy.exceptions import VisibleDeprecationWarning
from sunkit_spex.legacy.fitting.fitter import Fitter
from sunkit_spex.models.physical.thermal import thermal_emission
warnings.filterwarnings("ignore", category=RuntimeWarning)
try:
    warnings.filterwarnings("ignore", category=VisibleDeprecationWarning)
except AttributeError:
    warnings.filterwarnings("ignore", category=np.exceptions.VisibleDeprecationWarning)

#%% 
# Some global variables
CDTE_NUM = 1
CDTE_FILE = f"/Users/kris/Documents/umnPostdoc/projects/both/foxsi4/foxsi4-analysis/data/cdte-levels/cdte{CDTE_NUM}_level3_version2.fits"
ENERGIES = "Pt-side" # or "DOI"

SAVE_DIR = None #"/Users/kris/Documents/umnPostdoc/projects/both/foxsi4/foxsi4-analysis/figs/"
VARY_FE = False
VARY_GAIN = False

#%% 
# Function for the fitting later

def f_fec_vth(temp, em46, fea, energies=None):
    """Calculates optically thin thermal bremsstrahlung radiation as seen from Earth.

    [1] https://hesperia.gsfc.nasa.gov/ssw/packages/xray/idl/f_vth.pro

    Parameters
    ----------
    energies : 2d array
            Array of energy bins for the model to be calculated over.
            E.g., [[1,1.5],[1.5,2],[2,2.5],...].

    temperature : int or float
            Plasma temperature in megakelvin.

    emission_measure46 : int or float
            Emission measure in units of 1e46 cm^-3.

    Returns
    -------
    A 1d array of optically thin thermal bremsstrahlung radiation in units
    of ph s^-1 cm^-2 keV^-1.
    """
    from sunkit_spex.models.physical.thermal import thermal_emission
    import astropy.units as u

    # turn [[1,2],[2,3],[3,4]] into [1,2,3,4]
    energies = np.unique(np.array(energies).flatten()) << u.keV
    temperature = temp * 1e6 << u.K
    emission_measure = em46 * 1e46 << u.cm ** (-3)
    return thermal_emission(energies, temperature, emission_measure, relative_abundances=((26, fea),)).value

# the main code
if __name__=="__main__":
    # %%
    # CdTe data
    # ---------
    #
    # Start by looking at the flight data.
    #
    # Load in the CdTe data and filter
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #
    # let"s look at the data for FOXSI-4 telescope 2
    # - Thermal blanket -> Marshall 10-shell X-7 -> Al (0.015”) -> CdTe4

    # point to level CdTe N file
    cdte_hdul = load_fits(CDTE_FILE)
    cdte_data = Table(cdte_hdul[1].data)

    # only select nominal quality triggers
    cdte_data_ge = cdte_data[cdte_data["flag_quality"]==1]

    # %%
    # Observational selections
    # ~~~~~~~~~~~~~~~~~~~~~~~~
    #
    # Start making selections based on the features I want to analyse.
    #
    # Temporal
    # ~~~~~~~~
    #
    # I only want the data from the first flare pointing just now:

    # time of interest for the flight, let"s go for the first look at the flare
    ti_obs_start, ti_obs_end = cdte_hdul[0].header["TI2"], cdte_hdul[0].header["TI3"]
    cdte_data_toi = fce.ti_range_filter(cdte_data_ge, ti_obs_start, ti_obs_end)
    obs_time = (ti_obs_end - ti_obs_start) * 160e-9 << u.second

    ref_time = cdte_tools_py.contextCdTeInfo["launch_time"]["launch_utc"]

    binSize = 5 # seconds

    times_all = Time(cdte_data_ge["utc"], format="isot",scale="utc")
    binRange_all = [min(times_all.unix_tai),max(times_all.unix_tai)]
    n_bins_all = int((binRange_all[1] - binRange_all[0]) / binSize)
    bins_all = np.arange(n_bins_all)*binSize+binRange_all[0]
    hist_all, _bin_edges_all = np.histogram(times_all.unix_tai, bins=bins_all)
    bin_edges_all = Time(_bin_edges_all, format="unix_tai",scale="utc").datetime

    times = Time(cdte_data_toi["utc"], format="isot",scale="utc")
    hist, _bin_edges = np.histogram(times.unix_tai, bins=bins_all)
    bin_edges = Time(_bin_edges, format="unix_tai",scale="utc").datetime

    time_support(format="unix_tai")

    fig = plt.figure(figsize=(12, 6))
    gs = gridspec.GridSpec(2, 1)

    gs_ax0 = fig.add_subplot(gs[0, 0])
    gs_ax0.stairs(hist_all, bin_edges_all, label="Whole file data")
    date_format = DateFormatter("%H:%M:%S")
    gs_ax0.xaxis.set_major_formatter(date_format)
    gs_ax0.axvline(x=times[0].datetime, color="purple", linestyle="-.", linewidth=1.5, label="Time select start")
    gs_ax0.axvline(x=times[-1].datetime, color="purple", linestyle="-.", linewidth=1.5, label="Time select end")
    gs_ax0.axvline(x=Time("2024-04-17 22:18:19",scale="utc").datetime, color="red", linestyle="--", linewidth=1.5, label="Shift to AR")
    gs_ax0.set_ylabel("Counts")
    gs_ax0.set_title(f"CdTe{CDTE_NUM} Time Profile (T-bin:{binSize} s)")
    plt.xticks(rotation=10, ha="right")

    gs_ax0.stairs(hist, bin_edges, label="Time selected data")
    gs_ax0.legend()

    gs_ax1 = fig.add_subplot(gs[1, 0])
    gs_ax1.stairs(hist_all/binSize, bin_edges_all, label="Whole file data")
    date_format = DateFormatter("%H:%M:%S")
    gs_ax1.xaxis.set_major_formatter(date_format)
    gs_ax1.axvline(x=times[0].datetime, color="purple", linestyle="-.", linewidth=1.5, label="Time select start")
    gs_ax1.axvline(x=times[-1].datetime, color="purple", linestyle="-.", linewidth=1.5, label="Time select end")
    gs_ax1.axvline(x=Time("2024-04-17 22:18:19",scale="utc").datetime, color="red", linestyle="--", linewidth=1.5, label="Shift to AR")
    gs_ax1.set_ylabel("Counts/s")
    gs_ax1.set_xlabel("UTC Time [2024-04-17]")
    plt.xticks(rotation=10, ha="right")

    gs_ax1.stairs(hist/binSize, bin_edges, label="Time selected data")
    gs_ax1.legend()

    plt.tight_layout()
    if SAVE_DIR is not None:
        plt.savefig(SAVE_DIR+f"CdTe{CDTE_NUM} {ENERGIES} Time Profile (T-bin:{binSize} s).png", bbox_inches="tight")
    plt.show()
    
    # %%
    # Spatial
    # ~~~~~~~
    #
    # Plot the data loaded in and start making observational selections

    fig = plt.figure(figsize=(18, 6))
    gs = gridspec.GridSpec(1, 3)

    gs_ax0 = fig.add_subplot(gs[0, 0])
    gs_ax0.hist2d(cdte_data_toi["solar_x_unaligned"], cdte_data_toi["solar_y_unaligned"], bins=150, norm=LogNorm())
    gs_ax0.set_xlim(-1000,1000)
    gs_ax0.set_ylim(-1000,1000)
    gs_ax0.set_xlabel("Solar-X [Arcsec]")
    gs_ax0.set_ylabel("Solar-Y [Arcsec]")
    gs_ax0.set_title(f"CdTe{CDTE_NUM}")

    # zoom plot
    x_zoom = (-550, -300)
    y_zoom = (-200, 100)
    gs_ax1 = fig.add_subplot(gs[0, 1])
    gs_ax1.hist2d(cdte_data_toi["solar_x_unaligned"], cdte_data_toi["solar_y_unaligned"], bins=150, norm=LogNorm())
    gs_ax1.set_xlim(*x_zoom)
    gs_ax1.set_ylim(*y_zoom)
    gs_ax1.set_xlabel("Solar-X [Arcsec]")
    gs_ax1.set_ylabel("Solar-Y [Arcsec]")
    gs_ax1.set_title(f"CdTe{CDTE_NUM} (zoomed)")

    bool_x = (x_zoom[0]<=cdte_data_toi["solar_x_unaligned"])&(cdte_data_toi["solar_x_unaligned"]<x_zoom[-1])
    bool_y = (y_zoom[0]<=cdte_data_toi["solar_y_unaligned"])&(cdte_data_toi["solar_y_unaligned"]<y_zoom[-1])
    cdte_data_roi = cdte_data_toi[bool_x & bool_y]

    # confirm selection
    gs_ax2 = fig.add_subplot(gs[0, 2])
    gs_ax2.hist2d(cdte_data_roi["solar_x_unaligned"], cdte_data_roi["solar_y_unaligned"], bins=50, norm=LogNorm())
    gs_ax2.set_xlim(-750, -50)
    gs_ax2.set_ylim(-450, 250)
    gs_ax2.set_xlabel("Solar-X [Arcsec]")
    gs_ax2.set_ylabel("Solar-Y [Arcsec]")
    gs_ax2.set_title(f"CdTe{CDTE_NUM} (zoomed) - Region selection confirmation")

    plt.tight_layout()
    if SAVE_DIR is not None:
        plt.savefig(SAVE_DIR+f"CdTe{CDTE_NUM} {ENERGIES} Image.png", bbox_inches="tight")
    plt.show()

    # %%
    # Spectral
    # ~~~~~~~~
    #
    # Coming from the future, I know the Count axis binning for the CdTe
    # detector has edges: ``np.arange(2.95, 30.1, 0.1)<<u.keV``
    #
    # Plot the spectrum from the event list I"m left with:

    ebins = np.arange(2.95, 30.1, 0.1)<<u.keV
    ebinning = np.diff(ebins)

    fig = plt.figure(figsize=(18, 6))
    gs = gridspec.GridSpec(1, 2)

    if ENERGIES=="DOI":
        pc, pb = np.histogram(cdte_data_roi["doi_corrected_energy"], bins=ebins.value)
    elif ENERGIES=="Pt-side":
        pc, pb = np.histogram(cdte_data_roi["pt_gap_loss_corrected_merged_energy_list"][:, 0], bins=ebins.value)
    
    pc = pc << u.ct
    counts_error = np.sqrt(pc.value) << u.ct
    mid_bins = (pb[:-1]+pb[1:])/2

    gs_ax0 = fig.add_subplot(gs[0, 0])
    # gs_ax0.stairs(pc, pb, color="b")
    gs_ax0.errorbar(mid_bins, pc, yerr=counts_error, color="b", ls="")
    gs_ax0.set_xlim(ebins.value[0], ebins.value[-1])
    gs_ax0.set_ylim(0, np.max(pc.value)*1.05)
    gs_ax0.set_xlabel(f"Count Energy [{ebins.unit:latex}]")
    gs_ax0.set_ylabel(f"Spectrum [{pc.unit:latex}]")
    gs_ax0.set_title(f"CdTe{CDTE_NUM} Count Spectrum")

    flux = pc/obs_time/ebinning
    flux_error = counts_error/obs_time/ebinning
    
    gs_ax1 = fig.add_subplot(gs[0, 1])
    # gs_ax1.stairs(flux, pb, color="b")
    gs_ax1.errorbar(mid_bins, flux, yerr=flux_error, color="b", ls="")
    gs_ax1.set_xlim(ebins.value[0], ebins.value[-1])
    gs_ax1.set_ylim(0, np.max(flux.value)*1.05)
    gs_ax1.set_xlabel(f"Count Energy [{ebins.unit:latex}]")
    gs_ax1.set_ylabel(f"Spectrum (no lvt corr.) [{flux.unit:latex}]")
    gs_ax1.set_title(f"CdTe{CDTE_NUM} Flux Spectrum (T-int:{obs_time:.2f}, E-bin:{ebinning[0]:.2f})")

    plt.tight_layout()
    if SAVE_DIR is not None:
        plt.savefig(SAVE_DIR+f"CdTe{CDTE_NUM} {ENERGIES} Spectrum.png", bbox_inches="tight")
    plt.show()

    # %%
    # What response info do I need?
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #
    # Now use the selected data to work out the response that I want.
    #
    # Need to work out the:
    #
    # * detector region (RMF) and off-axis angle (optics)
    # * use the time selection (for atmospheric attenuation in ARF)

    # %%
    # CdTe region
    # ~~~~~~~~~~~

    # response is defined on the Pt side pitch regions so use:
    pt_strips = cdte_data_roi["pt_merged_position_list"][:,0]
    regions = [get_hit_region(strip) for strip in pt_strips]

    reg_hist, reg_bin_edges = np.histogram(regions, bins=[-2.5, -1.5, -0.5, 0.5, 1.5, 2.5])
    reg_bin_mids = np.array([-2, -1, 0, 1, 2])

    valid_counts = np.nonzero(reg_bin_mids>=0) # should be all but be safe
    total_valid = np.sum(reg_hist[valid_counts])

    reg0_frac = reg_hist[np.nonzero(reg_bin_mids==0)][0]/total_valid
    reg1_frac = reg_hist[np.nonzero(reg_bin_mids==1)][0]/total_valid
    reg2_frac = reg_hist[np.nonzero(reg_bin_mids==2)][0]/total_valid

    guard_ring = reg_hist[np.nonzero(reg_bin_mids==-1)][0]
    nonphysical = reg_hist[np.nonzero(reg_bin_mids==-2)][0]

    fig = plt.figure(figsize=(16, 6))
    gs = gridspec.GridSpec(1, 1)

    gs_ax0 = fig.add_subplot(gs[0, 0])
    gs_ax0.stairs(reg_hist, reg_bin_edges, label="Region distribution")

    gs_ax0.set_xlabel("Region Number")
    gs_ax0.set_ylabel("Counts")
    gs_ax0.set_title(f"CdTe{CDTE_NUM} Count Hit Location")

    gs_ax0.annotate(f"Valid region 0 Frac.: {reg0_frac}\nValid region 1 Frac.: {reg1_frac}\nValid region 2 Frac.: {reg2_frac}\nGuard ring counts (not inc.): {guard_ring}\nNon-physical counts (not inc.): {nonphysical}", 
                    (0.05,0.95), 
                    xycoords="axes fraction",
                    va="top",
                    ha="left",
                    )

    # Mapping:
    # Non-physical    : -2
    # Guard ring      : -1
    # Pitch of 60 um  :  0
    # Pitch of 80 um  :  1
    # Pitch of 100 um :  2
    sec = gs_ax0.secondary_xaxis(location="top")
    sec.set_xticks(reg_bin_mids, labels=["Non-physical", "Guard ring ", "60 um Pitch", "80 um Pitch", "100 um Pitch"])
    plt.tight_layout()
    if SAVE_DIR is not None:
        plt.savefig(SAVE_DIR+f"CdTe{CDTE_NUM} {ENERGIES} Count Hit Location.png", bbox_inches="tight")
    plt.show()

    # %%
    # Off-axis angle
    # ~~~~~~~~~~~~~~

    # middle of the CdTe detector is "strip" 63.5
    det_centre = strip_ref_arcsec(63.5, allow_float=True)
    ave_pt_strip = strip_ref_arcsec(np.mean(pt_strips), allow_float=True)
    al_strips = cdte_data_roi["al_merged_position_list"][:,0]
    ave_al_strip = strip_ref_arcsec(np.mean(al_strips), allow_float=True)

    pt_loc, al_loc = ave_pt_strip-det_centre , ave_al_strip-det_centre 

    off_axis_angle = np.sqrt(pt_loc**2 + al_loc**2) << u.arcmin

    # %%
    # Time
    # ~~~~
    #
    # This will change when the atmsopheric attenuation function 
    # includes time calibration but just use `obs_time` from the start
    # of the observation for now.
    #
    # For the `response_tools.attenuation.att_foxsi4_atmosphere` 
    # function, the start of the observation is 100 s in so range is:

    time_range = [100, obs_time.value+100] << u.second
    livetime = np.sum(cdte_data_roi["livetime"])*1e-8 << u.second

    new_flux = flux*(obs_time/livetime) 
    new_flux_error = flux_error*(obs_time/livetime) 

    # %%
    # Getting the response info
    # ~~~~~~~~~~~~~~~~~~~~~~~~~
    #
    # Now I really was the response information
    # set up the ARF with the RMF information then make the SRM

    # just some simple mapping so the correct telescope functions are chosen from `CDTE_NUM`
    tel_num = {1:5, 2:3, 3:4, 4:2}[CDTE_NUM]
    rmf_func = {2:responses.foxsi4_telescope2_rmf, 
                3:responses.foxsi4_telescope3_rmf, 
                4:responses.foxsi4_telescope4_rmf, 
                5:responses.foxsi4_telescope5_rmf}
    arf_func = {2:responses.foxsi4_telescope2_flight_arf, 
                3:responses.foxsi4_telescope3_flight_arf, 
                4:responses.foxsi4_telescope4_flight_arf, 
                5:responses.foxsi4_telescope5_flight_arf}

    # if pitch is defined instead, pass as ``responses.foxsi4_telescope2_rmf(pitch=user_pitch)``
    pos_rmf0 = rmf_func[tel_num](region=0)
    pos_rmf1 = rmf_func[tel_num](region=1)
    pos_rmf2 = rmf_func[tel_num](region=2)

    pos_rmf = responses.Response2DOutput(filename="No-File", 
                                         function_path="No-Path",
                                         input_energy_edges=pos_rmf0.input_energy_edges,
                                         output_energy_edges=pos_rmf0.output_energy_edges,
                                         response=reg0_frac*pos_rmf0.response\
                                                    +reg1_frac*pos_rmf1.response\
                                                        +reg2_frac*pos_rmf2.response,
                                         response_type="Combined-RFM",
                                         telescope=f"RMF0:{pos_rmf0.telescope},RMF1:{pos_rmf1.telescope},RMF2:{pos_rmf2.telescope}",
                                         elements=(pos_rmf0,
                                                   pos_rmf1,
                                                   pos_rmf2,
                                                   ),
                                        )

    rmf_mid_energies = (pos_rmf.input_energy_edges[:-1]+pos_rmf.input_energy_edges[1:])/2

    pos_arf = arf_func[tel_num](mid_energies=rmf_mid_energies,
                                off_axis_angle=off_axis_angle,
                                time_range=time_range)
    
    pos_srm = responses.foxsi4_telescope_spectral_response(pos_arf, pos_rmf)

    # %%
    # plot the above
    fig = plt.figure(figsize=(18, 5))
    gs = gridspec.GridSpec(1, 3)

    # the ARF result (1D) and general plotting code
    gs_ax0 = fig.add_subplot(gs[0, 0])
    gs_ax0.plot(pos_arf.mid_energies, pos_arf.response)
    gs_ax0.set_xlabel(f"Photon Energy [{pos_arf.mid_energies.unit:latex}]")
    gs_ax0.set_ylabel(f"Response [{pos_arf.response.unit:latex}]")
    gs_ax0.set_title(f"{pos_arf.response_type}:: {pos_arf.telescope}")

    # the RMF result (2D) and general plotting code
    gs_ax1 = fig.add_subplot(gs[0, 1])
    r = gs_ax1.imshow(pos_rmf.response.value,
                    origin="lower",
                    norm=LogNorm(vmin=0.001),
                    extent=[np.min(pos_rmf.output_energy_edges.value),
                            np.max(pos_rmf.output_energy_edges.value),
                            np.min(pos_rmf.input_energy_edges.value),
                            np.max(pos_rmf.input_energy_edges.value)]
                    )
    cbar = plt.colorbar(r)
    cbar.ax.set_ylabel(f"Response [{pos_rmf.response.unit:latex}]")
    gs_ax1.set_xlabel(f"Count Energy [{pos_rmf.output_energy_edges.unit:latex}]")
    gs_ax1.set_ylabel(f"Photon Energy [{pos_rmf.input_energy_edges.unit:latex}]")
    gs_ax1.set_title(f"{pos_rmf.response_type}:: {pos_rmf.telescope}\nCombined Response (reg0:{reg0_frac:.2f}, reg1:{reg1_frac:.2f} reg2:{reg2_frac:.2f})")

    # the SRM result (2D) and general plotting code
    gs_ax2 = fig.add_subplot(gs[0, 2])
    r = gs_ax2.imshow(pos_srm.response.value,
                    origin="lower",
                    norm=LogNorm(vmin=0.001),
                    extent=[np.min(pos_srm.output_energy_edges.value),
                            np.max(pos_srm.output_energy_edges.value),
                            np.min(pos_srm.input_energy_edges.value),
                            np.max(pos_srm.input_energy_edges.value)]
                    )
    cbar = plt.colorbar(r)
    cbar.ax.set_ylabel(f"Response [{pos_srm.response.unit:latex}]")
    gs_ax2.set_xlabel(f"Count Energy [{pos_srm.output_energy_edges.unit:latex}]")
    gs_ax2.set_ylabel(f"Photon Energy [{pos_srm.input_energy_edges.unit:latex}]")
    gs_ax2.set_title(f"{pos_srm.response_type}:: {pos_srm.telescope}")

    plt.tight_layout()
    if SAVE_DIR is not None:
        plt.savefig(SAVE_DIR+f"CdTe{CDTE_NUM} {ENERGIES} responses.png", bbox_inches="tight")
    plt.show()

    # %%
    # Onto spectral fitting
    # ~~~~~~~~~~~~~~~~~~~~~
    #
    # import sunkit-spex

    spec_plot_size = (16, 12)
    spec_font_size = 20
    default_text = 10
    ylims = [1, 2e4]
    spec_single_plot_size = (8, 10)

    # %%
    # Set up the data
    # ~~~~~~~~~~~~~~~
    #
    # {"photon_channel_bins":Photon Space Bins (e.g., [keV,keV],[keV,keV],...]),
    # "photon_channel_mids":Photon Space Bin Mid-points (e.g., [keV,...]),
    # "photon_channel_binning":Photon Space Binwidths (e.g., [keV,...]),
    # "count_channel_bins":Count Space Bins (e.g., [keV,keV],[keV,keV],...]),
    # "count_channel_mids":Count Space Bin Mid-points (e.g., [keV,...]),
    # "count_channel_binning":Count Space Binwidths (e.g., [keV,...]),
    # "counts":counts (e.g., cts),
    # "count_error":Count Error for `counts`,
    # "count_rate":Count Rate (e.g., cts/keV/s),
    # "count_rate_error":Count Rate Error for `count_rate`,
    # "effective_exposure":Effective Exposure (e.g., s),
    # "srm":Spectral Response Matrix (e.g., cts/ph * cm^2),
    # "extras":{"any_extra_info":or_empty_dict}
    # };

    phot_edges = np.concatenate((pos_srm.input_energy_edges.value[:-1][:,None],pos_srm.input_energy_edges.value[1:][:,None]), axis=1)
    count_edges = np.concatenate((pos_srm.output_energy_edges.value[:-1][:,None],pos_srm.output_energy_edges.value[1:][:,None]), axis=1)
    clean_srm = pos_srm.response
    clean_srm[np.isnan(clean_srm)] = 0 << pos_srm.response.unit
    custom_dict = {"photon_channel_bins":phot_edges,
                    "photon_channel_mids":np.mean(phot_edges, axis=1),
                    "photon_channel_binning":np.diff(phot_edges).flatten(),
                    "count_channel_bins":count_edges,
                    "count_channel_mids":np.mean(count_edges, axis=1),
                    "count_channel_binning":np.diff(count_edges).flatten(),
                    "counts":pc.value,
                    "count_error":counts_error.value,
                    "count_rate":new_flux.value,
                    "count_rate_error":new_flux_error.value,
                    "effective_exposure":livetime.value,
                    "srm":clean_srm.value,
                    "extras":{}
                    };
    
    # %%
    # Let's see how model fold thrhough the response and compare to data
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #

    temperature0, emission_measure0, colour0 = 15.8 << u.MK, 3.1e46 << u.cm**-3, "royalblue"
    thermal_photons0 = thermal_emission(pos_srm.input_energy_edges, temperature0, emission_measure0, relative_abundances=((26, 0),))#, abundance_type="sun_photospheric")
    temperature1, emission_measure1, colour1 = 21 << u.MK, 2e46 << u.cm**-3, "slateblue"
    thermal_photons1 = thermal_emission(pos_srm.input_energy_edges, temperature1, emission_measure1)#, abundance_type="sun_photospheric")
    temperature2, emission_measure2, colour2 = 10 << u.MK, 2e48 << u.cm**-3, "blueviolet"
    thermal_photons2 = thermal_emission(pos_srm.input_energy_edges, temperature2, emission_measure2)#, abundance_type="sun_photospheric")

    fig = plt.figure(figsize=(18, 6))
    gs = gridspec.GridSpec(1, 3)

    gs_ax0 = fig.add_subplot(gs[0, 0])
    # gs_ax0.stairs(pc, pb, color="b")
    gs_ax0.plot(custom_dict["photon_channel_mids"], thermal_photons0, c=colour0, label=f"T0:{temperature0:latex}, EM0:{emission_measure0:latex}")
    gs_ax0.plot(custom_dict["photon_channel_mids"], thermal_photons1, c=colour1, label=f"T1:{temperature1:latex}, EM1:{emission_measure1:latex}")
    gs_ax0.plot(custom_dict["photon_channel_mids"], thermal_photons2, c=colour2, label=f"T2:{temperature2:latex}, EM2:{emission_measure2:latex}")
    gs_ax0.set_xlim(custom_dict["photon_channel_mids"][0], custom_dict["photon_channel_mids"][-1])
    gs_ax0.set_ylim(np.min(thermal_photons2.value)*0.95, np.max(thermal_photons2.value)*1.05)
    gs_ax0.set_xlabel(f"Photon Energy [{pos_srm.input_energy_edges.unit:latex}]")
    gs_ax0.set_ylabel(f"Spectrum [{thermal_photons0.unit:latex}]")
    gs_ax0.set_title("Photon Spectrum")
    gs_ax0.set_xscale("log")
    gs_ax0.set_yscale("log")
    plt.legend()

    thermal_model0 = thermal_photons0@clean_srm
    thermal_model1 = thermal_photons1@clean_srm
    thermal_model2 = thermal_photons2@clean_srm
    
    gs_ax1 = fig.add_subplot(gs[0, 1])
    # gs_ax1.stairs(pc, pb, color="b")
    gs_ax1.errorbar(mid_bins, pc, yerr=counts_error, color="k", ls="")
    gs_ax1.plot(custom_dict["photon_channel_mids"], thermal_model0*livetime*ebinning, c=colour0, label=f"T0:{temperature0:latex}, EM0:{emission_measure0:latex}")
    gs_ax1.plot(custom_dict["photon_channel_mids"], thermal_model1*livetime*ebinning, c=colour1, label=f"T1:{temperature1:latex}, EM1:{emission_measure1:latex}")
    gs_ax1.plot(custom_dict["photon_channel_mids"], thermal_model2*livetime*ebinning, c=colour2, label=f"T2:{temperature2:latex}, EM2:{emission_measure2:latex}")
    gs_ax1.set_xlim(ebins.value[0], ebins.value[-1])
    gs_ax1.set_ylim(5e-1, np.max(pc.value)*1.05)
    gs_ax1.set_xlabel(f"Count Energy [{ebins.unit:latex}]")
    gs_ax1.set_ylabel(f"Spectrum [{pc.unit:latex}]")
    gs_ax1.set_title(f"CdTe{CDTE_NUM} Count Spectrum")
    gs_ax1.set_xscale("log")
    gs_ax1.set_yscale("log")
    plt.legend()
    
    gs_ax2 = fig.add_subplot(gs[0, 2])
    # gs_ax2.stairs(flux, pb, color="b")
    gs_ax2.errorbar(mid_bins, new_flux, yerr=new_flux_error, color="k", ls="")
    gs_ax2.plot(custom_dict["photon_channel_mids"], thermal_model0, c=colour0, label=f"T0:{temperature0:latex}, EM0:{emission_measure0:latex}")
    gs_ax2.plot(custom_dict["photon_channel_mids"], thermal_model1, c=colour1, label=f"T1:{temperature1:latex}, EM1:{emission_measure1:latex}")
    gs_ax2.plot(custom_dict["photon_channel_mids"], thermal_model2, c=colour2, label=f"T2:{temperature2:latex}, EM2:{emission_measure2:latex}")
    # gs_ax2.set_xlim(ebins.value[0], ebins.value[-1])
    # gs_ax2.set_ylim(1e-1, np.max(new_flux.value)*1.05)
    gs_ax2.set_xlim(4, 20)
    gs_ax2.set_ylim(1e-3, 1e2)
    gs_ax2.set_xlabel(f"Count Energy [{ebins.unit:latex}]")
    gs_ax2.set_ylabel(f"Spectrum (lvt corr.) [{new_flux.unit:latex}]")
    gs_ax2.set_title(f"CdTe{CDTE_NUM} Flux Spectrum (eff. T-int:{livetime:.2f}, E-bin:{ebinning[0]:.2f})")
    gs_ax2.set_xscale("log")
    gs_ax2.set_yscale("log")
    plt.legend()

    plt.tight_layout()
    if SAVE_DIR is not None:
        plt.savefig(SAVE_DIR+f"CdTe{CDTE_NUM} {ENERGIES} Models.png", bbox_inches="tight")
    plt.show()

    # %%
    # Pass data to be fit
    # ~~~~~~~~~~~~~~~~~~~

    custom_spec = Fitter(custom_dict)

    # add the model to be used in fitting
    custom_spec.add_photon_model(f_fec_vth)

    # assign the fitting code's active model to be a combination of ones you defined
    custom_spec.model = "f_fec_vth"

    print(custom_spec.params)

    custom_spec.energy_fitting_range = [[4, 6], [9,20]]
    custom_spec.params["temp1_spectrum1"] = {"Value": 16, "Bounds": (10, 30)}  # units MK
    custom_spec.params["em461_spectrum1"] = {"Value": 3, "Bounds": (1e-1, 1e2)}  # units 1e46 cm^-3
    if VARY_FE:
        custom_spec.params["fea1_spectrum1"] = {"Value": 0.1, "Bounds": (1e-3, 1)}  #  fraction
    else:
        custom_spec.params["fea1_spectrum1"] = {"Status":"fix", "Value": 1, "Bounds": (1e-3, 1)} #  fraction

    if VARY_GAIN:
        custom_spec.rParams["gain_slope_spectrum1"] = {"Status":"free", "Bounds": (0.5, 1.5)}


    # %%
    # Fit the data
    # ~~~~~~~~~~~~

    print(custom_spec.params)

    minimised_params = custom_spec.fit(tol=1e-6)

    print(custom_spec.params)

    plt.rcParams["font.size"] = spec_font_size
    plt.figure(figsize=spec_single_plot_size)

    # the only line needed to plot the result
    axes, res_axes = custom_spec.plot()
    axes[0].set_xscale("log")
    axes[0].set_yscale("log")
    res_axes[0].set_xscale("log")

    _ta0 = " vary-fe" if VARY_FE else " fix-fe"
    _ta1 = " vary-gain" if VARY_GAIN else " fix-gain"

    plt.tight_layout()
    if SAVE_DIR is not None:
        plt.savefig(SAVE_DIR+f"CdTe{CDTE_NUM} {ENERGIES} Min-fit{_ta0}{_ta1}.png", bbox_inches="tight")
    plt.show()
    plt.rcParams["font.size"] = default_text


    mcmc_result = custom_spec.run_mcmc(steps_per_walker=1_000)
    custom_spec.burn_mcmc = 500

    plt.figure()
    custom_spec.plot_log_prob_chain()
    plt.tight_layout()
    if SAVE_DIR is not None:
        plt.savefig(SAVE_DIR+f"CdTe{CDTE_NUM} {ENERGIES} MCMC-fit Chain{_ta0}{_ta1}.png", bbox_inches="tight")
    plt.show()

    corner_plot = custom_spec.corner_mcmc()
    plt.tight_layout()
    if SAVE_DIR is not None:
        plt.savefig(SAVE_DIR+f"CdTe{CDTE_NUM} {ENERGIES} MCMC-fit Corner{_ta0}{_ta1}.png", bbox_inches="tight")
    plt.show()

    plt.rcParams["font.size"] = spec_font_size
    plt.figure(figsize=spec_single_plot_size)
    axes, res_axes = custom_spec.plot()
    plt.tight_layout()
    if SAVE_DIR is not None:
        plt.savefig(SAVE_DIR+f"CdTe{CDTE_NUM} {ENERGIES} MCMC-fit{_ta0}{_ta1}.png", bbox_inches="tight")
    
    plt.show()
    plt.rcParams["font.size"] = default_text

    print(custom_spec.params)