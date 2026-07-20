"""
Script showing a quick example on how the level 2 files can be used.

This will produce a plot containing images over the first flare pointing 
time over all energies then two energy ranges. Spectra from this time 
range will also be included.

THIS SCRIPT IS AN EXAMPLE, DO NOT USE DIRECTLY. If you would like to run 
this script and play around with it then either use this file and be 
aware of not adding/commiting any changes you make and/or just make a 
copy of this file and play around with the copy.
"""

import os

from astropy.table import Table
from astropy.time import Time
from astropy.visualization import time_support
import matplotlib.colors as colors
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np

import foxsi4_science_tools_py.cdte.filter_convert_extract as fce
from foxsi4_science_tools_py.cdte.image import (plot_image, 
                                                image_array)
from foxsi4_science_tools_py.cdte.strip_positions import strip_edges_micrometers
from foxsi4_science_tools_py.io.fits_tools import load_fits

# for plots
plt.rcParams["axes.labelsize"] = "small"
plt.rcParams["axes.titlesize"] = "medium"
plt.rcParams["figure.labelsize"] = "small"
plt.rcParams["figure.titlesize"] = "small"
plt.rcParams["font.size"] = 7
plt.rcParams['xtick.major.pad'] = 0.1
plt.rcParams['ytick.major.pad'] = 0.1

# do you want to save the figure?
SAVE_FIG = False

# point to level 2 files
check_dir = "/Users/kris/Documents/umnPostdoc/projects/both/foxsi4/foxsi4-analysis/data/cdte-levels/"
files = os.listdir(check_dir)
level2s = sorted([f for f in files if ("level2" in f) & (f.endswith(".fits"))])

# set up figure
fig = plt.figure(figsize=(12,10))
gs = gridspec.GridSpec(4, 4)

# loop through the 4 dets
for cdte_no in range(1,5):
    # load in the level 2 file
    cdte_hdul = load_fits(check_dir+level2s[cdte_no-1])
    cdte_data = Table(cdte_hdul[1].data)

    # time of interest for the flight, let's go for the first look at the flare
    ## If you have access to `cdte-tools-py`, you can see these TI time-tags 
    ## here : https://github.com/foxsi/cdte-tools/tree/main/cdte-tools-py/cdte_tools_py/pipeline#level-1
    ## or just read the meta data in `cdte_hdul[0].header`
    ti_obs_start, ti_obs_end = cdte_hdul[0].header["TI2"], cdte_hdul[0].header["TI3"]
    cdte_data_toi = fce.ti_range_filter(cdte_data, ti_obs_start, ti_obs_end)

    # only select nominal quality triggers
    cdte_data_ge = cdte_data_toi[cdte_data_toi["flag_quality"]==1]

    # plot image of triggers after the above filtering
    gs_ax = fig.add_subplot(gs[0, cdte_no-1])
    det_edges = strip_edges_micrometers()
    pt_single_clump_locs = cdte_data_ge["pt_merged_position_list"][:,0]
    al_single_clump_locs = cdte_data_ge["al_merged_position_list"][:,0]
    plot_image(image_array(pt_single_clump_locs, al_single_clump_locs), 
               strip_edges=det_edges.value, 
               rotation=cdte_hdul[0].header["FP_ROT"], 
               axes=gs_ax)
    gs_ax.set_title(f"CdTe {cdte_no}: All energies, first microflare pointing")
    gs_ax.set_ylabel(f"y [{det_edges.unit}]")
    gs_ax.set_xlabel(f"x [{det_edges.unit}]")
    # add spectrum histogram inset plot
    e_range = (1.5, 18)
    axins = gs_ax.inset_axes([0.52, 0.12, 0.47, 0.35], xlim=e_range, yticklabels=[])
    mide = 10
    pc, pb = np.histogram(cdte_data_ge["pt_gap_loss_corrected_merged_energy_list"][:,0], bins=100)
    ac, ab = np.histogram(cdte_data_ge["al_gap_loss_corrected_merged_energy_list"][:,0], bins=100)
    axins.stairs(pc, pb, label="Pt")
    axins.stairs(ac, ab, label="Al")
    axins.legend(fontsize=5)
    axins.axvline(mide, color="k", ls=":")
    axins.set_ylim([0, np.max([np.max(pc), np.max(ac)])*1.05])
    axins.set_xlabel("Energy [keV]")

    # try and see the images in a lower/upper energy range
    le, ue = [1.5, mide], [mide, 32]
    cdte_data_ge_l = fce.energy_range_filter(cdte_data_ge, *le)
    cdte_data_ge_u = fce.energy_range_filter(cdte_data_ge, *ue)

    # plot the lower energy range in helioprojective
    gs_ax = fig.add_subplot(gs[1, cdte_no-1])
    gs_ax.hist2d(cdte_data_ge_l["solar_x_unaligned"], cdte_data_ge_l["solar_y_unaligned"], bins=150, norm=colors.LogNorm())
    gs_ax.set_xlim(-1000,1000)
    gs_ax.set_ylim(-1000,1000)
    gs_ax.set_xlabel("Arcsec")
    gs_ax.set_ylabel("Arcsec")
    gs_ax.set_title(f"CdTe {cdte_no} with SPARCS info, E-range: {le}")
    gs_ax.set_aspect("equal")

    # plot the upper energy range in helioprojective
    gs_ax = fig.add_subplot(gs[2, cdte_no-1])
    gs_ax.hist2d(cdte_data_ge_u["solar_x_unaligned"], cdte_data_ge_u["solar_y_unaligned"], bins=150, norm=colors.LogNorm())
    gs_ax.set_xlim(-1000,1000)
    gs_ax.set_ylim(-1000,1000)
    gs_ax.set_xlabel("Arcsec")
    gs_ax.set_ylabel("Arcsec")
    gs_ax.set_title(f"CdTe {cdte_no} with SPARCS info, E-range: {ue}")
    gs_ax.set_aspect("equal")

    # spectrogram plot set-up
    times = Time(cdte_data_ge['utc'], format='isot',scale='utc')
    tbinSize = 5. # second
    binRange = [min(times.unix_tai),max(times.unix_tai)]
    ntBins = int((binRange[1] - binRange[0]) / tbinSize)
    tbins = np.arange(ntBins)*tbinSize+binRange[0]
    ebinSize = 0.2 # keV
    neBins = int((e_range[1] - e_range[0]) / ebinSize)
    ebins = np.arange(neBins)*ebinSize+e_range[0]
    H, xedges, yedges = np.histogram2d(times.unix_tai, 
                                       cdte_data_ge["pt_gap_loss_corrected_merged_energy_list"][:,0], 
                                       bins=(tbins, ebins))
    # spectrogram plot
    gs_ax = fig.add_subplot(gs[3, cdte_no-1])
    time_support(format='unix_tai')
    gs_ax.pcolormesh(Time(tbins, format='unix_tai',scale='utc').datetime, ebins, H.T,
                   norm=colors.LogNorm())
    plt.xticks(rotation=30, ha='right')
    gs_ax.axhline(mide, color="k", ls=":")
    gs_ax.set_xlabel("Time")
    gs_ax.set_ylabel("Energy [keV]")

plt.tight_layout()
if SAVE_FIG:
    plt.savefig(f"{check_dir}cdte-energy-ims.pdf", bbox_inches="tight")
plt.show()
