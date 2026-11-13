"""

"""

import astropy.units as u
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import numpy as np
from astropy.time import Time
from astropy.visualization import time_support
from matplotlib import gridspec

from foxsi4_science_tools_py.cdte import FOXSI4_CDTE_TO_TEL_MAPPING
from foxsi4_science_tools_py.cdte.response import construct_cdte_response
from foxsi4_science_tools_py.io.detectors import load_cdte, load_cmos
from foxsi4_science_tools_py.list_filters import by_general_region
from foxsi4_science_tools_py.util.clean import cdte_usual_suspects, cmos_usual_suspects
from foxsi4_science_tools_py.util.maps import draw_grid
from foxsi4_science_tools_py.foxsi_core import disk_map, full_disk_map, coalign_event_lists, event_data_xy_shift, count_time_profile_array
from foxsi4_science_tools_py.shapes import (
        general_shape, 
        shape_plot, 
        ellipse_shape,
        )
from foxsi4_science_tools_py.regions_of_interest import (
        foxsi4_main_loop_roi, 
        foxsi4_flux_rope_roi, 
        foxsi4_south_east_loop_roi, 
        foxsi4_bonus_microflare_roi,
        )


cdte_dir = "/Users/kris/Documents/umnPostdoc/projects/both/foxsi4/foxsi4-analysis/data/cdte-levels/"
cdte1_file = f"{cdte_dir}foxsi4_cdte1_level3_version2.fits"
cdte2_file = f"{cdte_dir}foxsi4_cdte2_level3_version2.fits"
cdte3_file = f"{cdte_dir}foxsi4_cdte3_level3_version2.fits"
cdte4_file = f"{cdte_dir}foxsi4_cdte4_level3_version2.fits"
# pre-energy calibration
# cmos2_file = "/Users/kris/Documents/umnPostdoc/projects/both/foxsi4/foxsi4-analysis-tools/cmos-tools-main/cmos-tools/FOXSI-4/event_process/EVENT_FILES/level2/foxsi4_telescope-1_CMOS-2_PHOTON_EVENT_LEVEL2_V26SEP06.fits"
cmos2_file = "/Users/kris/Downloads/foxsi4_cmos_eventlists_20261107/foxsi4_cmos2_events_20261107.fits"

# load in the data
cdte1 = load_cdte(cdte1_file)
cmos2 = load_cmos(cmos2_file)

# extract the event lists
cdte1_data = cdte1[1]
cmos2_data = cmos2[1]

# get similar structured data lists to work with
cdte1_fields = cdte_usual_suspects(cdte1_data, name="cdte1")
cmos2_fields = cmos_usual_suspects(cmos2_data, name="cmos2")


# create maps from the event list
m_cdte1_fd = full_disk_map(cdte1_fields)
m_cmos2_fd = full_disk_map(cmos2_fields, number_of_bins=600)
bottom_left = [-800, -500] << u.arcsec
top_right = [0, 300] << u.arcsec
m_cdte1 = disk_map(cdte1_fields, bottom_left, top_right)
m_cmos2 = disk_map(cmos2_fields, bottom_left, top_right, number_of_bins=600)

# plot the maps
fig = plt.figure(figsize=(11, 9))
gs = gridspec.GridSpec(2, 2)
axes0 = fig.add_subplot(gs[0, 0], projection=m_cdte1_fd)
axes1 = fig.add_subplot(gs[0, 1], projection=m_cmos2_fd)
axes2 = fig.add_subplot(gs[1, 0], projection=m_cdte1)
axes3 = fig.add_subplot(gs[1, 1], projection=m_cmos2)

m_cdte1_fd.plot(axes=axes0, vmin=1e0)
draw_grid(m_cdte1_fd, axes0)
plt.colorbar(fraction=0.035, pad=0.03, label=m_cdte1_fd.meta["PIXLUNIT"])
axes0.set_title("CdTe1: Full Disk")

m_cmos2_fd.plot(axes=axes1, vmin=1e0)
draw_grid(m_cmos2_fd, axes1)
plt.colorbar(fraction=0.035, pad=0.03, label=m_cmos2_fd.meta["PIXLUNIT"])
axes1.set_title("CMOS2: Full Disk")

m_cdte1.plot(axes=axes2, vmin=1e0)
draw_grid(m_cdte1, axes2)
plt.colorbar(fraction=0.035, pad=0.03, label=m_cdte1.meta["PIXLUNIT"])
axes2.set_title(f"CdTe1: bl {bottom_left}, tr {top_right}")

m_cmos2.plot(axes=axes3, vmin=1e0)
draw_grid(m_cmos2, axes3)
plt.colorbar(fraction=0.035, pad=0.03, label=m_cmos2.meta["PIXLUNIT"])
axes3.set_title(f"CMOS2: bl {bottom_left}, tr {top_right}")

plt.tight_layout()
plt.show()

## filter by shape
# we could define a random polygon shape
# vertices = [[-550, 130] << u.arcsec, [-400, 300] << u.arcsec, [-350, 150] << u.arcsec]
# shape = general_shape(vertices)
# but let"s just stick to pre-defined ones
shape = foxsi4_flux_rope_roi()
# shape = foxsi4_south_east_loop_roi()

cdte1_fields_shaped = by_general_region(cdte1_fields, shape)
cmos2_fields_shaped = by_general_region(cmos2_fields, shape)
m_cdte1_shaped = disk_map(cdte1_fields_shaped, bottom_left, top_right)
m_cmos2_shaped = disk_map(cmos2_fields_shaped, bottom_left, top_right, number_of_bins=600)


# other detectors
# load in the data
cdte2 = load_cdte(cdte2_file)
cdte3 = load_cdte(cdte3_file)
cdte4 = load_cdte(cdte4_file)
# extract the event lists
cdte2_data = cdte2[1]
cdte3_data = cdte3[1]
cdte4_data = cdte4[1]
# get similar structured data lists to work with
cdte2_fields = cdte_usual_suspects(cdte2_data, name="cdte2")
cdte3_fields = cdte_usual_suspects(cdte3_data, name="cdte3")
cdte4_fields = cdte_usual_suspects(cdte4_data, name="cdte4")
# filter by shape
cdte2_fields_shaped = by_general_region(cdte2_fields, shape)
cdte3_fields_shaped = by_general_region(cdte3_fields, shape)
cdte4_fields_shaped = by_general_region(cdte4_fields, shape)
# maps
m_cdte2 = disk_map(cdte2_fields, bottom_left, top_right)
m_cdte3 = disk_map(cdte3_fields, bottom_left, top_right)
m_cdte4 = disk_map(cdte4_fields, bottom_left, top_right)
# region filtered maps
m_cdte2_shaped = disk_map(cdte2_fields_shaped, bottom_left, top_right)
m_cdte3_shaped = disk_map(cdte3_fields_shaped, bottom_left, top_right)
m_cdte4_shaped = disk_map(cdte4_fields_shaped, bottom_left, top_right)


# plot the maps
fig = plt.figure(figsize=(15, 9))
gs = gridspec.GridSpec(2, 6)
axes0 = fig.add_subplot(gs[0, 0:2], projection=m_cdte1)
axes1 = fig.add_subplot(gs[0, 2:4], projection=m_cdte2)
axes2 = fig.add_subplot(gs[0, 4:], projection=m_cdte3)
axes3 = fig.add_subplot(gs[1, 0:2], projection=m_cdte4)
axes4 = fig.add_subplot(gs[1, 2:4], projection=m_cmos2)

m_cdte1.plot(axes=axes0, cmap="Greys", vmin=1e0)
m_cdte1_shaped.plot(axes=axes0, vmin=1e0)
shape_plot(shape, m_cdte1_shaped.coordinate_frame, color="tab:red")
draw_grid(m_cdte1, axes0)
axes0.set_title(f"CdTe1\nbl {bottom_left}, tr {top_right}")

m_cdte2.plot(axes=axes1, cmap="Greys", vmin=1e0)
m_cdte2_shaped.plot(axes=axes1, vmin=1e0)
shape_plot(shape, m_cdte2_shaped.coordinate_frame, color="tab:red")
draw_grid(m_cdte2, axes1)
axes1.set_title(f"CdTe2\nbl {bottom_left}, tr {top_right}")

m_cdte3.plot(axes=axes2, cmap="Greys", vmin=1e0)
m_cdte3_shaped.plot(axes=axes2, vmin=1e0)
shape_plot(shape, m_cdte3_shaped.coordinate_frame, color="tab:red")
draw_grid(m_cdte3, axes2)
axes2.set_title(f"CdTe3\nbl {bottom_left}, tr {top_right}")

m_cdte4.plot(axes=axes3, cmap="Greys", vmin=1e0)
m_cdte4_shaped.plot(axes=axes3, vmin=1e0)
shape_plot(shape, m_cdte4_shaped.coordinate_frame, color="tab:red")
draw_grid(m_cdte4, axes3)
axes3.set_title(f"CdTe4\nbl {bottom_left}, tr {top_right}")

m_cmos2.plot(axes=axes4, cmap="Greys", vmin=1e0)
m_cmos2_shaped.plot(axes=axes4, vmin=1e0)
shape_plot(shape, m_cmos2_shaped.coordinate_frame, color="tab:red")
draw_grid(m_cmos2, axes4)
axes4.set_title(f"CMOS2\nbl {bottom_left}, tr {top_right}")

# co-align
number_of_bins = 300
cdte1_fields_coaligned, shiftxy = coalign_event_lists(cmos2_fields, cdte1_fields, map_kwargs={"bottom_left":bottom_left, "top_right":top_right, "number_of_bins":number_of_bins})
# shape filter
cdte1_fields_coaligned_shaped = by_general_region(cdte1_fields_coaligned, shape)
# maps
m_cdte1_coaligned = disk_map(cdte1_fields_coaligned, bottom_left, top_right)
# region filtered maps
m_cdte1_coaligned_shaped = disk_map(cdte1_fields_coaligned_shaped, bottom_left, top_right)

axes5 = fig.add_subplot(gs[1, 4:], projection=m_cdte1_coaligned)

m_cdte1_coaligned.plot(axes=axes5, cmap="Greys", vmin=1e0)
m_cdte1_coaligned_shaped.plot(axes=axes5, vmin=1e0)
shape_plot(shape, m_cdte1_coaligned_shaped.coordinate_frame, color="tab:red")
draw_grid(m_cdte1_coaligned, axes5)
axes5.set_title(f"CdTe1 co-aligned to CMOS2\nbl {bottom_left}, tr {top_right}")

plt.tight_layout()
plt.show()


# %%
# Co-aligning all CdTe
# --------------------

number_of_bins = 300
cdte2_fields_coaligned, shiftxy2 = coalign_event_lists(cdte1_fields, cdte2_fields, map_kwargs={"bottom_left":bottom_left, "top_right":top_right, "number_of_bins":number_of_bins})
cdte3_fields_coaligned, shiftxy3 = coalign_event_lists(cdte1_fields, cdte3_fields, map_kwargs={"bottom_left":bottom_left, "top_right":top_right, "number_of_bins":number_of_bins})
cdte4_fields_coaligned, shiftxy3 = coalign_event_lists(cdte1_fields, cdte4_fields, map_kwargs={"bottom_left":bottom_left, "top_right":top_right, "number_of_bins":number_of_bins})
# shape filter
cdte2_fields_coaligned_shaped = by_general_region(cdte2_fields_coaligned, shape)
cdte3_fields_coaligned_shaped = by_general_region(cdte3_fields_coaligned, shape)
cdte4_fields_coaligned_shaped = by_general_region(cdte4_fields_coaligned, shape)
# maps
m_cdte2_coaligned = disk_map(cdte2_fields_coaligned, bottom_left, top_right)
m_cdte3_coaligned = disk_map(cdte3_fields_coaligned, bottom_left, top_right)
m_cdte4_coaligned = disk_map(cdte4_fields_coaligned, bottom_left, top_right)
# region filtered maps
m_cdte2_coaligned_shaped = disk_map(cdte2_fields_coaligned_shaped, bottom_left, top_right)
m_cdte3_coaligned_shaped = disk_map(cdte3_fields_coaligned_shaped, bottom_left, top_right)
m_cdte4_coaligned_shaped = disk_map(cdte4_fields_coaligned_shaped, bottom_left, top_right)

# plot the maps
fig = plt.figure(figsize=(15, 9))
gs = gridspec.GridSpec(2, 6)
axes0 = fig.add_subplot(gs[0, 0:2], projection=m_cdte1)
axes1 = fig.add_subplot(gs[0, 2:4], projection=m_cdte2_coaligned)
axes2 = fig.add_subplot(gs[0, 4:], projection=m_cdte3_coaligned)
axes3 = fig.add_subplot(gs[1, 0:2], projection=m_cdte4_coaligned)
axes4 = fig.add_subplot(gs[1, 2:4], projection=m_cmos2)

m_cdte1.plot(axes=axes0, cmap="Greys", vmin=1e0)
m_cdte1_shaped.plot(axes=axes0, vmin=1e0)
shape_plot(shape, m_cdte1_shaped.coordinate_frame, color="tab:red")
draw_grid(m_cdte1, axes0)
axes0.set_title(f"CdTe1\nbl {bottom_left}, tr {top_right}")

m_cdte2_coaligned.plot(axes=axes1, cmap="Greys", vmin=1e0)
m_cdte2_coaligned_shaped.plot(axes=axes1, vmin=1e0)
shape_plot(shape, m_cdte2_coaligned_shaped.coordinate_frame, color="tab:red")
draw_grid(m_cdte2_coaligned, axes1)
axes1.set_title(f"CdTe2 co-aligned to CdTe1\nbl {bottom_left}, tr {top_right}")

m_cdte3_coaligned.plot(axes=axes2, cmap="Greys", vmin=1e0)
m_cdte3_coaligned_shaped.plot(axes=axes2, vmin=1e0)
shape_plot(shape, m_cdte3_coaligned_shaped.coordinate_frame, color="tab:red")
draw_grid(m_cdte3_coaligned, axes2)
axes2.set_title(f"CdTe3 co-aligned to CdTe1\nbl {bottom_left}, tr {top_right}")

m_cdte4_coaligned.plot(axes=axes3, cmap="Greys", vmin=1e0)
m_cdte4_coaligned_shaped.plot(axes=axes3, vmin=1e0)
shape_plot(shape, m_cdte4_coaligned_shaped.coordinate_frame, color="tab:red")
draw_grid(m_cdte4_coaligned, axes3)
axes3.set_title(f"CdTe4 co-aligned to CdTe1\nbl {bottom_left}, tr {top_right}")

m_cmos2.plot(axes=axes4, cmap="Greys", vmin=1e0)
m_cmos2_shaped.plot(axes=axes4, vmin=1e0)
shape_plot(shape, m_cmos2_shaped.coordinate_frame, color="tab:red")
draw_grid(m_cmos2, axes4)
axes4.set_title(f"CMOS2\nbl {bottom_left}, tr {top_right}")

plt.tight_layout()
plt.show()

# %%
# Now co-aligning all CdTe to CMOS2
# ---------------------------------

obs_time_range = Time(["2024-04-17T22:14:40", "2024-04-17T22:20:41"])

cdte1_fields_coaligned_to_cmos, shiftxy1 = coalign_event_lists(cmos2_fields, cdte1_fields, map_kwargs={"bottom_left":bottom_left, "top_right":top_right, "number_of_bins":number_of_bins})
cdte2_fields_coaligned_to_cmos = event_data_xy_shift(cdte2_fields_coaligned, dx=shiftxy1[0], dy=shiftxy1[1])
cdte3_fields_coaligned_to_cmos = event_data_xy_shift(cdte3_fields_coaligned, dx=shiftxy1[0], dy=shiftxy1[1])
cdte4_fields_coaligned_to_cmos = event_data_xy_shift(cdte4_fields_coaligned, dx=shiftxy1[0], dy=shiftxy1[1])

# shape filter
cdte1_fields_coaligned_to_cmos_shaped = by_general_region(cdte1_fields_coaligned_to_cmos, shape)
cdte2_fields_coaligned_to_cmos_shaped = by_general_region(cdte2_fields_coaligned_to_cmos, shape)
cdte3_fields_coaligned_to_cmos_shaped = by_general_region(cdte3_fields_coaligned_to_cmos, shape)
cdte4_fields_coaligned_to_cmos_shaped = by_general_region(cdte4_fields_coaligned_to_cmos, shape)
# maps
m_cdte1_coaligned_to_cmos = disk_map(cdte1_fields_coaligned_to_cmos, bottom_left, top_right)
m_cdte2_coaligned_to_cmos = disk_map(cdte2_fields_coaligned_to_cmos, bottom_left, top_right)
m_cdte3_coaligned_to_cmos = disk_map(cdte3_fields_coaligned_to_cmos, bottom_left, top_right)
m_cdte4_coaligned_to_cmos = disk_map(cdte4_fields_coaligned_to_cmos, bottom_left, top_right)
# region filtered maps
m_cdte1_coaligned_to_cmos_shaped = disk_map(cdte1_fields_coaligned_to_cmos_shaped, bottom_left, top_right)
m_cdte2_coaligned_to_cmos_shaped = disk_map(cdte2_fields_coaligned_to_cmos_shaped, bottom_left, top_right)
m_cdte3_coaligned_to_cmos_shaped = disk_map(cdte3_fields_coaligned_to_cmos_shaped, bottom_left, top_right)
m_cdte4_coaligned_to_cmos_shaped = disk_map(cdte4_fields_coaligned_to_cmos_shaped, bottom_left, top_right)

# plot the maps
fig = plt.figure(figsize=(15, 9))
gs = gridspec.GridSpec(2, 6)
axes0 = fig.add_subplot(gs[0, 0:2], projection=m_cdte1_coaligned_to_cmos)
axes1 = fig.add_subplot(gs[0, 2:4], projection=m_cdte2_coaligned_to_cmos)
axes2 = fig.add_subplot(gs[0, 4:], projection=m_cdte3_coaligned_to_cmos)
axes3 = fig.add_subplot(gs[1, 0:2], projection=m_cdte4_coaligned_to_cmos)
axes4 = fig.add_subplot(gs[1, 2:4], projection=m_cmos2)

m_cdte1_coaligned_to_cmos.plot(axes=axes0, cmap="Greys", vmin=1e0)
m_cdte1_coaligned_to_cmos_shaped.plot(axes=axes0, vmin=1e0)
shape_plot(shape, m_cdte1_coaligned_to_cmos_shaped.coordinate_frame, color="tab:red")
draw_grid(m_cdte1_coaligned_to_cmos, axes0)
axes0.set_title(f"CdTe1 co-aligned to CMOS2\nbl {bottom_left}, tr {top_right}")

m_cdte2_coaligned_to_cmos.plot(axes=axes1, cmap="Greys", vmin=1e0)
m_cdte2_coaligned_to_cmos_shaped.plot(axes=axes1, vmin=1e0)
shape_plot(shape, m_cdte2_coaligned_to_cmos_shaped.coordinate_frame, color="tab:red")
draw_grid(m_cdte2_coaligned_to_cmos, axes1)
axes1.set_title(f"CdTe2 co-aligned to CdTe1 then CMOS2\nbl {bottom_left}, tr {top_right}")

m_cdte3_coaligned_to_cmos.plot(axes=axes2, cmap="Greys", vmin=1e0)
m_cdte3_coaligned_to_cmos_shaped.plot(axes=axes2, vmin=1e0)
shape_plot(shape, m_cdte3_coaligned_to_cmos_shaped.coordinate_frame, color="tab:red")
draw_grid(m_cdte3_coaligned_to_cmos, axes2)
axes2.set_title(f"CdTe3 co-aligned to CdTe1 then CMOS2\nbl {bottom_left}, tr {top_right}")

m_cdte4_coaligned_to_cmos.plot(axes=axes3, cmap="Greys", vmin=1e0)
m_cdte4_coaligned_to_cmos_shaped.plot(axes=axes3, vmin=1e0)
shape_plot(shape, m_cdte4_coaligned_to_cmos_shaped.coordinate_frame, color="tab:red")
draw_grid(m_cdte4_coaligned_to_cmos, axes3)
axes3.set_title(f"CdTe4 co-aligned to CdTe1 then CMOS2\nbl {bottom_left}, tr {top_right}")

m_cmos2.plot(axes=axes4, cmap="Greys", vmin=1e0)
m_cmos2_shaped.plot(axes=axes4, vmin=1e0)
shape_plot(shape, m_cmos2_shaped.coordinate_frame, color="tab:red")
draw_grid(m_cmos2, axes4)
axes4.set_title(f"CMOS2\nbl {bottom_left}, tr {top_right}")

# can add a time series too
axes5 = fig.add_subplot(gs[1, 4:])
cdte1_times, cdte1_counts = count_time_profile_array(cdte1_fields, time_binning=10<<u.second, start_time=obs_time_range[0], end_time=obs_time_range[1])
cdte2_times, cdte2_counts = count_time_profile_array(cdte2_fields, time_bins=cdte1_times)
cdte3_times, cdte3_counts = count_time_profile_array(cdte3_fields, time_bins=cdte1_times)
cdte4_times, cdte4_counts = count_time_profile_array(cdte4_fields, time_bins=cdte1_times)
cmos2_times, cmos2_counts = count_time_profile_array(cmos2_fields, time_bins=cdte1_times)
cdte1_times_cs, cdte1_counts_cs = count_time_profile_array(cdte1_fields_coaligned_to_cmos_shaped, time_bins=cdte1_times)
cdte2_times_cs, cdte2_counts_cs = count_time_profile_array(cdte2_fields_coaligned_to_cmos_shaped, time_bins=cdte1_times)
cdte3_times_cs, cdte3_counts_cs = count_time_profile_array(cdte3_fields_coaligned_to_cmos_shaped, time_bins=cdte1_times)
cdte4_times_cs, cdte4_counts_cs = count_time_profile_array(cdte4_fields_coaligned_to_cmos_shaped, time_bins=cdte1_times)
cmos2_times_cs, cmos2_counts_cs = count_time_profile_array(cmos2_fields_shaped, time_bins=cdte1_times)
time_support(format="unix_tai")
cdte1_col = "blue"
axes5.stairs(cdte1_counts, cdte1_times.plot_date, color=cdte1_col , label="CdTe1")
axes5.stairs(cdte1_counts_cs, cdte1_times_cs.plot_date, color=cdte1_col, ls=":")
cdte2_col = "red"
axes5.stairs(cdte2_counts, cdte2_times.plot_date, color=cdte2_col , label="CdTe2")
axes5.stairs(cdte2_counts_cs, cdte2_times_cs.plot_date, color=cdte2_col, ls=":")
cdte3_col = "green"
axes5.stairs(cdte3_counts, cdte3_times.plot_date, color=cdte3_col , label="CdTe3")
axes5.stairs(cdte3_counts_cs, cdte3_times_cs.plot_date, color=cdte3_col, ls=":")
cdte4_col = "purple"
axes5.stairs(cdte4_counts, cdte4_times.plot_date, color=cdte4_col , label="CdTe4")
axes5.stairs(cdte4_counts_cs, cdte4_times_cs.plot_date, color=cdte4_col, ls=":")
cmos2_col = "k"
axes5.stairs(cmos2_counts, cmos2_times.plot_date, color=cmos2_col , label="CMOS2")
axes5.stairs(cmos2_counts_cs, cmos2_times_cs.plot_date, color=cmos2_col, ls=":")
axes5.set_title(f"{cmos2_fields.name} time profile (region=\":\", all counts=\"-\")")
plt.xticks(rotation=30, ha="right")
axes5.set_ylabel(f"{cmos2_counts.unit:latex}")
axes5.set_xlabel(f"Start Time {cdte1_times[0]}")
fmt = mdates.DateFormatter("%H:%M:%S")
axes5.xaxis.set_major_formatter(fmt)
# v0 = axes5.axvline(Time(obs_time_range[0]).plot_date, ls="-.", c="k")
# v0 = axes5.axvline(Time(obs_time_range[1]).plot_date, ls="-.", c="k")
axes5.set_yscale("log")
axes5.set_xlim([cdte1_times.plot_date[0], cdte1_times.plot_date[-1]])
plt.legend()

plt.tight_layout()
plt.show()

# now construct the response
# get distance from rough FOV center to source region
cdte1_deltax = np.mean(cdte1_fields_coaligned_to_cmos.solar_x) - np.mean(cdte1_fields_coaligned_to_cmos_shaped.solar_x)
cdte1_deltay = np.mean(cdte1_fields_coaligned_to_cmos.solar_y) - np.mean(cdte1_fields_coaligned_to_cmos_shaped.solar_y)
cdte1_oa_to_source = np.sqrt((cdte1_deltax) ** 2 + (cdte1_deltay) ** 2) << u.arcmin
cdte2_deltax = np.mean(cdte2_fields_coaligned_to_cmos.solar_x) - np.mean(cdte2_fields_coaligned_to_cmos_shaped.solar_x)
cdte2_deltay = np.mean(cdte2_fields_coaligned_to_cmos.solar_y) - np.mean(cdte2_fields_coaligned_to_cmos_shaped.solar_y)
cdte2_oa_to_source = np.sqrt((cdte2_deltax) ** 2 + (cdte2_deltay) ** 2) << u.arcmin
cdte3_deltax = np.mean(cdte3_fields_coaligned_to_cmos.solar_x) - np.mean(cdte3_fields_coaligned_to_cmos_shaped.solar_x)
cdte3_deltay = np.mean(cdte3_fields_coaligned_to_cmos.solar_y) - np.mean(cdte3_fields_coaligned_to_cmos_shaped.solar_y)
cdte3_oa_to_source = np.sqrt((cdte3_deltax) ** 2 + (cdte3_deltay) ** 2) << u.arcmin
cdte4_deltax = np.mean(cdte4_fields_coaligned_to_cmos.solar_x) - np.mean(cdte4_fields_coaligned_to_cmos_shaped.solar_x)
cdte4_deltay = np.mean(cdte4_fields_coaligned_to_cmos.solar_y) - np.mean(cdte4_fields_coaligned_to_cmos_shaped.solar_y)
cdte4_oa_to_source = np.sqrt((cdte4_deltax) ** 2 + (cdte4_deltay) ** 2) << u.arcmin
cmos2_deltax = np.mean(cmos2_fields.solar_x) - np.mean(cmos2_fields_shaped.solar_x)
cmos2_deltay = np.mean(cmos2_fields.solar_y) - np.mean(cmos2_fields_shaped.solar_y)
cmos2_oa_to_source = np.sqrt((cmos2_deltax) ** 2 + (cmos2_deltay) ** 2) << u.arcmin
# now responses
oa_not_implemented = 0 << u.arcmin
cdte1_resp = construct_cdte_response(cdte1_fields_coaligned_to_cmos_shaped, 
                                     off_axis_angle=oa_not_implemented, 
                                     time_range=obs_time_range, 
                                     det_to_tel_mapping=FOXSI4_CDTE_TO_TEL_MAPPING)
cdte2_resp = construct_cdte_response(cdte2_fields_coaligned_to_cmos_shaped, 
                                     off_axis_angle=oa_not_implemented, 
                                     time_range=obs_time_range, 
                                     det_to_tel_mapping=FOXSI4_CDTE_TO_TEL_MAPPING)
cdte3_resp = construct_cdte_response(cdte3_fields_coaligned_to_cmos_shaped, 
                                     off_axis_angle=oa_not_implemented, 
                                     time_range=obs_time_range, 
                                     det_to_tel_mapping=FOXSI4_CDTE_TO_TEL_MAPPING)
cdte4_resp = construct_cdte_response(cdte4_fields_coaligned_to_cmos_shaped, 
                                     off_axis_angle=oa_not_implemented, 
                                     time_range=obs_time_range, 
                                     det_to_tel_mapping=FOXSI4_CDTE_TO_TEL_MAPPING)

################################################################################################################################################
# get CMOS response: this is new so no support elsewhere for this yet
cmos2_rmf_file = "/Users/kris/Documents/umnPostdoc/projects/both/foxsi4/foxsi4-analysis/venvs/cmos-env/foxsi4_cmos2_response_LEVEL2_V26OCT22.rmf"
# old line would be `cmos2_rmf = cmos_det_resp(file=cmos2_rmf_file, telescope=1)`
from astropy.io import fits
with fits.open(cmos2_rmf_file) as hdul:
    print(hdul[1].header)
    ph_e_lo = hdul[1].data["ENERG_LO"]
    ph_e_hi = hdul[1].data["ENERG_HI"]
    ct_e_lo = hdul[2].data["E_MIN"]
    ct_e_hi = hdul[2].data["E_MAX"]
    matrix, counts, energy = hdul[1].data["MATRIX"]<<(u.DN/u.ph), np.append(ct_e_lo, ct_e_hi[-1])<<u.keV, np.append(ph_e_lo, ph_e_hi[-1])<<u.keV 
from response_tools.detector_response import DetectorResponseOutput
from response_tools.responses import Response2DOutput, foxsi4_telescope_spectral_response, foxsi4_telescope1_arf, foxsi4_telescope1_flight_arf
_cmos2_rmf = DetectorResponseOutput(filename=cmos2_rmf_file,
                                   function_path="None",
                                   input_energy_edges=energy,
                                   output_energy_edges=counts,
                                   detector_response=matrix,
                                   detector="CMOS2-Detector-Response"
                                   )
cmos2_rmf = Response2DOutput(filename="No-File",
                             function_path="",
                             input_energy_edges=_cmos2_rmf.input_energy_edges,
                             output_energy_edges=_cmos2_rmf.output_energy_edges,
                             response=_cmos2_rmf.detector_response,
                             response_type="RMF",
                             telescope="foxsi4-telescope1",
                             elements=(_cmos2_rmf,
                                       ),
                             )
cmos2_mid_energies = (cmos2_rmf.input_energy_edges[:-1]+cmos2_rmf.input_energy_edges[1:])/2
# print("cmos2_mid_energies", cmos2_mid_energies,  cmos2_oa_to_source, obs_time_range)
cmos2_arf = foxsi4_telescope1_flight_arf(cmos2_mid_energies, oa_not_implemented, obs_time_range)
cmos2_resp = foxsi4_telescope_spectral_response(cmos2_arf, cmos2_rmf)
# plt.figure()
# plt.plot(cmos2_arf.mid_energies, cmos2_arf.response)
# plt.show()
# plt.figure()
# plt.imshow(_cmos2_rmf.detector_response.value, origin="lower", extent=[counts.value[0], counts.value[-1], energy.value[0], energy.value[-1]])
# plt.show()
# plt.figure()
# plt.imshow(cmos2_resp.response.value, origin="lower", extent=[counts.value[0], counts.value[-1], energy.value[0], energy.value[-1]])
# plt.show()
################################################################################################################################################

## Spectral fitting
# The SRMS might need cleaning first
cdte1_resp.response[np.isnan(cdte1_resp.response)] = 0 << cdte1_resp.response.unit
cdte2_resp.response[np.isnan(cdte2_resp.response)] = 0 << cdte2_resp.response.unit
cdte3_resp.response[np.isnan(cdte3_resp.response)] = 0 << cdte3_resp.response.unit
cdte4_resp.response[np.isnan(cdte4_resp.response)] = 0 << cdte4_resp.response.unit
cmos2_resp.response[np.isnan(cmos2_resp.response)] = 0 << cmos2_resp.response.unit
# CdTe
# all CdTe will share some info
cdte_phot_edges = np.concatenate((cdte1_resp.input_energy_edges.value[:-1][:,None], cdte1_resp.input_energy_edges.value[1:][:,None]), axis=1)
cdte_count_edges = np.concatenate((cdte1_resp.output_energy_edges.value[:-1][:,None], cdte1_resp.output_energy_edges.value[1:][:,None]), axis=1)
cdte_count_width = np.diff(cdte_count_edges).flatten()
cdte_count_edges_for_hist = cdte1_resp.output_energy_edges.value
# some items will not be shared
_cdte1_pc, _ = np.histogram(cdte1_fields_coaligned_to_cmos_shaped.signal.value, bins=cdte_count_edges_for_hist)
cdte1_pc = _cdte1_pc << u.ct
cdte1_counts_error = np.sqrt(cdte1_pc.value) << u.ct
cdte1_lvt = np.sum(cdte1_fields_coaligned_to_cmos_shaped.extra_count_info["livetime"])
cdte1_flux = cdte1_pc/cdte1_lvt/cdte_count_width
cdte1_flux_error = cdte1_counts_error/cdte1_lvt/cdte_count_width
cdte1_custom_dict = {"photon_channel_bins":cdte_phot_edges,
                     "photon_channel_mids":np.mean(cdte_phot_edges, axis=1),
                     "photon_channel_binning":np.diff(cdte_phot_edges).flatten(),
                     "count_channel_bins":cdte_count_edges,
                     "count_channel_mids":np.mean(cdte_count_edges, axis=1),
                     "count_channel_binning":cdte_count_width,
                     "counts":cdte1_pc.value,
                     "count_error":cdte1_counts_error.value,
                     "count_rate":cdte1_flux.value,
                     "count_rate_error":cdte1_flux_error.value,
                     "effective_exposure":cdte1_lvt.value,
                     "srm":cdte1_resp.response.value,
                     "extras":{}
                     };
_cdte2_pc, _ = np.histogram(cdte2_fields_coaligned_to_cmos_shaped.signal.value, bins=cdte_count_edges_for_hist)
cdte2_pc = _cdte2_pc << u.ct
cdte2_counts_error = np.sqrt(cdte2_pc.value) << u.ct
cdte2_lvt = np.sum(cdte2_fields_coaligned_to_cmos_shaped.extra_count_info["livetime"])
cdte2_flux = cdte2_pc/cdte2_lvt/cdte_count_width
cdte2_flux_error = cdte2_counts_error/cdte2_lvt/cdte_count_width
cdte2_custom_dict = {"photon_channel_bins":cdte_phot_edges,
                     "photon_channel_mids":np.mean(cdte_phot_edges, axis=1),
                     "photon_channel_binning":np.diff(cdte_phot_edges).flatten(),
                     "count_channel_bins":cdte_count_edges,
                     "count_channel_mids":np.mean(cdte_count_edges, axis=1),
                     "count_channel_binning":cdte_count_width,
                     "counts":cdte2_pc.value,
                     "count_error":cdte2_counts_error.value,
                     "count_rate":cdte2_flux.value,
                     "count_rate_error":cdte2_flux_error.value,
                     "effective_exposure":cdte2_lvt.value,
                     "srm":cdte2_resp.response.value,
                     "extras":{}
                     };
_cdte3_pc, _ = np.histogram(cdte3_fields_coaligned_to_cmos_shaped.signal.value, bins=cdte_count_edges_for_hist)
cdte3_pc = _cdte3_pc << u.ct
cdte3_counts_error = np.sqrt(cdte3_pc.value) << u.ct
cdte3_lvt = np.sum(cdte3_fields_coaligned_to_cmos_shaped.extra_count_info["livetime"])
cdte3_flux = cdte3_pc/cdte3_lvt/cdte_count_width
cdte3_flux_error = cdte3_counts_error/cdte3_lvt/cdte_count_width
cdte3_custom_dict = {"photon_channel_bins":cdte_phot_edges,
                     "photon_channel_mids":np.mean(cdte_phot_edges, axis=1),
                     "photon_channel_binning":np.diff(cdte_phot_edges).flatten(),
                     "count_channel_bins":cdte_count_edges,
                     "count_channel_mids":np.mean(cdte_count_edges, axis=1),
                     "count_channel_binning":cdte_count_width,
                     "counts":cdte3_pc.value,
                     "count_error":cdte3_counts_error.value,
                     "count_rate":cdte3_flux.value,
                     "count_rate_error":cdte3_flux_error.value,
                     "effective_exposure":cdte3_lvt.value,
                     "srm":cdte3_resp.response.value,
                     "extras":{}
                     };
_cdte4_pc, _ = np.histogram(cdte4_fields_coaligned_to_cmos_shaped.signal.value, bins=cdte_count_edges_for_hist)
cdte4_pc = _cdte4_pc << u.ct
cdte4_counts_error = np.sqrt(cdte4_pc.value) << u.ct
cdte4_lvt = np.sum(cdte4_fields_coaligned_to_cmos_shaped.extra_count_info["livetime"])
cdte4_flux = cdte4_pc/cdte4_lvt/cdte_count_width
cdte4_flux_error = cdte4_counts_error/cdte4_lvt/cdte_count_width
cdte4_custom_dict = {"photon_channel_bins":cdte_phot_edges,
                     "photon_channel_mids":np.mean(cdte_phot_edges, axis=1),
                     "photon_channel_binning":np.diff(cdte_phot_edges).flatten(),
                     "count_channel_bins":cdte_count_edges,
                     "count_channel_mids":np.mean(cdte_count_edges, axis=1),
                     "count_channel_binning":cdte_count_width,
                     "counts":cdte4_pc.value,
                     "count_error":cdte4_counts_error.value,
                     "count_rate":cdte4_flux.value,
                     "count_rate_error":cdte4_flux_error.value,
                     "effective_exposure":cdte4_lvt.value,
                     "srm":cdte4_resp.response.value,
                     "extras":{}
                     };
# CMOS2
cmos_phot_edges = np.concatenate((cmos2_resp.input_energy_edges.value[:-1][:,None], cmos2_resp.input_energy_edges.value[1:][:,None]), axis=1)
cmos_count_edges = np.concatenate((cmos2_resp.output_energy_edges.value[:-1][:,None], cmos2_resp.output_energy_edges.value[1:][:,None]), axis=1)
cmos_count_width = np.diff(cmos_count_edges).flatten()
_cmos2_pc, _ = np.histogram(cmos2_fields_shaped.signal.value, bins=cmos2_resp.output_energy_edges.value)
cmos2_pc = _cmos2_pc << u.ct
cmos2_counts_error = np.sqrt(cmos2_pc.value) << u.ct
cmos2_eff_exp = np.array((cmos2_fields_shaped.utc[-1]-cmos2_fields_shaped.utc[0]).sec) << u.second
cmos2_flux = cmos2_pc/cmos2_eff_exp/cmos_count_width
cmos2_flux_error = cmos2_counts_error/cmos2_eff_exp/cmos_count_width

# plt.figure()
# plt.stairs(_cmos2_pc, cmos2_resp.output_energy_edges.value)
# plt.show()
# plt.figure()
# plt.errorbar(np.mean(cmos_count_edges, axis=1), cmos2_flux.value, yerr=cmos2_flux_error.value)
# plt.show()
# plt.figure()
# plt.imshow(cmos2_resp.response.value, origin="lower", extent=[cmos_count_edges[0,0], cmos_count_edges[1,-1], cmos_phot_edges[0,0], cmos_phot_edges[1,-1]])
# plt.show()
cmos2_custom_dict = {"photon_channel_bins":cmos_phot_edges,
                     "photon_channel_mids":np.mean(cmos_phot_edges, axis=1),
                     "photon_channel_binning":np.diff(cmos_phot_edges).flatten(),
                     "count_channel_bins":cmos_count_edges,
                     "count_channel_mids":np.mean(cmos_count_edges, axis=1),
                     "count_channel_binning":cmos_count_width,
                     "counts":cmos2_pc.value,
                     "count_error":cmos2_counts_error.value,
                     "count_rate":cmos2_flux.value,
                     "count_rate_error":cmos2_flux_error.value,
                     "effective_exposure":cmos2_eff_exp.value,
                     "srm":cmos2_resp.response.value,
                     "extras":{}
                     };


# print(cmos2_custom_dict)

import os
import pathlib
import warnings
from numpy.exceptions import VisibleDeprecationWarning
from sunkit_spex.legacy.fitting.fitter import Fitter
from sunkit_spex.legacy import thermal
warnings.filterwarnings("ignore", category=RuntimeWarning)
try:
    warnings.filterwarnings("ignore", category=VisibleDeprecationWarning)
except AttributeError:
    warnings.filterwarnings("ignore", category=np.exceptions.VisibleDeprecationWarning)
try:
    thermal.CONTINUUM_GRID = thermal.setup_continuum_parameters(os.path.join(pathlib.Path(__file__).parent, 
                                                                             "..",
                                                                             "foxsi4_science_tools_py",
                                                                             "data",
                                                                             "thermal_chianti_grid", 
                                                                             "chianti_cont_01_250_unity_v101_t41_kris.sav"))
    # my own file just extends the thermal grid used by the code so it doesn't break
    # for consistency, could change the file below to #"chianti_lines_01_12_unity_v101_t41_kris.sav"
    # but it doesn't look like it changes much
    thermal.LINE_GRID = thermal.setup_line_parameters(os.path.join(pathlib.Path(__file__).parent, 
                                                                   "..",
                                                                   "foxsi4_science_tools_py",
                                                                   "data",
                                                                   "thermal_chianti_grid", 
                                                                   "chianti_lines_01_12_unity_v901_t41.sav")) 
except FileNotFoundError:
    warnings.warn("Sunkit-spex thermal CHIANTI grid files do not exist in the correct folder.")

cdte1_custom_spec = Fitter(cdte1_custom_dict)
cdte1_custom_spec.model = "f_vth"
cdte1_custom_spec.energy_fitting_range = [[4,20]]
cdte1_custom_spec.params["T1_spectrum1"] = {"Value": 10, "Bounds": (5, 30)}  # units MK
cdte1_custom_spec.params["EM1_spectrum1"] = {"Value": 1e3, "Bounds": (1e1, 5e3)}  # units 1e46 cm^-3 (1e-1, 4e2)

cdte1_minimised_params = cdte1_custom_spec.fit(tol=1e-6)

print(cdte1_custom_spec.params)

plt.figure(figsize=(7,9))
# the only line needed to plot the result
axes, res_axes = cdte1_custom_spec.plot()
axes[0].set_xscale("log")
axes[0].set_yscale("log")
res_axes[0].set_xscale("log")
plt.tight_layout()
plt.show()

# Cmos
# all_custom_spec = Fitter(cdte1_custom_dict,
#                            cdte2_custom_dict,
#                            cdte3_custom_dict,
#                            cdte4_custom_dict,
#                            cmos2_custom_dict,
#                            )
all_custom_spec = Fitter(cmos2_custom_dict,
                        #    cdte2_custom_dict,
                        #    cdte3_custom_dict,
                        #    cdte4_custom_dict,
                        #    cmos2_custom_dict,
                           )
all_custom_spec.model = "C*(f_vth)"

print(all_custom_spec.params)

all_custom_spec.energy_fitting_range = [0.9,9]
all_custom_spec.params["T1_spectrum1"] = {"Value": 10, "Bounds": (5, 20)}  # units MK
all_custom_spec.params["EM1_spectrum1"] = {"Value": 1e0, "Bounds": (1e-1, 5e1)}  # units 1e46 cm^-3 (1e-1, 4e2)
# all_custom_spec.params["T2_spectrum1"] = {"Value": 20, "Bounds": (10, 30)}  # units MK
# all_custom_spec.params["EM2_spectrum1"] = {"Value": 1e-1, "Bounds": (1e-3, 1e0)}  # units 1e46 cm^-3 (1e-1, 4e2)
all_custom_spec.params["C_spectrum1"] = {"Status":"frozen"} 
# all_custom_spec.params["C_spectrum2"] = {"Status":"free", "Bounds":(0.1, 3)} 
# all_custom_spec.params["C_spectrum3"] = {"Status":"free", "Bounds":(0.1, 3)} 
# all_custom_spec.params["C_spectrum4"] = {"Status":"free", "Bounds":(0.1, 3)} 
# all_custom_spec.params["C_spectrum5"] = {"Status":"free", "Bounds":(0.1, 3)} 

cdte1_minimised_params = all_custom_spec.fit(tol=1e-6)

print(all_custom_spec.params)

plt.figure(figsize=(12,9))
# the only line needed to plot the result
axes, res_axes = all_custom_spec.plot()
axes[0].set_xscale("log")
axes[0].set_yscale("log")
res_axes[0].set_xscale("log")
plt.tight_layout()
plt.show()


all_custom_spec = Fitter(cdte1_custom_dict,
                           cdte2_custom_dict,
                           cdte3_custom_dict,
                           cdte4_custom_dict,
                           cmos2_custom_dict,
                           )
# all_custom_spec.model = "C*(f_vth+f_vth)"
all_custom_spec.model = "C*(f_vth)"

print(all_custom_spec.params)

all_custom_spec.energy_fitting_range = {"spectrum1":[6, 15],
                                        "spectrum2":[6, 15],
                                        "spectrum3":[6, 15],
                                        "spectrum4":[6, 15],
                                        "spectrum5":[0.9,9]}
all_custom_spec.params["T1_spectrum1"] = {"Value": 10, "Bounds": (6, 30)}  # units MK
all_custom_spec.params["EM1_spectrum1"] = {"Value": 1e2, "Bounds": (1e0, 2e2)}  # units 1e46 cm^-3 (1e-1, 4e2)
# all_custom_spec.params["T2_spectrum1"] = {"Value": 20, "Bounds": (5, 30)}  # units MK
# all_custom_spec.params["EM2_spectrum1"] = {"Value": 1e-1, "Bounds": (1e-3, 1e1)}  # units 1e46 cm^-3 (1e-1, 4e2)
all_custom_spec.params["C_spectrum1"] = {"Status":"free", "Value":1e6, "Bounds":(1e4, 1e7)} 
all_custom_spec.params["C_spectrum2"] = {"Status":"free", "Value":1e6, "Bounds":(1e4, 1e7)} 
all_custom_spec.params["C_spectrum3"] = {"Status":"free", "Value":1e6, "Bounds":(1e4, 1e7)} 
all_custom_spec.params["C_spectrum4"] = {"Status":"free", "Value":1e6, "Bounds":(1e4, 1e7)} 
all_custom_spec.params["C_spectrum5"] = {"Status":"fix", "Value":1} 

cdte1_minimised_params = all_custom_spec.fit(tol=1e-6)

print(all_custom_spec.params)

plt.figure(figsize=(12,9))
# the only line needed to plot the result
axes, res_axes = all_custom_spec.plot()
# axes[0].set_xscale("log")
# axes[0].set_yscale("log")
# res_axes[0].set_xscale("log")
axes[-1].set_xlim([0.1, 15])
axes[-1].set_ylim([0.1, 4e2])
res_axes[-1].set_xlim([0.1, 15])
plt.tight_layout()
plt.show()

mcmc_result = all_custom_spec.run_mcmc(steps_per_walker=10000)

all_custom_spec.burn_mcmc = 1000

plt.figure()
all_custom_spec.plot_log_prob_chain()
plt.show()

corner_plot = all_custom_spec.corner_mcmc()

plt.figure(figsize=(12,9))
axes, res_axes = all_custom_spec.plot()
axes[-1].set_xlim([0.1, 15])
axes[-1].set_ylim([0.1, 4e2])
res_axes[-1].set_xlim([0.1, 15])
plt.show()