"""
WHAT'S HERE?
============

This example shows using a bunch of functions in the repository to work
with CdTe and CMOS data.

The topics covered are:

1. Loading in the files
2. Plotting maps and submaps
3. Filtering by shapes and co-alignment
4. Now showing time profiles
5. Filtering by more shapes

Most of the code here is plotting code as usual. To obtain the actual 
products, like event lists, etc., only takes a few lines.

A USER SHOULD CHANGE
--------------------
The following variables located roughly after the import statements 
below:

- `cdte_dir`
- `cdte1_file`
- `cmos2_file`

The code needs to find the files you want.

FOXSI DATA CAN BE FOUND
-----------------------

`FOXSI/FOXSI-4 Flight 36.370/Flight Science/CdTe data analysis/team_level_data` on the Google Drive
`FOXSI/FOXSI-4 Flight 36.370/Flight Science/CMOS data analysis/team_level_data` on the Google Drive
"""

import astropy.units as u
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
from astropy.time import Time
from astropy.visualization import time_support
from matplotlib import gridspec

from foxsi4_science_tools_py.list_filters import by_general_region
from foxsi4_science_tools_py.util.maps import draw_grid
from foxsi4_science_tools_py.foxsi_core import (
    cdte_events,
    cmos_events,
    disk_map, 
    full_disk_map, 
    coalign_event_lists, 
    count_time_profile_array,
    )
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


# %%
# 1. Loading in the files
# -----------------------

# Find CdTe data at `FOXSI/FOXSI-4 Flight 36.370/Flight Science/CdTe data analysis/team_level_data` on the Google Drive
cdte_dir = "/Users/kris/Documents/umnPostdoc/projects/both/foxsi4/foxsi4-analysis/data/cdte-levels/"
cdte1_file = f"{cdte_dir}foxsi4_cdte1_level3_version2.fits"
# Find CMOS data at `FOXSI/FOXSI-4 Flight 36.370/Flight Science/CMOS data analysis/team_level_data` on the Google Drive
cmos2_file = "/Users/kris/Downloads/foxsi4_cmos_eventlists_20261107/foxsi4_cmos2_events_20261107.fits"

# get similar structured data lists to work with
cdte1_events = cdte_events(cdte1_file, name="cdte1")
cmos2_events = cmos_events(cmos2_file, name="cmos2")


# %%
# 2. Plotting maps and submaps
# ----------------------------

# create maps from the event list
m_cdte1_fd = full_disk_map(cdte1_events)
m_cmos2_fd = full_disk_map(cmos2_events, number_of_bins=600)
bottom_left = [-800, -500] << u.arcsec
top_right = [0, 300] << u.arcsec
m_cdte1 = disk_map(cdte1_events, bottom_left, top_right)
m_cmos2 = disk_map(cmos2_events, bottom_left, top_right, number_of_bins=600)

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

# %%
# 3. Filtering by shapes and co-alignment
# ---------------------------------------
#
# The co-alignment only needs to happen once if the user wants.
#
# The `coalign_event_lists` function returns the shift found to align 
# the second map with the first. 
# 
# This shift can then be applied to another map using the 
# `foxsi4_science_tools_py.foxsi_core.event_data_xy_shift` function 
# which will return the original maps shifted.

## filter by shape
# let's look at the flux rope region
shape = foxsi4_flux_rope_roi()

cdte1_events_shaped = by_general_region(cdte1_events, shape)
cmos2_events_shaped = by_general_region(cmos2_events, shape)
m_cdte1_shaped = disk_map(cdte1_events_shaped, bottom_left, top_right)
m_cmos2_shaped = disk_map(cmos2_events_shaped, bottom_left, top_right, number_of_bins=600)

# plot the maps
fig = plt.figure(figsize=(15, 6))
gs = gridspec.GridSpec(1, 3)
axes0 = fig.add_subplot(gs[0, 0], projection=m_cdte1)
axes1 = fig.add_subplot(gs[0, 1], projection=m_cmos2)

m_cdte1.plot(axes=axes0, cmap="Greys", vmin=1e0)
m_cdte1_shaped.plot(axes=axes0, vmin=1e0)
shape_plot(shape, m_cdte1_shaped.coordinate_frame, color="tab:red")
draw_grid(m_cdte1, axes0)
axes0.set_title(f"CdTe1\nbl {bottom_left}, tr {top_right}")

m_cmos2.plot(axes=axes1, cmap="Greys", vmin=1e0)
m_cmos2_shaped.plot(axes=axes1, vmin=1e0)
shape_plot(shape, m_cmos2_shaped.coordinate_frame, color="tab:red")
draw_grid(m_cmos2, axes1)
axes1.set_title(f"CMOS2\nbl {bottom_left}, tr {top_right}")

# co-align
number_of_bins = 300
cdte1_events_coaligned, shiftxy = coalign_event_lists(cmos2_events, cdte1_events, map_kwargs={"bottom_left":bottom_left, "top_right":top_right, "number_of_bins":number_of_bins})
# shape filter
cdte1_events_coaligned_shaped = by_general_region(cdte1_events_coaligned, shape)
# maps
m_cdte1_coaligned = disk_map(cdte1_events_coaligned, bottom_left, top_right)
# region filtered maps
m_cdte1_coaligned_shaped = disk_map(cdte1_events_coaligned_shaped, bottom_left, top_right)

axes2 = fig.add_subplot(gs[0, 2], projection=m_cdte1_coaligned)

m_cdte1_coaligned.plot(axes=axes2, cmap="Greys", vmin=1e0)
m_cdte1_coaligned_shaped.plot(axes=axes2, vmin=1e0)
shape_plot(shape, m_cdte1_coaligned_shaped.coordinate_frame, color="tab:red")
draw_grid(m_cdte1_coaligned, axes2)
axes2.set_title(f"CdTe1 co-aligned to CMOS2\nShift-xy: {shiftxy}")

plt.tight_layout()
plt.show()

# %%
# 4. Now showing time profiles
# ----------------------------

# plot the maps
fig = plt.figure(figsize=(15, 6))
gs = gridspec.GridSpec(1, 3)
axes0 = fig.add_subplot(gs[0, 0], projection=m_cmos2)
axes1 = fig.add_subplot(gs[0, 1], projection=m_cdte1)

m_cmos2.plot(axes=axes0, cmap="Greys", vmin=1e0)
m_cmos2_shaped.plot(axes=axes0, vmin=1e0)
shape_plot(shape, m_cmos2_shaped.coordinate_frame, color="tab:red")
draw_grid(m_cmos2, axes0)
axes0.set_title(f"CMOS2\nbl {bottom_left}, tr {top_right}")

m_cdte1_coaligned.plot(axes=axes1, cmap="Greys", vmin=1e0)
m_cdte1_coaligned_shaped.plot(axes=axes1, vmin=1e0)
shape_plot(shape, m_cdte1_coaligned_shaped.coordinate_frame, color="tab:red")
draw_grid(m_cdte1_coaligned, axes1)
axes1.set_title(f"CdTe1 co-aligned to CMOS2\nShift-xy: {shiftxy}")

# can add a time series too
axes2 = fig.add_subplot(gs[0, 2])
obs_time_range = Time(["2024-04-17T22:14:40", "2024-04-17T22:20:41"])
cdte1_times, cdte1_counts = count_time_profile_array(cdte1_events, time_binning=10<<u.second, start_time=obs_time_range[0], end_time=obs_time_range[1])
cmos2_times, cmos2_counts = count_time_profile_array(cmos2_events, time_bins=cdte1_times)
cdte1_times_cs, cdte1_counts_cs = count_time_profile_array(cdte1_events_shaped, time_bins=cdte1_times)
cmos2_times_cs, cmos2_counts_cs = count_time_profile_array(cmos2_events_shaped, time_bins=cdte1_times)
time_support(format="unix_tai")
cdte1_col = "blue"
axes2.stairs(cdte1_counts, cdte1_times.plot_date, color=cdte1_col , label="CdTe1")
axes2.stairs(cdte1_counts_cs, cdte1_times_cs.plot_date, color=cdte1_col, ls=":")
cmos2_col = "k"
axes2.stairs(cmos2_counts, cmos2_times.plot_date, color=cmos2_col , label="CMOS2")
axes2.stairs(cmos2_counts_cs, cmos2_times_cs.plot_date, color=cmos2_col, ls=":")
axes2.set_title(f"{cmos2_events.name} time profile (region=\":\", all counts=\"-\")")
plt.xticks(rotation=30, ha="right")
axes2.set_ylabel(f"{cmos2_counts.unit:latex}")
axes2.set_xlabel(f"Start Time {cdte1_times[0]}")
fmt = mdates.DateFormatter("%H:%M:%S")
axes2.xaxis.set_major_formatter(fmt)
axes2.set_yscale("log")
axes2.set_xlim([cdte1_times.plot_date[0], cdte1_times.plot_date[-1]])
plt.legend()

plt.tight_layout()
plt.show()



# %%
# 5. Filtering by more shapes
# ---------------------------

# start with pre-defined shapes
shape0 = foxsi4_flux_rope_roi()
shape1 = foxsi4_main_loop_roi() 
shape2 = foxsi4_south_east_loop_roi()
shape3 = foxsi4_bonus_microflare_roi()

cmos2_events_shaped0 = by_general_region(cmos2_events, shape0)
cmos2_events_shaped1 = by_general_region(cmos2_events, shape1)
cmos2_events_shaped2 = by_general_region(cmos2_events, shape2)
cmos2_events_shaped3 = by_general_region(cmos2_events, shape3)

m_cmos2_shaped0 = disk_map(cmos2_events_shaped0, bottom_left, top_right)
m_cmos2_shaped1 = disk_map(cmos2_events_shaped1, bottom_left, top_right)
m_cmos2_shaped2 = disk_map(cmos2_events_shaped2, bottom_left, top_right)
m_cmos2_shaped3 = disk_map(cmos2_events_shaped3, bottom_left, top_right)

# we can also add some other random shapes 
# we could define a random polygon shape
vertices = [[-630, -270] << u.arcsec, [-440, -195] << u.arcsec, [-640, -20] << u.arcsec]
shape4 = general_shape(vertices)
# or a customisable ellipse
shape5 = ellipse_shape(center=[-365, -50]<<u.arcsec, semi_axis=[90, 35]<<u.arcsec, rotation=-30<<u.deg)
# there is also `circle_shape` and `rectangle_shape`
# these shapes can be combines (e.g., see the `foxsi4_flux_rope_roi` function source code)
cmos2_events_shaped4 = by_general_region(cmos2_events, shape4)
cmos2_events_shaped5 = by_general_region(cmos2_events, shape5)
m_cmos2_shaped4 = disk_map(cmos2_events_shaped4, bottom_left, top_right)
m_cmos2_shaped5 = disk_map(cmos2_events_shaped5, bottom_left, top_right)

fig = plt.figure(figsize=(10, 6))
gs = gridspec.GridSpec(1, 2)
axes0 = fig.add_subplot(gs[0, 0], projection=m_cmos2)
axes1 = fig.add_subplot(gs[0, 1], projection=m_cmos2)

m_cmos2.plot(axes=axes0, cmap="Greys", vmin=1e0)
m_cmos2_shaped0.plot(axes=axes0, vmin=1e0)
m_cmos2_shaped1.plot(axes=axes0, vmin=1e0)
m_cmos2_shaped2.plot(axes=axes0, vmin=1e0)
m_cmos2_shaped3.plot(axes=axes0, vmin=1e0)
shape_plot(shape0, m_cmos2_shaped.coordinate_frame, color="tab:red")
shape_plot(shape1, m_cmos2_shaped.coordinate_frame, color="purple")
shape_plot(shape2, m_cmos2_shaped.coordinate_frame, color="green")
shape_plot(shape3, m_cmos2_shaped.coordinate_frame, color="orange")
draw_grid(m_cmos2, axes0)
axes0.set_title(f"CMOS2 with\npre-defined regions from\n`foxsi4_science_tools_py.regions_of_interest`")

m_cmos2.plot(axes=axes1, cmap="Greys", vmin=1e0)
m_cmos2_shaped4.plot(axes=axes1, vmin=1e0)
m_cmos2_shaped5.plot(axes=axes1, vmin=1e0)
shape_plot(shape4, m_cmos2_shaped.coordinate_frame, color="black")
shape_plot(shape5, m_cmos2_shaped.coordinate_frame, color="pink")
draw_grid(m_cmos2, axes1)
axes1.set_title(f"CMOS2 with\ncustom shapes")

plt.tight_layout()
plt.show()