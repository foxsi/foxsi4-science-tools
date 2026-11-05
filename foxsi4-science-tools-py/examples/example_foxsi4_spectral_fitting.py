"""

"""

import astropy.units as u
import matplotlib.pyplot as plt
from matplotlib import gridspec

from foxsi4_science_tools_py.io.detectors import load_cdte, load_cmos
from foxsi4_science_tools_py.list_filters import by_general_region
from foxsi4_science_tools_py.util.clean import cdte_usual_suspects, cmos_usual_suspects
from foxsi4_science_tools_py.util.maps import draw_grid
from foxsi4_science_tools_py.foxsi_core import disk_map, full_disk_map, coalign_event_lists, event_data_xy_shift
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
cmos2_file = "/Users/kris/Documents/umnPostdoc/projects/both/foxsi4/foxsi4-analysis-tools/cmos-tools-main/cmos-tools/FOXSI-4/event_process/EVENT_FILES/level2/foxsi4_telescope-1_CMOS-2_PHOTON_EVENT_LEVEL2_V26SEP06.fits"

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

plt.tight_layout()
plt.show()
