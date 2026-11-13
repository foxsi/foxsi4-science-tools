import astropy.units as u

from foxsi4_science_tools_py.shapes import ellipse_shape

__all__ = ["foxsi4_main_loop_roi",
           "foxsi4_flux_rope_roi",
           "foxsi4_south_east_loop_roi",
           "foxsi4_bonus_microflare_roi",
           ]


def foxsi4_main_loop_roi():
    """A Shapely polygon containing the FOXSI-4 core flare region."""
    return ellipse_shape(center=[-365, -60]<<u.arcsec, semi_axis=[90, 35]<<u.arcsec, rotation=-16<<u.deg)

def foxsi4_flux_rope_roi():
    """A Shapely polygon containing the FOXSI-4 flux rope region."""
    return ellipse_shape(center=[-350, -30]<<u.arcsec, semi_axis=[77, 50]<<u.arcsec, rotation=-20<<u.deg).difference(foxsi4_main_loop_roi())

def foxsi4_south_east_loop_roi():
    """A Shapely polygon ellipse area for the south-east loop region."""
    return ellipse_shape(center=[-415, -118]<<u.arcsec, semi_axis=[57, 30]<<u.arcsec, rotation=-40<<u.deg)

def foxsi4_bonus_microflare_roi():
    """A Shapely polygon ellipse area for the bonus microflare region."""
    return ellipse_shape(center=[-270, -210]<<u.arcsec, semi_axis=[43, 25]<<u.arcsec, rotation=-25<<u.deg)
