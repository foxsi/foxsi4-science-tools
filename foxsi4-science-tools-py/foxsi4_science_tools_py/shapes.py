import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np
import shapely
from astropy.coordinates import SkyCoord
from shapely.geometry import Point
from shapely.geometry.multipolygon import MultiPolygon
from shapely.geometry.polygon import Polygon


__all__ = [
    "circle_shape",
    "rectangle_shape",
    "general_shape",
    "shape_plot",
    "points_plot",
]

def general_shape(*args):
    """A function to turn a list of (x, y) points into a polygon.

    The coordinates will be the vertices of the desired shape and passed
    to ``shapely.geometry.polygon.Polygon``.

    Parameters
    ----------
    *args :
        The vertices of the polygon and should be in or convertible to
        ``~atropy.units.Quantity``.
    """
    # make sure arcsecond
    # need to index [0] as it weirdly becomes a list in a list for some reason
    unit_check_list = next(a << u.arcsec for a in args)
    unitless_for_poly = [ucl.value for ucl in unit_check_list]
    return Polygon(unitless_for_poly)


def rectangle_shape(bottom_left, top_right):
    """Return a ``shapely.geometry.polygon.Polygon`` rectangle shape.

    Parameters
    ----------
    bottom_left : `~astropy.units.Quantity`
        The bottom left corner (X, Y) of the desired box. Should be
        convertible to arcseconds.

    top_right : `~astropy.units.Quantity`
        The top right corner (X, Y) of the desired box. Should be
        convertible to arcseconds.
    """
    vertices = [
        bottom_left << u.arcsec,
        [bottom_left[0], top_right[1]] << u.arcsec,
        top_right << u.arcsec,
        [top_right[0], bottom_left[1]] << u.arcsec,
    ]
    return general_shape(vertices)


def circle_shape(center, radius):
    """Return a ``shapely.geometry.polygon.Polygon`` circle shape.

    Parameters
    ----------
    center : `~astropy.units.Quantity`
        Thecenter (X, Y) of the desired circle. Should be convertible to
        arcseconds.

    radius : `~astropy.units.Quantity`
        The radius R of the desired circle. Should be convertible to
        arcseconds.
    """
    _center = center << u.arcsec
    _radius = radius << u.arcsec
    return Point(*_center.value).buffer(_radius.value)


def ellipse_shape(center, semi_axis, rotation):
    """Return a ``shapely.geometry.polygon.Polygon`` ellipse shape.

    See https://gis.stackexchange.com/questions/243459/drawing-ellipse-with-shapely
    for solution.

    Parameters
    ----------
    center : `~astropy.units.Quantity`
        Thecenter (X, Y) of the desired circle. Should be convertible to
        arcseconds.

    semi_axis : `~astropy.units.Quantity`
        The semi-axis along the x-axis and the y-axis. Should be convertible to 
        arcseconds.

    rotation : `~astropy.units.Quantity`
        Rotation angle, clockwise is -ve and anti-clockwise is +ve. Should be 
        convertible to degrees.
    """
    _center = center << u.arcsec
    _semi_axis = semi_axis << u.arcsec
    _rotation = rotation << u.deg

    # create a uit circle
    circle = circle_shape(_center, 1)

    # now stretch that circle
    ellipse  = shapely.affinity.scale(circle, *_semi_axis.value)

    # now rotate and return
    return shapely.affinity.rotate(ellipse, _rotation.value)


def shape_plot(
    shapely_polygon, coordinate_frame, axes=None, sky_coord_kwargs=None, **kwargs
):
    """Plot a ``shapely.geometry.polygon.Polygon``.

    Parameters
    ----------
    shapely_polygon :  `shapely.geometry.polygon.Polygon` | `geometry.multipolygon.MultiPolygon`
        Will assume ``shapely_polygon.exterior.coords`` is in arcseconds
        by default (see `sky_coord_kwargs` for more). If ``geometry.multipolygon.MultiPolygon``
        (say a symmetric different has been used, i.e. multiple separate
        polygons exist) then these will be plotted with colours cycling
        through Matplotib's color cycle. Although, passing a line
        property such as ``c`` or ``color`` as a ``kwargs`` will apply
        to all polygons.

    coordinate_frame : Sunpy coordinate map frame

    axes : axis object
        The axis to be plotted on. If None then ``plt.gca()`` is used.
        Default: None

    sky_coord_kwargs : `dict` | `None
        Passed to ``astropy.coordinates.SkyCoord`` If None then the
        default unit for the ``shapely.geometry.polygon.Polygon`` will
        be arcseconds.
        Default: None

    **kwargs :
        Passed to ``plot_coord``.

    Returns
    -------
    Axes input.
    """
    axes = plt.gca() if axes is None else axes
    sky_coord_kwargs = {} if sky_coord_kwargs is None else sky_coord_kwargs
    _sky_coord_kwargs = {"unit": u.arcsec} | sky_coord_kwargs

    if isinstance(shapely_polygon, Polygon):
        vertices = np.array(shapely_polygon.exterior.coords)
        points_plot(
            vertices[:, 0],
            vertices[:, 1],
            coordinate_frame,
            axes=axes,
            sky_coord_kwargs=_sky_coord_kwargs,
            **kwargs,
        )

    elif isinstance(shapely_polygon, MultiPolygon):
        for shape in shapely_polygon.geoms:
            vertices = np.array(shape.exterior.coords)
            points_plot(
                vertices[:, 0],
                vertices[:, 1],
                coordinate_frame,
                axes=axes,
                sky_coord_kwargs=_sky_coord_kwargs,
                **kwargs,
            )
    else:
        raise TypeError(
            "Input for ``shapely_polygon`` is not of type ``geometry.multipolygon.MultiPolygon`` or ``geometry.polygon.Polygon``. "
        )
    return axes


def points_plot(x, y, coordinate_frame, axes=None, sky_coord_kwargs=None, **kwargs):
    """Plot a set of coordinates.

    Parameters
    ----------
    x, y :  `list` | array-like
        The x and y coordinates to be plotted. Given in solar pixel
        coordinates and plotted in helioprojective.

    coordinate_frame : Sunpy coordinate map frame

    axes : axis object
        The axis to be plotted on. If None then ``plt.gca()`` is used.
        Default: None

    sky_coord_kwargs : `dict` | `None
        Passed to ``astropy.coordinates.SkyCoord`` If None then the
        default unit for the ``shapely.geometry.polygon.Polygon`` will
        be arcseconds.
        Default: None

    **kwargs :
        Passed to ``plot_coord``.

    Returns
    -------
    Axes input.
    """
    axes = plt.gca() if axes is None else axes
    sky_coord_kwargs = {} if sky_coord_kwargs is None else sky_coord_kwargs
    _sky_coord_kwargs = {"unit": u.arcsec} | sky_coord_kwargs
    sc = SkyCoord(x, y, frame=coordinate_frame, **_sky_coord_kwargs)
    axes.plot_coord(sc, **kwargs)
    return axes