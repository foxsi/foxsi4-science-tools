"""Module to host functions used to create image products from CdTe data. """

from matplotlib import transforms
import matplotlib.pyplot as plt
import numpy as np

def channel_bins():
    """ Define the strip and ADC bins. """
    strip_bins = np.arange(257)-0.5
    side_strip_bins = np.arange(129)-0.5
    adc_bins = np.arange(1025)-0.5

    return strip_bins, side_strip_bins, adc_bins

def strip_edges_micrometers():
    """ 
    Function to define the physical sizes of the different pitch 
    widths and what to do as they transition. 

    First edge starts at 0 um.
    
    Pitch widths are 100, 80, 60 um and the spaced between are 90 
    and 70 um (100/2+80/2 and 80/2+60/2, respectively).
    """
    C = np.arange(29)*100 # ignore channel 28 just now

    B = np.arange(20)*80 

    A = np.arange(16)*60
    
    new_b = B + C[-1] + 50 + 40 # add last value from (C) then half a bin in (C) and half one in (B)
    
    new_a = A + new_b[-1] + 40 + 30
    
    right_a = A[:-1] + new_a[-1] + 60
    right_b = B + right_a[-1] + 30 + 40
    right_c = C + right_b[-1] + 40 + 50
    
    edges = np.concatenate((C,new_b,new_a,right_a,right_b,right_c))
    
    return edges

def strip_edges_arcminutes():
    """ 
    Function to define the physical sizes of the different pitch 
    widths and what to do as they transition. 
    
    Returns the edges as arcminutes from centre.
    """
    edges = strip_edges_micrometers()

    cdte_fov = 18.7 # arc-minutes

    frac_dist_centred = edges/np.max(edges)-0.5
    
    return frac_dist_centred * cdte_fov
    
def pixel_areas():
    """ From the pitch widths (in um), get the strip-pixel areas (um^2). """
    strip_width_edges = strip_edges_micrometers()
    return np.diff(strip_width_edges)[:,None]@np.diff(strip_width_edges)[None,:]

def filter_max_adc(event_dataframe):
    """
    Function to filer the count data and extracting the strip values with 
    the maximum ADC. 

    The data frame must have keys: 
    - \"index_pt\"
    - \"index_al\"
    - \"adc_cmn_pt\"
    - \"adc_cmn_al\"
    where the indices are the original (NOT re-mapped) values.

    Paramters
    ---------
    event_dataframe : numpy structured array
        The numpy structured array returned from the parser containing 
        the CdTe data.

    Returns
    -------
    `tuple[list]`:
        First list containing the Pt strips corresponding to the maximum 
        ADC value in the allowed strips then the same of the Al side.
        Third list gives the indices of the events that survive selection 
        and were used to select the events in the first two lists.
    """
    
    pt_ind = event_dataframe['index_pt']
    al_ind = event_dataframe['index_al']+128
    pt_adc = event_dataframe['adc_cmn_pt']
    al_adc = event_dataframe['adc_cmn_al']

    # filter Pt side
    pt_strip_selection = ((pt_ind<59) | (pt_ind>68))
    pt_adc[~pt_strip_selection] = 0
    pt_adc_max = np.max(pt_adc, axis=1)
    pt_adc_selection = (pt_adc==pt_adc_max[:,None])
    pt_single_trigs = (np.sum(pt_adc_selection, axis=1)==1) & (pt_adc_max>0)

    # filter Al side
    al_strip_selection = ((al_ind>131) & (al_ind<252))
    al_adc[~al_strip_selection] = 0
    al_adc_max = np.max(al_adc, axis=1)
    al_adc_selection = (al_adc==al_adc_max[:,None])
    al_single_trigs = (np.sum(al_adc_selection, axis=1)==1) & (al_adc_max>0)

    # get single events common to both detector sides
    single_events = (al_single_trigs & pt_single_trigs)
    good_evts = np.nonzero(single_events)[0] 

    return (pt_ind[single_events][pt_adc_selection[single_events]], 
            al_ind[single_events][al_adc_selection[single_events]]-128,
            good_evts)

def get_new_limits_post_rotation(axes, affine_transform, **kwargs):
    """ 
    So `imshow()` does not rescale the limits of the plot after performing 
    an affine transformation to it so need to find the new plot limits
    after a rotation has been performed.

    Note: `pcolormesh` rescales its limts by default.
    """
    if affine_transform is None:
         return

    new_data_corners = get_new_corners(axes, affine_transform, **kwargs)

    axes.set_xlim([np.min(new_data_corners[:,0]), 
                    np.max(new_data_corners[:,0])])
    axes.set_ylim([np.min(new_data_corners[:,1]), 
                    np.max(new_data_corners[:,1])])
        
def get_new_corners(axes, affine_transform, **kwargs):
    """ After any rotation, where are the image corners. """
    x1, x2, y1, y2 = get_image_extent(**kwargs)

    return points_post_rotation(axes, [(x1,y1), (x2,y1), (x1,y2), (x2,y2)], affine_transform)

def get_image_extent(imshow=None, pcolormesh=None):
    """ Get the extent of the image. """
    if pcolormesh is not None:
        x1, x2 = pcolormesh["x_bins"][0], pcolormesh["x_bins"][-1]
        y1, y2 = pcolormesh["y_bins"][0], pcolormesh["y_bins"][-1]
    elif imshow is not None:
        # extent is for non-rotated array
        x1, x2, y1, y2 = imshow.get_extent() 
    return x1, x2, y1, y2

def points_post_rotation(axes, data_points, affine_transform):
    """ Given data-points before any rotation, get their new coordinates. """

    # maybe I will return and see if there is a more direct way, but for now...
    # get new display coordniates of the image corners after the transform
    new_display_corners = affine_transform.transform(data_points)
    # can now convert the display coords to data coords
    return np.array(axes.transData.inverted().transform(new_display_corners))

def image_array(pt_strips, al_strips):
    """Function to make 2D histogram from strip indices. 
    
    This function will remap the CdTe strips to their physical location. 
    """
    _, strip_range, _ = channel_bins()
    im, _, _ = np.histogram2d(pt_strips, 
                              al_strips, 
                              bins=[strip_range, 
                                    strip_range])
    return im[:,::-1] # flip for plotting

def plot_image(image_array, strip_edges=None, rotation=0, reflection=180, figure_kwargs=None, plotting_kwargs=None, axes=None):
    """
    Method to plot the image of the CdTe file.

    Parameters
    ----------
    image_array : `numpy.array`
        A 2d array representing an image.
    
    strip_edges : `numpy.array` or `None`
        The strip edges to be used in the image creation. If given then 
        matplotlib's pcolormesh is used, if `None` then imshow is used. 
        For edge options, some examples from this module are:
        - strip_edges_arcminutes()
        - strip_edges_micrometers()
        The benefit of providing edges is that the true pixel size will 
        be plotted instead of imshow's uniform pixel sizes.
        Default: None
    
    rotation : `int` or `float`
        The rotation of the image. Positive is anti-clockwise and 
        negative is clockwise. 
    
    reflection : `int` or `float`
        An additional 180 rotation/reflection is  applied due to 
        the optics. The affine transform rotation applied in this 
        function is then actually `rotation + reflection`.
        Default: 180 
    
    figure_kwargs : `dict`
        Inputs expanded and passed to the Figure object.
        Default: None
    
    plotting_kwargs: `dict`
        Inputs expanded and passed to either pcolormesh or imshow.
        Default: None 
    
    axes : `matplotlib.axes._axes.Axes`
        The axes to be plotted on. If `None` then `plt.gca()` is used.
        Default: None

    Returns
    -------
    `matplotlib.collections.QuadMesh` or `matplotlib.pyplot.imshow`:
        The object which holds the plot depending if 
        `true_pixel_size=True` or `False`, respectivley.
    """
    figure_kwargs = {} if figure_kwargs is None else figure_kwargs
    plotting_kwargs = {} if plotting_kwargs is None else plotting_kwargs

    axes = plt.gca() if axes is None else axes

    tr = transforms.Affine2D().rotate_deg(rotation+reflection) #rotation_in_degrees
    _plotting_kwargs = {"origin":"lower", "interpolation":"nearest", "rasterized":True, "transform":tr+axes.transData} | plotting_kwargs
    
    if strip_edges is not None:
        _plotting_kwargs.pop("origin", None)
        _plotting_kwargs.pop("interpolation", None)
        i = axes.pcolormesh(strip_edges, 
                            strip_edges, 
                            image_array, 
                            **_plotting_kwargs)
        get_new_limits_post_rotation(axes, 
                                     _plotting_kwargs["transform"], 
                                     pcolormesh={"x_bins":strip_edges, 
                                                 "y_bins":strip_edges})
    else:
        i = axes.imshow(image_array, 
                        **_plotting_kwargs)
        get_new_limits_post_rotation(axes,
                                     _plotting_kwargs["transform"], 
                                     imshow=i)
    
    plt.xlabel("Al")
    plt.ylabel("Pt")

    axes.set_aspect("equal")
    
    return axes, i