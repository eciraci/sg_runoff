"""
Enrico Ciraci 03/2022
Set of utility functions that can be used to generate figures with matplotlib.
"""
# - python dependencies
from __future__ import print_function
import os
import numpy as np
import geopandas as gpd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.ticker as mticker
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib_scalebar.scalebar import ScaleBar
import cartopy.crs as ccrs
from cartopy.io.shapereader import Reader
from cartopy.feature import ShapelyFeature
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from utils.utility_functions_rasterio import load_raster


def add_colorbar(fig: plt.figure, ax: plt.Axes,
                 im: plt.imshow) -> plt.colorbar:
    """
    Add colorbar to the selected plt.Axes.
    :param fig: plt.figure object
    :param ax: plt.Axes object.
    :param im: plt.pcolormesh object.
    :return: plt.colorbar
    """
    divider = make_axes_locatable(ax)
    cax = divider.new_vertical(size='5%', pad=0.6, pack_start=True)
    fig.add_axes(cax)
    cb = fig.colorbar(im, cax=cax, orientation='horizontal')
    return cb
