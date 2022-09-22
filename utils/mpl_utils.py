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


def plot_dhdt_map(data_dir: str, img: np.array, xx: np.array, yy: np.array,
                  out_path: str, ref_crs: int, extent: int = 1,
                  land_color: str = "black", ice_color: str = "grey",
                  grnd_ln_color: str = "k", grnd_zn_color: str = "g",
                  grnd_zn_buffer: float = 0, title: str = "",
                  annotate: str = "", cmap=plt.get_cmap("bwr_r"),
                  vmin: int = -10, vmax: int = 10,
                  fig_format: str = "jpeg") -> None:
    """
    Plot Elevation Change DhDt Map - Wide Domain
    :param data_dir: Project Data Directory
    :param img: elevation change [numpy array]
    :param xx: xx m-grid - centroids [numpy array]
    :param yy: yy m-grid - centroids [numpy array]
    :param out_path: absolute path to output file
    :param ref_crs: coordinates reference system
    :param extent: figure extent [1,2,3]
    :param land_color: land edges color
    :param ice_color: ice edges color
    :param grnd_ln_color: Grounding Line Color
    :param grnd_zn_color: Grounding Zone Color
    :param grnd_zn_buffer: Grounding Zone Buffer
    :param title: figure title
    :param annotate: add annotate object
    :param fig_format: figure format [jpeg]
    :param cmap: imshow/pcolormesh color map
    :param vmin: imshow/pcolormesh vmin
    :param vmax: imshow/pcolormesh vmax
    :return: None
    """
    # - text size
    txt_size = 14  # - main text size
    leg_size = 13  # - legend text size
    label_size = 12  # - label text size

    if extent == 1:
        # - Map Extent 1 - wide
        map_extent = [-68, -55, 80, 82]
        figsize = (9, 9)
    elif extent == 2:
        # - Map Extent 2 - zoom 1
        map_extent = [-61.5, -57, 79.5, 81.5]
        figsize = (5, 8)
    elif extent == 3:
        # - Map Extent 3 - zoom2
        map_extent = [-61.5, -57.6, 80.2, 81.2]
        figsize = (9, 9)
    else:
        # - Map Extent 4 - zoom3
        map_extent = [-61.1, -59.9, 80.4, 81.2]
        figsize = (5, 8)
        txt_size = 12  # - main text size
        leg_size = 10   # - legend text size
        label_size = 10  # - label text size

    # - Path to Ice and Land Masks
    ics_shp = os.path.join("..", "esri_shp", "GIMP",
                           "Petermann_Domain_glaciers_wgs84.shp")
    land_shp = os.path.join("..", "esri_shp", "GIMP",
                            "GSHHS_i_L1_Petermann_clip.shp")

    # - Petermann Grounding Line - 2021/2021
    gnd_ln_shp = os.path.join(data_dir, "coco_petermann_grnd_lines_2020-2021",
                              "grnd_lines_shp_to_share",
                              "coco20200501_20200502-20200517_20200518",
                              "coco20200501_20200502-20200517_20200518"
                              "_grnd_line.shp")
    gnd_ln_df = gpd.read_file(gnd_ln_shp).to_crs(epsg=ref_crs)
    # -
    xg, yg = gnd_ln_df["geometry"].geometry[0].xy

    # - Petermann Grounding Zone - 2011/2021
    gnd_zn_shp = os.path.join(data_dir, "GIS_Data",
                              "Petermann_features_extraction",
                              "Petermann_grounding_line_migration_"
                              "range_epsg3413.shp")
    if grnd_zn_buffer:
        # - Clip Output Raster
        clip_shp_mask_path \
            = os.path.join(data_dir, "GIS_Data",
                           "Petermann_features_extraction",
                           "Petermann_iceshelf_clipped_epsg3413.shp")
        clip_mask = gpd.read_file(clip_shp_mask_path).to_crs(epsg=ref_crs)

        gnd_zn_shp_buff = os.path.join(data_dir, "GIS_Data",
                                       "Petermann_features_extraction",
                                       "Petermann_grounding_line_migration_"
                                       f"range_buff{grnd_zn_buffer}"
                                       f"_epsg3413.shp")
        if not os.path.isfile(gnd_zn_shp_buff):
            gnd_zn_to_bf = gpd.read_file(gnd_zn_shp).to_crs(epsg=ref_crs)
            gnd_zn_to_bf["geometry"] = gnd_zn_to_bf.geometry\
                .buffer(grnd_zn_buffer)
            # - clip the obtained buffered mask with the
            # - ice shelf perimeter mask.
            gnd_zn_to_bf = gpd.overlay(gnd_zn_to_bf, clip_mask,
                                       how="intersection")
            # - save buffered mask to file
            gnd_zn_to_bf.to_file(gnd_zn_shp_buff)

        # -
        gnd_zn_shp = gnd_zn_shp_buff

    # - set Coordinate Reference System
    ref_crs = ccrs.NorthPolarStereo(central_longitude=-45,
                                    true_scale_latitude=70)
    # - initialize matplotlib figure object
    fig = plt.figure(figsize=figsize)
    # - initialize legend labels
    leg_label_list = []
    ax = fig.add_subplot(1, 1, 1, projection=ref_crs)
    ax.set_extent(map_extent, crs=ccrs.PlateCarree())
    # - Plot Coastlines
    shape_feature = ShapelyFeature(Reader(land_shp).geometries(),
                                   ccrs.PlateCarree())
    ax.add_feature(shape_feature, facecolor="None", edgecolor=land_color)
    # - Plot Glaciers Mask
    shape_feature = ShapelyFeature(Reader(ics_shp).geometries(),
                                   ccrs.PlateCarree())
    ax.add_feature(shape_feature, facecolor="None", edgecolor=ice_color)

    # - Plot Grounding Line 2020/2021
    l1, = ax.plot(xg, yg, color=grnd_ln_color, lw=2, zorder=10, ls="-.")
    leg_label_list.append("Grounding Line 2020")
    # - Plot Grounding Zone 2011/2021
    shape_feature = ShapelyFeature(Reader(gnd_zn_shp).geometries(), ref_crs)
    ax.add_feature(shape_feature, facecolor="None",
                   edgecolor=grnd_zn_color, linestyle="--",
                   linewidth=2)
    leg_label_list.append("Grounding Zone 2011-2020")
    l2 = mpatches.Rectangle((0, 0), 1, 0.1, linewidth=2,
                            edgecolor=grnd_zn_color, facecolor="none",
                            linestyle="--")

    # - Set Map Grid
    gl = ax.gridlines(draw_labels=True, dms=True, x_inline=False,
                      y_inline=False)
    gl.top_labels = False
    gl.bottom_labels = True
    gl.right_labels = False
    gl.xlocator = mticker.FixedLocator(np.arange(map_extent[0] - 5,
                                                 map_extent[1] + 5, 2))
    gl.ylocator = mticker.FixedLocator(
        np.arange(np.floor(map_extent[2]) - 5,
                  np.floor(map_extent[3]) + 5,
                  0.5))
    gl.xlabel_style = {"rotation": 0, "weight": "bold", "size": label_size}
    gl.ylabel_style = {"rotation": 0, "weight": "bold", "size": label_size}
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER

    # - Figure title
    ax.set_title(title, weight="bold", loc="left", size=txt_size)
    # - Add Figure Annotation
    ax.annotate(annotate, xy=(0.03, 0.03), xycoords="axes fraction",
                size=label_size, zorder=100,
                bbox=dict(boxstyle="square", fc="w", alpha=0.8))

    # - Plot DH/DT map
    im = ax.pcolormesh(xx, yy, img, cmap=cmap,
                       zorder=0, vmin=vmin, vmax=vmax)

    # - colorbar - solution specific for cartopy figures
    divider = make_axes_locatable(ax)
    ax_cb = divider.new_horizontal(size="5%", pad=0.1, axes_class=plt.Axes)
    fig.add_axes(ax_cb)
    cb = plt.colorbar(im, cax=ax_cb)
    cb.set_label(label="[m/year]", weight="bold", size=label_size)
    cb.ax.tick_params(labelsize=leg_size)

    # - Add Legend to DhDt Maps
    ax.legend([l1, l2], leg_label_list, loc="upper right",
              fontsize=leg_size, framealpha=0.8,
              facecolor="w", edgecolor="k")

    # - Add ScaleBar
    ax.add_artist(ScaleBar(1, units="m", location="lower right",
                           border_pad=1, pad=0.5, box_color="w",
                           frameon=True))

    # - save output figure
    plt.savefig(out_path, dpi=200, format=fig_format)
    plt.close()


def plot_dhdt_map_zoom(data_dir: str, dhdt_path: str, out_path: str,
                       ref_crs: int, extent: int = 1, land_color: str = "black",
                       ice_color: str = "grey", grnd_ln_color: str = "k",
                       grnd_zn_color: str = "g", grnd_zn_buffer: float = 0,
                       title: str = "", annotate: str = "",
                       fig_format: str = "jpeg",
                       cmap=plt.get_cmap("bwr_r"),
                       vmin: int = -10, vmax: int = 10
                       ) -> None:
    """
    Plot Ice Elevation Change DhDt Map - Zoom Domain
    :param data_dir: absolute path to project data directory
    :param dhdt_path: absolute path to dhdt in GeoTiff format
    :param out_path: absolute path to output figure
    :param ref_crs: reference CRS
    :param extent: map extent
    :param land_color: land edges color
    :param ice_color: ice edges color
    :param grnd_ln_color: Grounding Line Color
    :param grnd_zn_color: Grounding Zone Color
    :param grnd_zn_buffer: Grounding Zone Buffer
    :param title: figure title
    :param annotate: add annotate object
    :param fig_format: figure format [jpeg]
    :param cmap: imshow/pcolormesh color map
    :param vmin: imshow/pcolormesh vmin
    :param vmax: imshow/pcolormesh vmax
    :return:
    """
    # - text size
    txt_size = 14       # - main text size
    leg_size = 13       # - legend text size
    label_size = 12     # - label text siz

    # - Plot DhDt Map
    dhdt_plot = load_raster(dhdt_path)
    dhdt_fig = dhdt_plot["data"]
    dhdt_fig[dhdt_fig == dhdt_plot["nodata"]] = np.nan
    x_coords = dhdt_plot["x_centroids"]
    y_coords = dhdt_plot["y_centroids"]
    xx, yy = np.meshgrid(x_coords, y_coords)

    # - Map Extent
    if extent == 1:
        map_extent = [-61.1, -59.9, 80.4, 81.2]
        figsize = (6, 9)
    else:
        map_extent = [-60.8, -59.1, 80.4, 80.7]
        figsize = (9, 9)

    # - Path to Ice and Land Masks
    ics_shp = os.path.join("..", "esri_shp", "GIMP",
                           "Petermann_Domain_glaciers_wgs84.shp")
    land_shp = os.path.join("..", "esri_shp", "GIMP",
                            "GSHHS_i_L1_Petermann_clip.shp")
    # - Petermann Grounding Line - 2021/2021
    gnd_path = grounding_line_path(data_dir)
    gnd_ln_shp = gnd_path["gnd_ln_shp"]
    gnd_ln_df = gpd.read_file(gnd_ln_shp).to_crs(epsg=ref_crs)
    # -
    xg, yg = gnd_ln_df["geometry"].geometry[0].xy

    # - Petermann Grounding Zone - 2011/2021
    gnd_zn_shp = gnd_path["gnd_zn_shp"]

    if grnd_zn_buffer:
        # - Clip Output Raster
        clip_shp_mask_path \
            = os.path.join(data_dir, "GIS_Data",
                           "Petermann_features_extraction",
                           "Petermann_iceshelf_clipped_epsg3413.shp")
        clip_mask = gpd.read_file(clip_shp_mask_path).to_crs(epsg=ref_crs)

        gnd_zn_shp_buff = os.path.join(data_dir, "GIS_Data",
                                       "Petermann_features_extraction",
                                       "Petermann_grounding_line_migration_"
                                       f"range_buff{grnd_zn_buffer}"
                                       f"_epsg3413.shp")
        if not os.path.isfile(gnd_zn_shp_buff):
            gnd_zn_to_bf = gpd.read_file(gnd_zn_shp).to_crs(epsg=ref_crs)
            gnd_zn_to_bf["geometry"] = gnd_zn_to_bf.geometry\
                .buffer(grnd_zn_buffer)
            # - clip the obtained buffered mask with the
            # - ice shelf perimeter mask.
            gnd_zn_to_bf = gpd.overlay(gnd_zn_to_bf, clip_mask,
                                       how="intersection")
            # - save buffered mask to file
            gnd_zn_to_bf.to_file(gnd_zn_shp_buff)

        # -
        gnd_zn_shp = gnd_zn_shp_buff

    # - set Coordinate Reference System
    ref_crs = ccrs.NorthPolarStereo(central_longitude=-45,
                                    true_scale_latitude=70)
    fig = plt.figure(figsize=figsize)
    # - initialize legend labels
    leg_label_list = []

    # - Plot DhDt Map
    ax = fig.add_subplot(1, 1, 1, projection=ref_crs)
    ax.set_extent(map_extent, crs=ccrs.PlateCarree())
    # - Plot Coastlines
    shape_feature = ShapelyFeature(Reader(land_shp).geometries(),
                                   ccrs.PlateCarree())
    ax.add_feature(shape_feature, facecolor="None", edgecolor=land_color)

    # - Plot Glaciers Mask
    shape_feature = ShapelyFeature(Reader(ics_shp).geometries(),
                                   ccrs.PlateCarree())
    ax.add_feature(shape_feature, facecolor="None", edgecolor=ice_color)

    # - Plot Grounding Line 2020/2021
    l1, = ax.plot(xg, yg, color=grnd_ln_color, lw=2, zorder=10, ls="-.")
    leg_label_list.append("Grounding Line 2020")
    # - Plot Grounding Zone 2011/2021
    shape_feature = ShapelyFeature(Reader(gnd_zn_shp).geometries(), ref_crs)
    ax.add_feature(shape_feature, facecolor="None",
                   edgecolor=grnd_zn_color, linestyle="--",
                   linewidth=2)
    leg_label_list.append("Grounding Zone 2011-2020")
    l2 = mpatches.Rectangle((0, 0), 1, 0.1, linewidth=2,
                            edgecolor=grnd_zn_color, facecolor="none",
                            linestyle="--")

    # - Set Map Grid
    gl = ax.gridlines(draw_labels=True, dms=True, x_inline=False,
                      y_inline=False, color="k", linestyle="dotted",
                      alpha=0.3)
    gl.top_labels = False
    gl.bottom_labels = True
    gl.right_labels = False
    gl.xlocator \
        = mticker.FixedLocator(np.arange(np.floor(map_extent[0]) - 3.5,
                                         np.floor(map_extent[1]) + 3, 1))
    gl.ylocator \
        = mticker.FixedLocator(np.arange(np.floor(map_extent[2]) - 5,
                                         np.floor(map_extent[3]) + 5, 0.2))
    gl.xlabel_style = {"rotation": 0, "weight": "bold", "size": label_size}
    gl.ylabel_style = {"rotation": 0, "weight": "bold", "size": label_size}
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER

    # - Plot DH/DT map
    im = ax.pcolormesh(xx, yy, dhdt_fig, cmap=cmap,
                       zorder=0, vmin=vmin, vmax=vmax, rasterized=True)

    # add an axes above the main axes.
    divider = make_axes_locatable(ax)
    ax_cb = divider.new_horizontal(size="5%", pad=0.1, axes_class=plt.Axes)
    fig.add_axes(ax_cb)
    cb = plt.colorbar(im, cax=ax_cb)
    cb.set_label(label="[m/year]", weight="bold", size=label_size)
    cb.ax.tick_params(labelsize="medium")

    # - Figure title
    ax.set_title(title, weight="bold", loc="left", size=txt_size)
    # - Add Figure Annotation
    ax.annotate(annotate, xy=(0.03, 0.03), xycoords="axes fraction",
                size=label_size, zorder=100,
                bbox=dict(boxstyle="square", fc="w", alpha=0.8))

    # - Add Legend to DhDt Maps
    ax.legend([l1, l2], leg_label_list, loc="upper right",
              fontsize=leg_size, framealpha=0.8,
              facecolor="w", edgecolor="k")

    # - Add ScaleBar
    ax.add_artist(ScaleBar(1, units="m", location="lower right",
                           border_pad=1, pad=0.5, box_color="w",
                           frameon=True))

    # - save output figure
    plt.savefig(out_path, dpi=200, format=fig_format)
    plt.close()