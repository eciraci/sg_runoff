#!/usr/bin/env python
"""
Written by Enrico Ciraci' - (09/2022)

Test Scripts:
Compute sub-glacial discharge flow by employing ice thickness and bed topography
from BedMachine v4.0. Discharge flow is determined by the spatial gradient of the
Hydraulic Potential computed at the interface between ice and bedrock.

For more Info about tha Hydraulic Potential evaluation see:
Cuffey & Peterson, 2006 - "The Physics of Glaciers" - Chapter #6.
https://www.elsevier.com/books/the-physics-of-glaciers/ cuffey/978-0-12-369461-4

Hydraulic Potential:
    phi_x = Pw + rho_water * gravity * Z
        - Where:
        - Pw is Water Pressure calculated as:
        -     Pw = rho_ice * gravity * [S - Z]
        -          units: [kg / m3] * [m / s2] * [m] = [kg / m / s2]
        -
        - Z is Ice Bed Elevation
        - S is Ice Surface Elevation


Once the Hydraulic Potential has been computed, the accumulated flow is
calculated by employing the method D8 described int:

O’Callaghan, J.F., Mark, D.M., 1984. The Extraction of Drainage Networks from
Digital Elevation Data. Computer vision, graphics, and image
processing 28, 323–344.

The D8 method assigns flow from a focal cell to one and only one of its 8
neighbouring cells. The chosen neighbour is the one accessed via the steepest
slope. When such a neighbour does not exist, no flow direction is assigned.
When two or more neighbours have the same slope, the chosen neighbour is the
first one considered by the algorithm.

-> This is a convergent, deterministic flow method.

-> More details here: https://richdem.readthedocs.io/en/latest/flow_metrics.html

NOTE: In this case, the implementation of the D8 algorithm provided by the
      Pysheds project is used -> https://github.com/mdbartos/pysheds
      The current implementation of this routine assigns directions arbitrarily
      if multiple steepest paths exist.

To Run this script digit:

$ python compute_flowacc_pysheds.py


PYTHON DEPENDENCIES:
pysheds: Simple and fast watershed delineation in python.
    https://github.com/mdbartos/pysheds
numpy: The fundamental package for scientific computing with Python.
    https://numpy.org/
matplotlib: Visualization with Python.
    https://matplotlib.org/
fiona: Fiona is GDAL’s neat and nimble vector API for Python programmers.
    https://github.com/Toblerity/Fiona
seaborn: statistical data visualization in Python.
    https://seaborn.pydata.org/
rasterio: access to geospatial raster data in Python.
    https://rasterio.readthedocs.io
"""

import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors
import seaborn as sns
import fiona
from pysheds.grid import Grid
# - Library Dependencies
from utility_functions_rasterio import load_raster, save_raster, clip_raster


def create_dir(abs_path: str, dir_name: str) -> str:
    """
    Create directory
    :param abs_path: absolute path to the output directory
    :param dir_name: new directory name
    :return: absolute path to the new directory
    """
    dir_to_create = os.path.join(abs_path, dir_name)
    if not os.path.exists(dir_to_create):
        os.mkdir(dir_to_create)
    return dir_to_create


def main() -> None:
    # - Processing Parameters
    # - PySheds implements two flow routing algorithms:
    # - d8 : D8 flow directions
    # - dinf : D-infinity flow directions
    routing = 'dinf'
    # - Project Data Directory
    project_dir = os.path.join('/', 'Volumes', 'Extreme Pro', 'BedMachine',
                               'output.dir')
    output_dir = create_dir(project_dir, 'subglacial_flow_pysheds')
    output_dir = create_dir(output_dir, routing)
    # - Processing Parameters
    # - Accumulation value used as threshold for River Network determination.
    rn_thresh = 100

    # - Figure Parameters
    fig_format = 'jpeg'
    dpi = 200

    # - Absolute Path to BedMachine Data
    # - Bedrock
    bed_path = \
        os.path.join('/', 'Volumes', 'Extreme Pro', 'BedMachine',
                     'output.dir', 'BedMachineGreenland-version_05_bed',
                     'Petermann_Domain_Velocity_Stereo_bed'
                     '_EPSG-3413_res150.tiff')
    # - Ice Elevation
    elev_path \
        = os.path.join('/', 'Volumes', 'Extreme Pro', 'BedMachine',
                       'output.dir', 'BedMachineGreenland-version_05_surface',
                       'Petermann_Domain_Velocity_Stereo_surface'
                       '_EPSG-3413_res150.tiff')

    # - Ice Thickness
    thick_path \
        = os.path.join('/', 'Volumes', 'Extreme Pro', 'BedMachine',
                       'output.dir', 'BedMachineGreenland-version_05_thickness',
                       'Petermann_Domain_Velocity_Stereo_thickness'
                       '_EPSG-3413_res150.tiff')

    # - Compute Hydraulic Potential
    print('# - Compute Hydraulic Potential.')
    bed_input = load_raster(bed_path)
    elev_input = load_raster(elev_path)
    thick_input = load_raster(thick_path)
    # - compute ice mask from ice thickness
    ice_mask = np.zeros(np.shape(thick_input['data']))
    ice_mask[thick_input['data'] > 0] = 1

    # - Conversion Parameters
    rho_water = 1000        # Water Density [kg / m3]
    rho_ice = 917           # Ice Density [kg / m3]
    gravity = 9.81          # Gravity Acceleration [m / s2]

    # For more Info see:
    # Cuffey & Peterson, 2006 - "The Physics of Glaciers" - Chapter #6.
    # https://www.elsevier.com/books/the-physics-of-glaciers/
    #         cuffey/978-0-12-369461-4
    #
    # - Hydraulic Potential:
    # phi_x = Pw + rho_water * gravity * Z
    # - Where:
    # - Pw is Water Pressure calculated as:
    # -     Pw = rho_ice * gravity * [S - Z]
    # -          units: [kg / m3] * [m / s2] * [m] = [kg / m / s2]
    # -
    # - Z is Ice Bed Elevation
    # - S is Ice Surface Elevation

    surf_elev = elev_input['data']
    bed_elev = bed_input['data']
    x_vect = elev_input['x_coords']
    y_vect = elev_input['y_coords']
    crs = elev_input['crs']
    res = elev_input['res']

    # - Hydraulic Potential
    hydro_pot = ((rho_ice * gravity * surf_elev)
                 + ((rho_water - rho_ice) * gravity * bed_elev))
    hydro_pot[ice_mask == 0] = np.nan

    # -
    out_path = os.path.join(output_dir, 'Petermann_Domain_Velocity'
                                        '_Stereo_hydro_pot_EPSG-3413'
                                        '_res150.tiff')
    save_raster(hydro_pot, res, x_vect, y_vect, out_path, crs)

    # - Need to clip the Ice Shelf surface before computing the
    # - accumulated flow.
    ice_shelf_mask = os.path.join('/', 'Volumes', 'Extreme Pro', 'GIS_Data',
                                  'Petermann_Ice_Shelf',
                                  'Petermann_Domain-Ice_Shelf_Mask_GL-'
                                  'ERS2011.shp')
    out_path_clip = out_path.replace('_Stereo_hydro_pot',
                                     '_Stereo_hydro_pot_clipped')
    clip_raster(out_path, ice_shelf_mask, out_path_clip, nodata=np.nan)
    # - The clipping reduces the raster dimension by one row

    # - Load Clipped Mask
    hydro_pot_in = load_raster(out_path_clip)
    ice_mask = np.zeros(np.shape(hydro_pot_in['data']))
    ice_mask[np.isfinite(hydro_pot_in['data'])] = 1
    x_vect = hydro_pot_in['x_coords']
    y_vect = hydro_pot_in['y_coords']
    crs = hydro_pot_in['crs']
    res = hydro_pot_in['res']
    hydro_pot_in['data'][ice_mask == 0] = np.nan
    save_raster(hydro_pot_in['data'], res, x_vect, y_vect, out_path_clip, crs)

    print('# - Compute Accumulated Sub-Glacial Flow using Pysheds.')
    # - Read Datasets using pysheds.grid
    grid = Grid.from_raster(out_path_clip, data_name='hydro_pot')
    dem = grid.read_raster(out_path_clip)

    # - Show Hydraulic Potential Map
    fig, ax = plt.subplots(figsize=(8, 6))
    fig.patch.set_alpha(0)
    plt.imshow(dem, extent=grid.extent, cmap='Greys', zorder=1)
    plt.colorbar(label='[kg / m / s2]')
    plt.grid(zorder=0)
    plt.title('Hydraulic Potential Map', size=14)
    plt.xlabel('Longitude')
    plt.ylabel('Latitude')
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir,
                             f'hydraulic_potential.{fig_format}'),
                dpi=dpi, format=fig_format)

    # - Condition DEM -> Hydraulic Potential
    # - Fill pits in DEM -> Hydraulic Potential
    pit_filled_dem = grid.fill_pits(dem)
    # - Fill depressions in DEM -> Hydraulic Potential
    flooded_dem = grid.fill_depressions(pit_filled_dem)
    # - Resolve flats in DEM -> Hydraulic Potential
    inflated_dem = grid.resolve_flats(flooded_dem)
    inflated_dem[np.flipud(ice_mask) == 0] = np.nan

    # - Compute Elevation Flow Direction
    # - Determine D8 flow directions from DEM -> Hydraulic Potential
    dirmap = (64, 128, 1, 2, 4, 8, 16, 32)      # -  Specify directional mapping

    # - Compute flow directions
    if routing == 'd8':
        fdir = grid.flowdir(inflated_dem, dirmap=dirmap, nodata_in=np.nan,
                            routing=routing)
    else:
        fdir = grid.flowdir(inflated_dem, dirmap=dirmap, nodata_in=np.nan,
                            routing=routing, flats=np.nan, pits=np.nan)

    # - Show Flow Direction
    fig = plt.figure(figsize=(8, 6))
    fig.patch.set_alpha(0)
    plt.imshow(fdir, extent=grid.extent, cmap='viridis', zorder=2)
    boundaries = ([0] + sorted(list(dirmap)))
    plt.colorbar(boundaries=boundaries,
                 values=sorted(dirmap))
    plt.xlabel('Longitude')
    plt.ylabel('Latitude')
    plt.title('Flow direction grid', size=14)
    plt.grid(zorder=-1)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir,
                             f'flow_direction_grid_{routing}.{fig_format}'),
                dpi=dpi, format=fig_format)
    plt.close()

    # - Calculate flow accumulation
    acc = grid.accumulation(fdir, dirmap=dirmap, routing=routing,
                            nodata_in=np.nan)
    acc[np.flipud(ice_mask) == 0] = np.nan

    # - Save the obtained raster
    out_path = os.path.join(output_dir, f'flow_accumulation_{routing}.tiff')
    save_raster(np.flipud(acc), res, x_vect, y_vect, out_path, crs)
    out_path = os.path.join(output_dir, f'flow_accumulation_{routing}_log.tiff')
    save_raster(np.log(np.flipud(acc)), res, x_vect, y_vect, out_path, crs)

    # - Show Flow Accumulation
    fig, ax = plt.subplots(figsize=(8, 6))
    fig.patch.set_alpha(0)
    plt.grid('on', zorder=0)
    im = ax.imshow(acc, extent=grid.extent, zorder=2,
                   cmap='cubehelix',
                   norm=colors.LogNorm(1, float(np.nanmax(acc))),
                   interpolation='bilinear')
    plt.colorbar(im, ax=ax, label='Upstream Cells')
    plt.title('Flow Accumulation', size=14)
    plt.xlabel('Easting')
    plt.ylabel('Northing')
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir,
                             f'flow_accumulation_{routing}.{fig_format}'),
                dpi=dpi, format=fig_format)
    plt.close()

    # - Extract river network
    if routing == 'd8':
        branches \
            = grid.extract_river_network(fdir, acc > rn_thresh,
                                         dirmap=dirmap, nodata_in=np.nan)
        sns.set_palette('husl')
        fig, ax = plt.subplots(figsize=(8.5, 6.5))
        plt.xlim(grid.bbox[0], grid.bbox[2])
        plt.ylim(grid.bbox[1], grid.bbox[3])
        ax.set_aspect('equal')

        for branch in branches['features']:
            line = np.asarray(branch['geometry']['coordinates'])
            plt.plot(line[:, 0], line[:, 1])

        _ = plt.title('D8 channels', size=14)
        plt.savefig(os.path.join(output_dir,
                                 f'flow_network_{routing}.{fig_format}'),
                    dpi=dpi, format=fig_format)
        plt.close()

        # - Define output index shapefile schema
        schema = {
            'geometry': 'LineString',
            'properties': [('id', 'int')]
        }
        out_f_name = os.path.join(output_dir, f'flow_network_{routing}.shp')

        with fiona.open(out_f_name, mode='w', driver='ESRI Shapefile',
                        schema=schema, crs="EPSG:3413") as poly_shp:
            for cnt, feat in enumerate(branches['features']):
                # -
                row_dict = {
                    # - Geometry [Polygon]
                    'geometry': {'type': 'LineString',
                                 'coordinates': feat["geometry"]["coordinates"]},
                    # - Properties [based on the schema defined above]
                    'properties': {'id': cnt},
                }
                poly_shp.write(row_dict)


if __name__ == '__main__':
    main()
