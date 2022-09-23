#!/usr/bin/env python
"""
Written by Enrico Ciraci' - (09/2022)

Test Scripts:
Compute sub-glacial discharge flow by employing ice thickness and bed topography
from BedMachine v4.0. Discharge flow is determined by the spatial gradient of
the Hydraulic Potential computed at the interface between ice and bedrock.

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
      RichDEM project is used -> https://richdem.readthedocs.io

To Run this script digit:

$ python compute_flowacc_richdem.py


PYTHON DEPENDENCIES:
RichDEM:High-Performance Terrain Analysis
    https://richdem.readthedocs.io
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
# - Python Dependencies
from __future__ import print_function
import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors
import richdem as rd
# - Library Dependencies
from utils.make_dir import make_dir
from utils.utility_functions_rasterio import \
    load_raster, save_raster, clip_raster


def main() -> None:
    # - Processing Parameters
    routing = 'Dinf'     # - Selected Routing Algorithm
    # - For other accumulation routing schemes:
    # - https://richdem.readthedocs.io/en/latest/python_api.html

    # - Figure Parameters
    fig_format = 'jpeg'
    dpi = 200

    # - Project Data Directory
    project_dir = os.path.join('/', 'Volumes', 'Extreme Pro', 'BedMachine',
                               'output.dir')
    output_dir = make_dir(project_dir, 'subglacial_flow_rich')
    output_dir = make_dir(output_dir, routing)

    # - Absolute Path to BedMachine Data
    # - Bedrock
    bed_path \
        = os.path.join('/', 'Volumes', 'Extreme Pro', 'BedMachine',
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
    extent = elev_input['extent']

    # - Hydraulic Potential
    hydro_pot = ((rho_ice * gravity * surf_elev)
                 + ((rho_water - rho_ice) * gravity * bed_elev))

    # - Consider only ice-covered regions
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

    # - Load Clipped Mask
    hydro_pot_in = load_raster(out_path_clip)
    hydro_pot_clp = np.flipud(hydro_pot_in['data'])
    ice_mask_clp = np.where(np.isnan(hydro_pot_clp))

    print('# - Compute Accumulated Sub-Glacial Flow using RichDEM.')
    dem = rd.LoadGDAL(out_path_clip, no_data=-9999)
    dem[ice_mask_clp] = np.nan
    # - Show Hydraulic Potential
    fig, ax = plt.subplots(figsize=(8, 6))
    fig.patch.set_alpha(0)
    plt.imshow(dem, extent=extent, cmap='plasma', zorder=1,
               interpolation='bilinear')
    plt.colorbar(label='[kg / m / s2]')
    plt.grid(zorder=0)
    plt.title('Hydraulic Potential Map', size=14)
    plt.xlabel('Longitude')
    plt.ylabel('Latitude')
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir,
                             f'hydraulic_potential.{fig_format}'),
                dpi=dpi, format=fig_format)
    plt.show()
    plt.close()

    # - Fill depressions with epsilon gradient to ensure drainage
    rd.FillDepressions(dem, epsilon=True, in_place=True)
    # - Get flow accumulation with no explicit weighting. The default will be 1.
    accum_d = rd.FlowAccumulation(dem, method=routing)
    d_fig = rd.rdShow(accum_d, zxmin=450, zxmax=550, zymin=550, zymax=450,
                      figsize=(8, 5.5), axes=False, cmap='jet')
    plt.close()

    # - Compute terrain attributed
    slope = rd.TerrainAttribute(dem, attrib='aspect')
    rd.rdShow(slope, axes=False, cmap='jet', figsize=(8, 5.5))

    # ---------------------------
    # Calculate flow accumulation
    # ---------------------------
    accum_d[ice_mask_clp] = np.nan
    # - Save Accumulated Flow Map
    out_path = os.path.join(output_dir, f'flow_accumulation_{routing}.tiff')
    save_raster(np.flipud(accum_d), res, x_vect, y_vect, out_path, crs,
                nodata=np.nan)
    # - Save the Logarithm of Accumulated Flow Map
    out_path = os.path.join(output_dir, f'flow_accumulation_{routing}_log.tiff')
    save_raster(np.flipud(np.log(accum_d)), res, x_vect, y_vect, out_path, crs)

    # - Plot Accumulated Flow
    fig, ax = plt.subplots(figsize=(8, 6))
    fig.patch.set_alpha(0)
    plt.grid('on', zorder=0)
    im = ax.imshow(accum_d, extent=extent, zorder=2,
                   cmap='cubehelix',
                   norm=colors.LogNorm(1, float(np.nanmax(accum_d))),
                   interpolation='bilinear')
    plt.colorbar(im, ax=ax, label='# Upstream Cells')
    plt.title(f'Log(Flow Accumulation) - {routing}', size=14)
    plt.xlabel('Easting')
    plt.ylabel('Northing')
    plt.savefig(os.path.join(output_dir,
                             f'flow_accumulation_{routing}_log.{fig_format}'),
                dpi=dpi, format=fig_format)
    plt.close()


if __name__ == '__main__':
    main()
