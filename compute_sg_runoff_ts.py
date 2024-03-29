#!/usr/bin/env python
"""
Written by Enrico Ciraci' - (09/2022)

CALCULATE SUB-GLACIAL DISCHARGE TIME SERIES FOR THE SELECTED ICE-COVERED BASIN

Monthly total discharge maps are obtained by following the steps reported below:

1.  Evaluate discharge flow direction based on the Hydraulic Potential gradient
    evaluated at the interface between ice and bedrock.
    The definition of Hydraulic Potential and its calculation is presented in:

    Cuffey & Peterson, 2006 - "The Physics of Glaciers" - Chapter #6.
    https://www.elsevier.com/books/the-physics-of-glaciers/
                        cuffey/978-0-12-369461-4

    Hydraulic Potential:
        phi_x = Pw + rho_water * gravity * Z
            - Where:
            - Pw is Water Pressure calculated as:
            -     Pw = rho_ice * gravity * [S - Z]
            -          units: [kg / m3] * [m / s2] * [m] = [kg / m / s2]
            -
            - Z is Ice Bed Elevation
                - S is Ice Surface Elevation

    In this case, we used ice thickness and bedrock elevation estimates
    from BedMachine version 5. Bedmachine outputs are available at this
    link: https://nsidc.org/data/idbmg4/versions/5

2.  Once the Hydraulic Potential has been computed, the accumulated flow is
    calculated by employing the python package RichDEM. Find more details at
    this link: https://richdem.readthedocs.io

    The List of flow routing algorithms available with RichDEM is available
    at this link:  https://richdem.readthedocs.io/en/latest/python_api.html

    NOTE: by default, the D-Infinite (Dinf) routing algorithm is used to
        evaluate the discharge flow direction presented in:
        Taroboton at al. 1997
        https://agupubs.onlinelibrary.wiley.com/doi/abs/10.1029/96WR03137


3.  Compute the runoff generated within each domain pixel that will be routed
    according the obtained flow routing map. The components of the discharge
    generated at each location are:

        - Runoff from ice melt:
            . Estimates using output from RACMO 2.3p2;
            . All the runoff generated is available to generate
              sub-glacial discharge;

        - Geothermal Heat Melt Production calculated as:
            . melt = G x area / rho_i / L / 2.
            . Where:
                + G - Geothermal heat flux cont. [J/s/m2];
                + area - pixel area [m2];
                + rho_i - density of ice Ice Density [kg / m3];
                + L - latent heat fusion of ice [J/kg]

        - Basal Friction Melt Production calculates as:
            .melt = tau * Vice x area / rho_i / L / 2.
            . Where:
                + tau - basal shear stress equal to 100,000 kPa (1 bar);
                + Vice - Ice velocity expressed in meters per second;
                  -> Annual velocity maps from the NASA MeaSUREs project
                     are used to estimate this component. Note that this means
                     that we are assuming ice velocity at the interface between
                     ice and bedrock equal to the ice surface velocity.

        Note: In both cases, the factor 2 is included to the formula based on
            the assumption that only 50% of heat produced by these two processes
            is available for ice melt production.

usage: compute_sg_runoff_ts.py [-h] [--directory DIRECTORY] [--routing ROUTING]
        [--domain DOMAIN] [--f_year F_YEAR] [--l_year L_YEAR]

Compute sub-glacial discharge monthly time series over the selected domain.

optional arguments:
  -h, --help            show this help message and exit
  --directory DIRECTORY, -S DIRECTORY
                        Data Directory.
  --routing ROUTING, -R ROUTING
                        Flow routing method.
  --domain DOMAIN, -D DOMAIN
                        Region of Interest [Integration Domain].
  --f_year F_YEAR, -F F_YEAR
                        First Year of the selected period [def. 1958].
  --l_year L_YEAR, -L L_YEAR
                        Last Year of the selected period [def. 2021].
  --save_monthly, -M    Save estimated sub-glacial discharge monthly maps.
  --netcdf, -N          Save estimated sub-glacial discharge monthly maps
                        in NetCDF format.

To Run this script digit:

$ python compute_sg_runoff_ts.py ...


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
geopandas: open source project to make working with geospatial data in python.
    https://geopandas.org
netcdf4: netcdf4-python is a Python interface to the netCDF C library.
    https://unidata.github.io/netcdf4-python/
xarray: N-D labeled arrays and datasets in Python.
    https://docs.xarray.dev/en/stable/
"""
# - Python Dependencies
from __future__ import print_function
import os
import argparse
import numpy as np
import matplotlib.pyplot as plt
import geopandas as gpd
import netCDF4 as nC4
import richdem as rd
import xarray as xr
import datetime
from pyproj import CRS, Transformer
# - Library Dependencies
from utils.make_dir import make_dir
from utils.utility_functions_rasterio import \
    load_raster, save_raster, clip_raster, sample_in_memory_dataset
from utils.mpl_utils import add_colorbar
from utils.load_velocity_map import load_velocity_map_nearest
from utils.xyscale_north import xyscale_north
# -
plt.rc('font', family='monospace')
plt.rc('font', weight='bold')
plt.style.use('seaborn-v0_8-deep')


def main() -> None:
    # - Read the system arguments listed after the program
    parser = argparse.ArgumentParser(
        description="""Compute sub-glacial discharge monthly time series 
        over the selected domain."""
    )
    # - Optional Arguments
    project_dir = os.environ['PYTHONDATA']
    parser.add_argument('--directory', '-S', type=str,
                        default=project_dir, help='Data Directory.')

    parser.add_argument('--routing', '-R', type=str, default='Dinf',
                        help='Flow routing method.')

    parser.add_argument('--domain', '-D', type=str,
                        default='Petermann_Drainage_Basin_EPSG3413',
                        help='Region of Interest [Integration Domain].')

    parser.add_argument('--f_year', '-F', type=int, default=1958,
                        help='First Year of the selected period [def. 1958].')

    parser.add_argument('--l_year', '-L', type=int, default=2021,
                        help='Last Year of the selected period [def. 2021].')

    parser.add_argument('--save_monthly', '-M', action='store_true',
                        help='Save estimated sub-glacial discharge '
                             'monthly maps.')

    parser.add_argument('--netcdf', '-N', action='store_true',
                        help='Save estimated sub-glacial discharge monthly maps'
                             ' in NetCDF format.')

    # - Processing Parameters
    args = parser.parse_args()

    # - Conversion Parameters
    rho_water = 1000        # - Water Density [kg / m3]
    rho_ice = 917           # - Ice Density [kg / m3]
    gravity = 9.81          # - Gravity Acceleration [m / s2]
    # - Other Parameters
    # - Geothermal Heat Melt Production - melt = G x area / rho_i / L / 2.
    geothermal_const = 51e-3    # - geothermal heat flux cont. [J/s/m**2]
    latent_heat_ice = 334e3     # - latent heat fusion of ice [J/kg]

    # - Basal Friction Melt Production
    # - melt = 1e5 x v (in meters per second) x area / rho_i / L / 2.
    tau = 1e5
    n_sec_year = 365 * 24 * 60 * 60
    n_sec_month = 30 * 24 * 60 * 60

    # - Pixel Area [m2]
    pixel_area = 150 * 150

    # - Figure Parameters
    fig_format = 'jpeg'
    dpi = 300

    # - Project Data Directory
    project_dir = args.directory
    bdmch_dir = os.path.join(project_dir, 'BedMachine', 'output.dir')
    output_dir = make_dir(bdmch_dir, 'subglacial_runoff_petermann_ts')
    output_dir = make_dir(output_dir, args.domain)
    output_dir = make_dir(output_dir, args.routing)

    # - Domain Mask
    clip_mask = os.path.join(project_dir, 'GIS_Data', 'Petermann_Domain',
                             f'{args.domain}.shp')

    # - Shapefile containing the points used to sample the total
    # - discharge at the grounding line.
    try:
        sample_pts_path\
            = os.path.join(project_dir, 'GIS_Data',
                           'Petermann_features_extraction',
                           f'subglacial_runoff_sample_{args.routing}.shp')
        sample_pts_in = gpd.read_file(sample_pts_path)
    except FileNotFoundError:
        sample_pts_path\
            = os.path.join(project_dir, 'GIS_Data',
                           'Petermann_features_extraction',
                           f'subglacial_runoff_sample.shp')
        sample_pts_in = gpd.read_file(sample_pts_path)

    sample_ps_iter = []
    for pt in sample_pts_in.geometry:
        sample_ps_iter.append((pt.xy[0][0], pt.xy[1][0]))

    # - Absolute Path to BedMachine Data
    # - NOTE - these estimates must be computed before running the script.
    # - Bedrock
    bed_path \
        = os.path.join(bdmch_dir, 'BedMachineGreenland-version_05_bed',
                       f'{args.domain}_bed_EPSG-3413_res150_bilinear.tiff')
    # - Ice Elevation
    elev_path \
        = os.path.join(bdmch_dir, 'BedMachineGreenland-version_05_surface',
                       f'{args.domain}_surface_EPSG-3413_res150_bilinear.tiff')

    # - Ice Thickness
    thick_path \
        = os.path.join(bdmch_dir, 'BedMachineGreenland-version_05_thickness',
                       f'{args.domain}_thickness_EPSG-3413_res150'
                       f'_bilinear.tiff')

    # - Compute Hydraulic Potential
    print('# - Compute Hydraulic Potential.')
    bed_input = load_raster(bed_path)
    elev_input = load_raster(elev_path)
    thick_input = load_raster(thick_path)
    # - compute ice mask from ice thickness
    ice_mask = np.zeros(np.shape(thick_input['data']))
    ice_mask[thick_input['data'] > 0] = 1

    # For more Info about this derivation, see:
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
    crs_epsg = crs.to_epsg()
    res = elev_input['res']

    # - Hydraulic Potential
    hydro_pot = ((rho_ice * gravity * surf_elev)
                 + ((rho_water - rho_ice) * gravity * bed_elev))

    # - Consider only ice-covered regions
    hydro_pot[ice_mask == 0] = np.nan
    # - Hydraulic Potential in GeoTiff format
    out_path = os.path.join(output_dir, f'{args.domain}_hydro_pot_res150.tiff')
    save_raster(hydro_pot, res, x_vect.copy(), y_vect.copy(), out_path, crs)

    # - Clip Hydraulic Potential over the actual drainage basin
    out_path_clip = out_path.replace('_hydro_pot',
                                     '_hydro_pot_clipped')
    clip_raster(out_path, clip_mask, out_path_clip, nodata=np.nan)

    # - Load Clipped Mask
    hydro_pot_in = load_raster(out_path_clip)
    hydro_pot_clp = hydro_pot_in['data']
    x_vect = hydro_pot_in['x_coords']
    y_vect = hydro_pot_in['y_coords']
    extent = elev_input['extent']
    hydro_pot = hydro_pot_clp
    # - Need to use flipud here because when loaded using RichDEM
    # - y-axis direction of is inverted.
    ice_mask_clp = np.where(np.isnan(np.flipud(hydro_pot_clp)))

    print('# - Compute Accumulated Sub-Glacial Flow using RichDEM.')
    dem = rd.LoadGDAL(out_path_clip, no_data=-9999)
    dem[ice_mask_clp] = np.nan

    # - Plot Clipped Hydraulic Potential
    # - Show Hydraulic Potential
    fig = plt.figure(figsize=(6, 9), dpi=dpi)
    fig.patch.set_alpha(0)
    ax = fig.add_subplot(1, 1, 1)
    im = ax.imshow(dem, extent=extent, cmap='plasma', zorder=1,
                   interpolation='bilinear')
    ax.grid(color='m', linestyle='dotted', alpha=0.3)
    ax.set_title(f'Hydraulic Potential Map', size=14)
    ax.set_xlabel('Easting')
    ax.set_ylabel('Northing')

    cb_1 = add_colorbar(fig, ax, im)
    cb_1.set_label(label=r'[$kg/m/s^2$]', weight='bold')
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir,
                             f'hydraulic_potential_map_'
                             f'{args.routing}.{fig_format}'),
                dpi=dpi, format=fig_format)
    plt.close()

    # - Initialize Raster Domain Corners
    x_min = None
    x_max = None
    y_min = None
    y_max = None
    ice_mask_crp = None
    x_vect_crp = None
    y_vect_crp = None
    # - Initialize RichDEM Object
    dem_c = None
    # - Initialize Melt Component
    out_dir_mnth = None

    # - Create Time Series vector and time axis
    time_ax = []
    gl_discharge = []

    # - Crate directory that will contain discharge monthly estimates
    if args.save_monthly:
        out_dir_mnth = make_dir(output_dir, 'monthly_est')

    for year in range(args.f_year, args.l_year+1):

        # - Load Ice Velocity Maps for the considered year
        if year > 2013:
            # - Load Interpolate velocity Map
            # - If selected, apply smoothing filter to
            # - the interpolated maps.
            v_map = load_velocity_map_nearest(year, 1, 1,
                                              args.directory,
                                              domain=args.domain,
                                              verbose=True)
        else:
            # - Velocity maps for years before 2014 are characterized
            # - by high noise level and discontinuities. Do not use
            # - these maps. Use velocity from 2014 based on the
            # - assumption that ice velocity does not change
            # - significantly over time.
            v_map = load_velocity_map_nearest(2014, 6, 1,
                                              args.directory,
                                              domain=args.domain,
                                              verbose=True)
        # - Velocity x and y components
        velocity_x = v_map['vx_out']
        velocity_y = v_map['vy_out']
        velocity_x_vect = v_map['x']
        velocity_y_vect = v_map['y']

        # - Compute velocity magnitude
        velocity_mag = np.sqrt((velocity_x ** 2) + (velocity_y ** 2))

        for month in range(1, 13):
            print(f'# - {year} - {month}')
            # - Load Runoff estimates from RACMO 2.3p2
            racmo_path \
                = os.path.join(project_dir, 'SMB', 'RACMO2.3p2', 'output.dir',
                               f'{args.domain}', 'runoff',
                               f'{args.domain}_runoff_EPSG-{crs_epsg}_res150',
                               f'{year}')
            racmo_f_name \
                = [x for x in os.listdir(racmo_path)
                   if x.startswith(f'runoff.{year}-{month:02}')
                   and x.endswith('150m_nearest.tiff')][0]
            runoff_input = load_raster(os.path.join(racmo_path, racmo_f_name))
            runoff = runoff_input['data']
            runoff_x_vect = runoff_input['x_coords']
            runoff_y_vect = runoff_input['y_coords']

            if (year == year) & (month == 1):
                # - Crop the two datasets over the overlapping domain
                x_min = np.max([runoff_x_vect[0], x_vect[0],
                                velocity_x_vect[0]])
                x_max = np.min([runoff_x_vect[-1], x_vect[-1],
                                velocity_x_vect[-1]])
                y_min = np.max([runoff_y_vect[0], y_vect[0],
                                velocity_y_vect[0]])
                y_max = np.min([runoff_y_vect[-1], y_vect[-1],
                                velocity_y_vect[-1]])
                # - Hydraulic Potential
                ind_hp_x = np.where((x_vect >= x_min) & (x_vect <= x_max))
                ind_hp_y = np.where((y_vect >= y_min) & (y_vect <= y_max))
                ind_hp_xx, ind_hp_yy = np.meshgrid(ind_hp_x, ind_hp_y)
                hydro_pot = hydro_pot[ind_hp_yy, ind_hp_xx]
                x_vect_crp = x_vect[ind_hp_x]
                y_vect_crp = y_vect[ind_hp_y]

                # - Save the Obtained Rasters
                # - Hydraulic potential
                out_hp_crp \
                    = os.path.join(output_dir,
                                   f'{args.domain}_hydro_pot'
                                   f'_EPSG-{crs_epsg}_res150_crop.tiff')
                # - NOTE: Hydraulic Potential Map is saved here with the y-axis
                # -       orientation inverted. This operation is necessary to
                # -       correctly compute discharge using RichDEM Python API.
                # -       It is not clear to me why this step is necessary.
                # -       This raster is overwritten with a correctly oriented
                # -       one at the end if the script, once the computation of
                # -       discharge has been completed.
                save_raster(np.flipud(hydro_pot), res,
                            x_vect_crp.copy(), y_vect_crp.copy(),
                            out_hp_crp, crs)
                ice_mask_crp = np.where(np.isnan(hydro_pot))

                # -------------------------------------------
                # - Initialize RichDEM Object
                # -------------------------------------------
                # - Calculate flow accumulation with weights
                # -------------------------------------------
                dem_c = rd.LoadGDAL(out_hp_crp, no_data=-9999)
                # - Fill depressions with epsilon gradient to ensure drainage
                rd.FillDepressions(dem_c, epsilon=True, in_place=True)

                if month == 1:
                    # - Ice Velocity Magnitude
                    ind_v_x = np.where(
                        (velocity_x_vect >= x_min) & (velocity_x_vect <= x_max))
                    ind_v_y = np.where(
                        (velocity_y_vect >= y_min) & (velocity_y_vect <= y_max))
                    ind_v_xx, ind_v_yy = np.meshgrid(ind_v_x, ind_v_y)

                    velocity_mag = velocity_mag[ind_v_yy, ind_v_xx]

            # - Calculates the scaling factor for a polar stereographic
            # - projection to correct area calculations.

            # - compute raster coordinates mesh-grids
            m_xx, m_yy = np.meshgrid(x_vect_crp.copy(), y_vect_crp.copy())
            crs_4326 = CRS.from_epsg(4326)  # - input projection
            crs_3413 = CRS.from_epsg(3413)  # - input projection
            transformer = Transformer.from_crs(crs_3413, crs_4326)
            # - transform raster coordinates to lat-lon
            m_lat, m_lon = transformer.transform(m_yy, m_xx)
            # - Pixel Area scaling factor
            px_area_sc = xyscale_north(m_lat * np.pi / 180.0)

            # - Compute Melt Production components
            # - Geothermal Heat Melt Production -
            # -   > melt = G x area / rho_i / L / 2.
            geothermal_melt \
                = (((geothermal_const * pixel_area) / rho_ice)
                   / latent_heat_ice) / 2.
            geothermal_melt \
                = np.ones(np.shape(velocity_mag)) * geothermal_melt

            # - Basal Friction Melt Production
            friction_melt \
                = tau * ((((velocity_mag / n_sec_year) * pixel_area)
                          / rho_ice) / latent_heat_ice) / 2.

            # - Runoff
            ind_rf_x = np.where(
                (runoff_x_vect >= x_min) & (runoff_x_vect <= x_max))
            ind_rf_y = np.where(
                (runoff_y_vect >= y_min) & (runoff_y_vect <= y_max))
            ind_rf_xx, ind_rf_yy = np.meshgrid(ind_rf_x, ind_rf_y)
            runoff = runoff[ind_rf_yy, ind_rf_xx]

            # - Runoff - convert from mm.WE/month to m.W.E./sec
            runoff_melt = ((runoff / 1e3) / n_sec_month) * pixel_area

            # - Total discharge produced within each pixel - [m3/sec]
            weights \
                = (geothermal_melt + friction_melt + runoff_melt) * px_area_sc

            # - Get flow accumulation using runoff generated per pixel
            # - as weight.
            accum_dw \
                = rd.FlowAccumulation(dem_c, method=args.routing,
                                      weights=weights)

            accum_dw[ice_mask_crp] = np.nan
            if args.save_monthly:
                if args.netcdf:
                    # - Save  monthly discharge Raster in NetCDF4 format
                    # - Calculate Latitude/Longitude coordinates for
                    # - each grid point centroid.
                    # - More Info here:
                    # -     https://daac.ornl.gov/submit/netcdfrequirements/
                    x_coords_c = x_vect_crp.copy() + (res[0] / 2.)
                    y_coords_c = y_vect_crp.copy() + (res[1] / 2.)

                    ds_flow = xr.Dataset(data_vars=dict(
                        disch=(["y", "x"], accum_dw)),
                        coords=dict(x=(["x"], x_coords_c),
                                    y=(["y"], y_coords_c))
                    )
                    # - Dataset Attributes
                    ds_flow.attrs['long_name'] = 'Sub-glacial Discharge'
                    ds_flow.attrs['standard_name'] = 'discharge'
                    ds_flow.attrs['unit'] = 'm3/sec'
                    ds_flow.attrs['routing_algorithm'] = args.routing
                    # - Add Coordinates Reference System Information
                    ds_flow.attrs['ellipsoid'] = 'WGS84'
                    ds_flow.attrs['false_easting'] = 0.0
                    ds_flow.attrs['false_northing'] = 0.0
                    ds_flow.attrs['grid_mapping_name'] = "polar_stereographic"
                    ds_flow.attrs['longitude_of_projection_origin'] = 45.0
                    ds_flow.attrs['latitude_of_projection_origin'] = 90.0
                    ds_flow.attrs['standard_parallel'] = 70.0
                    ds_flow.attrs['straight_vertical_longitude_from_pole'] = 0.0
                    ds_flow.attrs['EPSG'] = crs_epsg
                    ds_flow.attrs['spatial_resolution'] = '150m'

                    # - save the cropped velocity field
                    out_path \
                        = os.path.join(out_dir_mnth,
                                       f'sub_glacial_discharge_map_{year}'
                                       f'-{month:02}_{args.routing}.nc4')
                    ds_flow.to_netcdf(out_path, format='NETCDF4',
                                      engine='netcdf4',
                                      encoding={
                                          'disch': dict(zlib=True, complevel=9,
                                                        chunksizes=(100, 100),
                                                        dtype='float64',
                                                        )})
                else:
                    # - Save  monthly discharge Raster in GeoTiff format
                    out_path \
                        = os.path.join(out_dir_mnth,
                                       f'sub_glacial_discharge_map_{year}'
                                       f'-{month:02}_{args.routing}.tiff')
                    save_raster(accum_dw, res,
                                x_vect_crp.copy(), y_vect_crp.copy(),
                                out_path, crs, nodata=np.nan)

            # - Sample the computed Accumulated floe
            # - Sample the computed Accumulated flow at the points of interest
            sampled_pts = sample_in_memory_dataset(accum_dw, res,
                                                   x_vect_crp.copy(),
                                                   y_vect_crp.copy(),
                                                   crs, sample_ps_iter,
                                                   nodata=np.nan)
            total_discharge = np.sum(sampled_pts)
            print(f'# - Total Discharge at GL [m3/sec]: {total_discharge:.2f}')

            time_ax.append(datetime.datetime(year=year, month=month, day=1))
            gl_discharge.append(total_discharge)

    # - save the obtained discharge time series as a netcdf archive
    file_to_save \
        = os.path.join(output_dir, f'{args.domain}_sub_glacial_discharge_'
                                   f'[m3*sec-1]_routing-{args.routing}'
                                   f'_{args.f_year}'
                                   f'_{args.l_year}.nc')
    time_ax_out = []
    for tm in range(len(time_ax)):
        time_ax_out.append(
            (datetime.date(int(time_ax[tm].year), int(time_ax[tm].month),
                           int(time_ax[tm].day)) - datetime.date(2000, 1,
                                                                 1)).days)

    # - Writing the Output file
    rootgrp = nC4.Dataset(file_to_save, mode='w', format='NETCDF4')
    time_unit = 'days since 2000-01-01 00:00'
    calendar = 'standard'
    # - Create Variable dimensions
    rootgrp.createDimension('time', len(time_ax_out))
    # - Create output variables:
    var_time = rootgrp.createVariable('time', 'f4', 'time')
    var_time.units = time_unit
    var_time.calendar = calendar
    var_disch_ts = rootgrp.createVariable('Discharge', 'f4', 'time')
    var_disch_ts.units = 'm3/sec'
    # -
    var_time[:] = time_ax_out
    var_disch_ts[:] = np.array(gl_discharge)
    var_disch_ts.actual_range = [np.min(gl_discharge), np.max(gl_discharge)]
    var_disch_ts.standard_name = 'Total Sub-glacial Discharge at ' \
                                 'the Grounding Line'
    # -
    rootgrp.close()

    # - Overwrite previously saved Hydraulic potential Map.
    out_hp_crp \
        = os.path.join(output_dir, f'{args.domain}_hydro_pot'
                                   f'_EPSG-{crs_epsg}_res150_crop.tiff')
    save_raster(hydro_pot, res, x_vect_crp.copy(), y_vect_crp.copy(),
                out_hp_crp, crs)


if __name__ == '__main__':
    main()
