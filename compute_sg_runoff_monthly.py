#!/usr/bin/env python
"""
Written by Enrico Ciraci' - (09/2022)

CALCULATE SUB-GLACIAL DISCHARGE FOR THE SELECTED ICE-COVERED BASIN

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

usage: compute_sg_runoff.py [-h] [--directory DIRECTORY] [--routing ROUTING]
            [--domain DOMAIN] [--year YEAR] [--month MONTH]

Compute sub-glacial discharge flow over the selected domain.

options:
  -h, --help            show this help message and exit
  --directory DIRECTORY, -S DIRECTORY
                        Data Directory.
  --routing ROUTING, -R ROUTING
                        Flow routing method.
  --domain DOMAIN, -D DOMAIN
                        Region of Interest [Integration Domain].
  --year YEAR, -Y YEAR  Year of interest [def. 2017].
  --month MONTH, -m MONTH
                        Month of interest [def. 7 - July].
  --netcdf, -N          Save estimated sub-glacial discharge in NetCDF format.
  --show                Show intermediate results.

To Run this script digit:

$ python compute_sg_runoff_monthly.py ...

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
scikit-image: Image processing in Python.
    https://scikit-image.org/
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
import richdem as rd
import fiona
import xarray as xr
from skimage.morphology import skeletonize
# - Library Dependencies
from utils.make_dir import make_dir
from utils.utility_functions_rasterio import \
    load_raster, save_raster, clip_raster, sample_in_memory_dataset
from utils.mpl_utils import add_colorbar
from utils.load_velocity_map import load_velocity_map_nearest
# -
plt.rc('font', family='monospace')
plt.rc('font', weight='bold')
plt.style.use('seaborn-deep')


def main() -> None:
    # - Read the system arguments listed after the program
    parser = argparse.ArgumentParser(
        description="""Compute sub-glacial discharge flow over the 
        selected domain."""
    )
    # - Optional Arguments
    project_dir = os.path.join('/', 'Volumes', 'Extreme Pro')
    parser.add_argument('--directory', '-S', type=str,
                        default=project_dir,  help='Data Directory.')

    parser.add_argument('--routing', '-R', type=str, default='Dinf',
                        help='Flow routing method.')

    parser.add_argument('--domain', '-D', type=str,
                        default='Petermann_Drainage_Basin_EPSG3413',
                        help='Region of Interest [Integration Domain].')

    parser.add_argument('--year', '-Y', type=int, default=2017,
                        help='Year of interest [def. 2017].')

    parser.add_argument('--month', '-M', type=int, default=7,
                        help='Month of interest [def. 7 - July].')

    parser.add_argument('--netcdf', '-N', action='store_true',
                        help='Save estimated sub-glacial discharge '
                             'in NetCDF format.')

    parser.add_argument('--show', action='store_true',
                        help='Show intermediate results.')

    # - Processing Parameters
    args = parser.parse_args()

    # - Conversion Parameters
    rho_water = 1000            # - Water Density [kg / m3]
    rho_ice = 917               # - Ice Density [kg / m3]
    gravity = 9.81              # - Gravity Acceleration [m / s2]
    # - Other Parameters
    # - Geothermal Heat Melt Production - melt = G x area / rho_i / L / 2.
    geothermal_const = 51e-3     # - geothermal heat flux cont. [J/s/m**2]
    latent_heat_ice = 334e3      # - latent heat fusion of ice [J/kg]

    # - Basal Friction Melt Production
    # - melt = 1e5 x v (in meters per second) x area / rho_i / L / 2.
    tau = 1e5
    n_sec_year = 365 * 24 * 60 * 68
    n_sec_month = 30 * 24 * 60 * 68

    # - Pixel Area [m2]
    pixel_area = 150 * 150

    # - Figure Parameters
    fig_format = 'jpeg'
    dpi = 300
    # - Thresholds used to:
    bin_thresh = 0.05        # - generate discharge binary mask
    sk_thresh = 1            # - generate discharge skeletonized map

    # - Project Data Directory
    project_dir = args.directory
    bdmch_dir = os.path.join(project_dir, 'BedMachine', 'output.dir')
    output_dir = make_dir(bdmch_dir, f'subglacial_runoff_petermann')
    output_dir = make_dir(output_dir, args.domain)
    output_dir = make_dir(output_dir, f'{args.month:02}-{args.year}')
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
    except fiona.errors.DriverError:
        sample_pts_path\
            = os.path.join(project_dir, 'GIS_Data',
                           'Petermann_features_extraction',
                           f'subglacial_runoff_sample.shp')
        sample_pts_in = gpd.read_file(sample_pts_path)

    sample_ps_iter = []
    for pt in sample_pts_in.geometry:
        sample_ps_iter.append((pt.xy[0][0], pt.xy[1][0]))

    # - Absolute Path to Interpolated and Cropped BedMachine Data
    # - NOTE - Need to compute the data listed below before running the script.
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
    # -
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
    # - Add Colorbar
    cb_1 = add_colorbar(fig, ax, im)
    cb_1.set_label(label=r'[$kg/m/s^2$]', weight='bold')
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir,
                             f'hydraulic_potential_map_'
                             f'{args.routing}.{fig_format}'),
                dpi=dpi, format=fig_format)
    plt.close()

    # - Load Runoff estimates from RACMO 2.3p2
    racmo_path \
        = os.path.join(project_dir, 'SMB', 'RACMO2.3p2', 'output.dir',
                       f'{args.domain}', 'runoff',
                       f'{args.domain}_runoff_EPSG-{crs_epsg}_res150',
                       f'{args.year}')
    racmo_f_name \
        = [x for x in os.listdir(racmo_path)
           if x.startswith(f'runoff.{args.year}-{args.month:02}')
           and x.endswith('150m_nearest.tiff')][0]
    racmo_path = os.path.join(racmo_path, racmo_f_name)

    runoff_input = load_raster(racmo_path)
    runoff = runoff_input['data']
    runoff_x_vect = runoff_input['x_coords']
    runoff_y_vect = runoff_input['y_coords']

    # - Ice Velocity Maps
    if args.year > 2013:
        # - Load Interpolate velocity Map
        # - If selected, apply smoothing filter to the interpolated maps.
        v_map = load_velocity_map_nearest(args.year, args.month, 1,
                                          args.directory,
                                          domain=args.domain, verbose=True)
    else:
        # - Velocity maps for years before 2014 are characterized
        # - by high noise level and discontinuities. Do not use these maps.
        # - Use velocity from 2014 based on the assumption that ice velocity
        # - does not change significantly over time.
        v_map = load_velocity_map_nearest(2014, 6, 1,
                                          args.directory,
                                          domain=args.domain, verbose=True)
    # - Velocity x and y components
    velocity_x = v_map['vx_out']
    velocity_y = v_map['vy_out']
    velocity_x_vect = v_map['x']
    velocity_y_vect = v_map['y']

    # - Compute velocity magnitude
    velocity_mag = np.sqrt((velocity_x**2) + (velocity_y**2))

    # - Crop the two datasets over the overlapping domain
    x_min = np.max([runoff_x_vect[0], x_vect[0], velocity_x_vect[0]])
    x_max = np.min([runoff_x_vect[-1], x_vect[-1], velocity_x_vect[-1]])
    y_min = np.max([runoff_y_vect[0], y_vect[0], velocity_y_vect[0]])
    y_max = np.min([runoff_y_vect[-1], y_vect[-1], velocity_y_vect[-1]])

    # - Hydraulic Potential
    ind_hp_x = np.where((x_vect >= x_min) & (x_vect <= x_max))
    ind_hp_y = np.where((y_vect >= y_min) & (y_vect <= y_max))
    ind_hp_xx, ind_hp_yy = np.meshgrid(ind_hp_x, ind_hp_y)
    hydro_pot = hydro_pot[ind_hp_yy, ind_hp_xx]
    x_vect_crp = x_vect[ind_hp_x]
    y_vect_crp = y_vect[ind_hp_y]

    # - Ice Velocity Magnitude
    ind_v_x = np.where((velocity_x_vect >= x_min) & (velocity_x_vect <= x_max))
    ind_v_y = np.where((velocity_y_vect >= y_min) & (velocity_y_vect <= y_max))
    ind_v_xx, ind_v_yy = np.meshgrid(ind_v_x, ind_v_y)
    velocity_mag = velocity_mag[ind_v_yy, ind_v_xx]

    # - Runoff
    ind_rf_x = np.where((runoff_x_vect >= x_min) & (runoff_x_vect <= x_max))
    ind_rf_y = np.where((runoff_y_vect >= y_min) & (runoff_y_vect <= y_max))
    ind_rf_xx, ind_rf_yy = np.meshgrid(ind_rf_x, ind_rf_y)
    runoff = runoff[ind_rf_yy, ind_rf_xx]

    print(f'# - Average Velocity Magnitude over '
          f'the considered area: {np.nanmean(velocity_mag ):.3f}')

    # - Save the Obtained Rasters
    # - Hydraulic potential
    out_hp_crp \
        = os.path.join(output_dir, f'{args.domain}_hydro_pot'
                                   f'_EPSG-{crs_epsg}_res150_crop.tiff')
    # - NOTE: Hydraulic Potential Map is saved here with the y-axis orientation
    # -       inverted. This operation is necessary to correctly compute
    # -       discharge using RichDEM Python API.
    # -       It is not clear to me why this step is necessary.
    # -       This raster is overwritten with a correctly oriented one at the
    # -       end if the script, once the computation of discharge has been
    # -       completed.
    save_raster(np.flipud(hydro_pot), res, x_vect_crp.copy(), y_vect_crp.copy(),
                out_hp_crp, crs)

    ice_mask_crp = np.where(np.isnan(hydro_pot))

    # - runoff
    out_rf_crp \
        = os.path.join(output_dir, f'{args.domain}_runoff'
                                   f'_EPSG-{crs_epsg}_res150_crop.tiff')
    save_raster(runoff, res, x_vect_crp.copy(), y_vect_crp.copy(),
                out_rf_crp, crs)

    # - velocity magnitude
    out_vel_crp \
        = os.path.join(output_dir, f'{args.domain}_velocity_magnitude'
                                   f'_EPSG-{crs_epsg}_res150_crop.tiff')
    save_raster(velocity_mag, res, x_vect_crp.copy(), y_vect_crp.copy(),
                out_vel_crp, crs)

    # - Compute Melt Production components
    # - Geothermal Heat Melt Production - melt = G x area / rho_i / L / 2.
    geothermal_melt \
        = (((geothermal_const * pixel_area) / rho_ice) / latent_heat_ice) / 2.
    geothermal_melt = np.ones(np.shape(runoff)) * geothermal_melt

    # - Basal Friction Melt Production
    friction_melt \
        = tau * ((((velocity_mag/n_sec_year) * pixel_area)
                  / rho_ice) / latent_heat_ice) / 2.

    # - Runoff - convert from mm.WE/month to m.W.E./sec
    runoff_melt = ((runoff/1e3)/n_sec_month) * pixel_area

    # - Total discharge produced within each pixel - [m3/sec]
    weights = geothermal_melt + friction_melt + runoff_melt

    # - Weights
    out_weights \
        = os.path.join(output_dir, f'{args.domain}_weights'
                                   f'_EPSG-{crs_epsg}_res150_crop.tiff')
    save_raster(weights.copy(), res, x_vect_crp.copy(), y_vect_crp.copy(),
                out_weights, crs)
    out_weights_clp \
        = os.path.join(output_dir, f'{args.domain}_weights'
                                   f'_EPSG-{crs_epsg}_res150_clip.tiff')
    clip_raster(out_weights, clip_mask, out_weights_clp, nodata=np.nan)

    # - Origin Flow before accumulation
    fig = plt.figure(figsize=(6, 10), dpi=dpi)
    ax = fig.add_subplot(1, 1, 1)
    im = ax.imshow(np.flipud(weights), extent=extent, zorder=2,
                   cmap='jet', interpolation='bilinear',
                   vmin=0, vmax=0.008)

    ax.grid(color='m', linestyle='dotted', alpha=0.3)
    ax.set_title(f'Sub-Glacial Discharge - {args.routing}', size=14)
    ax.set_xlabel('Easting')
    ax.set_ylabel('Northing')
    cb_1 = add_colorbar(fig, ax, im)
    cb_1.set_label(label=r'[${m^3}/sec$]', weight='bold')
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir,
                             f'ice_melt_map_'
                             f'{args.routing}.{fig_format}'),
                dpi=dpi, format=fig_format)
    plt.close()

    # -------------------------------------------
    # - Calculate flow accumulation with weights
    # -------------------------------------------
    dem_c = rd.LoadGDAL(out_hp_crp, no_data=-9999)
    weights = rd.rdarray(weights, no_data=0)
    # - Fill depressions with epsilon gradient to ensure drainage
    rd.FillDepressions(dem_c, epsilon=True, in_place=True)
    # - Get flow accumulation with no explicit weighting. The default will be 1.
    accum_dw = rd.FlowAccumulation(dem_c, method=args.routing, weights=weights)

    if args.show:
        rd.rdShow(np.flipud(accum_dw), zxmin=450, zxmax=550, zymin=550,
                  zymax=450, figsize=(6, 9), axes=False, cmap='jet')
        plt.close()

    accum_dw[ice_mask_crp] = np.nan
    # - Save Sub-Glacier Discharges Map
    out_path \
        = os.path.join(output_dir, f'sub_glacial_discharge_map_'
                                   f'{args.routing}.tiff')
    save_raster(accum_dw, res, x_vect_crp.copy(), y_vect_crp.copy(),
                out_path, crs, nodata=np.nan)

    # - Sample the computed Accumulated flow at the points of interest
    sampled_pts = sample_in_memory_dataset(accum_dw, res, x_vect_crp.copy(),
                                           y_vect_crp.copy(),
                                           crs, sample_ps_iter, nodata=np.nan)
    total_discharge = np.sum(sampled_pts)
    print(f'# - Total Discharge at GL [m3/sec]: {total_discharge:.2f}')

    # - Plot Accumulated Flow
    fig = plt.figure(figsize=(6, 10), dpi=dpi)
    ax = fig.add_subplot(1, 1, 1)
    im = ax.imshow(np.flipud(accum_dw), extent=extent, zorder=2,
                   cmap='Reds',  interpolation='bilinear',
                   vmin=0, vmax=0.008)

    ax.grid(color='m', linestyle='dotted', alpha=0.3)
    ax.set_title(f'Sub-Glacial Discharge - {args.routing}', size=14)
    ax.set_xlabel('Easting')
    ax.set_ylabel('Northing')
    ann_txt = f'Total Discharge at GL [m3/sec]: {total_discharge:.2f}'
    ax.annotate(ann_txt, xy=(0.03, 0.03), xycoords="axes fraction",
                size=12, zorder=100,
                bbox=dict(boxstyle="square", fc="w"))

    cb_1 = add_colorbar(fig, ax, im)
    cb_1.set_label(label=r'[${m^3}/sec$]', weight='bold')
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir,
                             f'sub_glacial_discharge_map_'
                             f'{args.routing}.{fig_format}'),
                dpi=dpi, format=fig_format)
    plt.close()

    # - Overwrite previously saved Hydraulic potential Map.
    out_hp_crp \
        = os.path.join(output_dir, f'{args.domain}_hydro_pot'
                                   f'_EPSG-{crs_epsg}_res150_crop.tiff')
    save_raster(hydro_pot, res, x_vect_crp.copy(), y_vect_crp.copy(),
                out_hp_crp, crs)

    # - Generate a binary map with values different from zero only
    # - for channels with discharge values above a certain threshold.
    dich_bin = np.full(np.shape(accum_dw), np.nan)
    dich_bin[accum_dw >= bin_thresh] = 1
    # -
    out_hp_crp \
        = os.path.join(output_dir, f'sub_glacial_discharge_map_binary_'
                                   f'{args.routing}.tiff')
    save_raster(dich_bin, res, x_vect_crp.copy(), y_vect_crp.copy(),
                out_hp_crp, crs)

    # - Skeletonize the binary map
    dich_bin = np.zeros(np.shape(accum_dw))
    dich_bin[accum_dw >= sk_thresh] = 1
    skeleton = skeletonize(dich_bin)
    skeleton_bin = np.full(np.shape(accum_dw), np.nan)
    skeleton_bin[np.array(skeleton) == 1] = 1

    # - Export discharge binary map
    out_hp_crp \
        = os.path.join(output_dir, f'sub_glacial_discharge_map_binary_skeleton_'
                                   f'{args.routing}.tiff')
    save_raster(skeleton_bin, res, x_vect_crp.copy(), y_vect_crp.copy(),
                out_hp_crp, crs)

    if args.netcdf:
        # - Save  monthly discharge Raster in NetCDF4 format
        # - Calculate Latitude/Longitude coordinates for
        # - each grid point centroid.
        # - More Info here: https://daac.ornl.gov/submit/netcdfrequirements/
        x_coords_c = x_vect_crp.copy() + (res[0]/2.)
        y_coords_c = y_vect_crp.copy() + (res[1]/2.)

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
            = os.path.join(output_dir,
                           f'sub_glacial_discharge_map_{args.routing}.nc4')
        ds_flow.to_netcdf(out_path, format='NETCDF4', engine='netcdf4',
                          encoding={'disch': dict(zlib=True, complevel=9,
                                                  chunksizes=(100, 100),
                                                  dtype='float64',
                                                  )})


if __name__ == '__main__':
    main()
