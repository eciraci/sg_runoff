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


To Run this script digit:

$ python compute_sg_runoff.py ...


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
import argparse
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors
import geopandas as gpd
import richdem as rd
import rasterio
# - Library Dependencies
from utils.make_dir import make_dir
from utils.utility_functions_rasterio import \
    load_raster, save_raster, clip_raster
from utils.mpl_utils import add_colorbar
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

    # - Processing Parameters
    domain_name = 'Petermann_Drainage_Basin'       # - integration domain
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
    n_sec_year = 365*24*60*68
    n_sec_month = 30*24*60*68

    # - Pixel Area [m2]
    pixel_area = 150*150

    # - Figure Parameters
    fig_format = 'jpeg'
    dpi = 300

    # - Project Data Directory
    project_dir = args.directory
    bdmch_dir = os.path.join(project_dir, 'BedMachine', 'output.dir')
    output_dir = make_dir(bdmch_dir, f'subglacial_runoff_petermann')
    output_dir = make_dir(output_dir, args.domain)
    output_dir = make_dir(output_dir, f'{args.month:02}-{args.year}')
    output_dir = make_dir(output_dir, args.routing)

    # - Domain Mask
    clip_mask = os.path.join(project_dir, 'GIS_Data', 'Petermann_Domain',
                             f'Petermann_Drainage_Basin_EPSG3413.shp')

    # - Shapefile containing the points used to sample the total
    # - discharge at the grounding line.
    try:
        sample_pts_path\
            = os.path.join(project_dir, 'GIS_Data', 'Petermann_features_extraction',
                           f'subglacial_runoff_sample_{args.routing}.shp')
        sample_pts_in = gpd.read_file(sample_pts_path)
    except:
        sample_pts_path\
            = os.path.join(project_dir, 'GIS_Data', 'Petermann_features_extraction',
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
    extent = elev_input['extent']

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

    # - Load Runoff estimates from RACMO 2.3p2
    racmo_path \
        = os.path.join(project_dir, 'SMB', 'RACMO2.3p2', 'output.dir',
                       f'{args.domain}', 'runoff',
                       f'{args.domain}_runoff_EPSG-{crs_epsg}_res150',
                       f'{args.year}',
                       f'runoff.{args.year}-{args.month:02}.BN_RACMO2.3p2'
                       f'_ERA5_3h_FGRN055.150m_nearest.tiff')
    runoff_input = load_raster(racmo_path)
    runoff = runoff_input['data']
    runoff_x_vect = runoff_input['x_coords']
    runoff_y_vect = runoff_input['y_coords']

    # - Ice Velocity Maps
    velocity_path \
        = os.path.join(project_dir, 'Greenland_Ice_Velocity_MEaSUREs',
                       f'{args.domain}', 'interp_vmaps_res150',
                       'vel_2017-07-01_2018-06-31',
                       'vel_2017-07-01_2018-06-31-rio_EPSG-3413_'
                       'res-150_average.tiff')
    velocity_input = load_raster(velocity_path, nbands=2)
    velocity = velocity_input['data']
    velocity_x_vect = velocity_input['x_coords']
    velocity_y_vect = velocity_input['y_coords']

    # - Compute velocity magnitude
    velocity_mag = np.sqrt((velocity[0, :, :]**2) + (velocity[1, :, :]**2))

    # - Crop the two datasets over the overlapping domain
    x_min = np.max([runoff_x_vect[0], x_vect[0], velocity_x_vect[0]])
    x_max = np.min([runoff_x_vect[-1], x_vect[-1], velocity_x_vect[-1]])
    y_min = np.max([runoff_y_vect[0], y_vect[0], velocity_y_vect[0]])
    y_max = np.min([runoff_y_vect[-1], y_vect[-1], velocity_y_vect[-1]])

    # - Hydraulic Potential
    ind_hp_x = np.where((x_vect >= x_min) & (x_vect <= x_max))
    ind_hp_y = np.where((y_vect >= y_min) & (y_vect <= y_max))
    ind_hp_xx, ind_hp_yy = np.meshgrid(ind_hp_x, ind_hp_y)
    hydro_pot = np.flipud(hydro_pot[ind_hp_yy, ind_hp_xx])
    x_vect_crp = x_vect[ind_hp_x]
    y_vect_crp = y_vect[ind_hp_y]

    # - Runoff
    ind_rf_x = np.where((runoff_x_vect >= x_min) & (runoff_x_vect <= x_max))
    ind_rf_y = np.where((runoff_y_vect >= y_min) & (runoff_y_vect <= y_max))
    ind_rf_xx, ind_rf_yy = np.meshgrid(ind_rf_x, ind_rf_y)
    runoff = runoff[ind_rf_yy, ind_rf_xx]

    # - Ice Velocity Magnitude
    ind_v_x = np.where((velocity_x_vect >= x_min) & (velocity_x_vect <= x_max))
    ind_v_y = np.where((velocity_y_vect >= y_min) & (velocity_y_vect <= y_max))
    ind_v_xx, ind_v_yy = np.meshgrid(ind_v_x, ind_v_y)
    velocity_mag = velocity_mag[ind_v_yy, ind_v_xx]

    # - Save the Obtained Rasters
    # - Hydraulic potential
    out_hp_crp \
        = os.path.join(output_dir, f'{args.domain}_hydro_pot'
                                   f'_EPSG-{crs_epsg}_res150_crop.tiff')
    save_raster(hydro_pot, res, x_vect_crp.copy(), y_vect_crp.copy(),
                out_hp_crp, crs)
    ice_mask_crp = np.where(np.isnan(np.flipud(hydro_pot)))

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
    save_raster(weights, res, x_vect_crp.copy(), y_vect_crp.copy(),
                out_weights, crs)
    out_weights_clp \
        = os.path.join(output_dir, f'{args.domain}_weights'
                                   f'_EPSG-{crs_epsg}_res150_clip.tiff')
    clip_raster(out_weights, clip_mask, out_weights_clp, nodata=np.nan)

    # -------------------------------------------
    # - Calculate flow accumulation with weights
    # -------------------------------------------
    dem_c = rd.LoadGDAL(out_hp_crp, no_data=-9999)
    # - Fill depressions with epsilon gradient to ensure drainage
    rd.FillDepressions(dem_c, epsilon=True, in_place=True)
    # - Get flow accumulation with no explicit weighting. The default will be 1.
    accum_dw = rd.FlowAccumulation(dem_c, method=args.routing, weights=weights)
    rd.rdShow(np.flipud(accum_dw), zxmin=450, zxmax=550, zymin=550, zymax=450,
              figsize=(6, 9), axes=False, cmap='jet')
    plt.close()

    # - Compute terrain attributed
    slope = rd.TerrainAttribute(dem, attrib='aspect')
    rd.rdShow(slope, axes=False, cmap='jet', figsize=(6, 9))
    plt.close()

    accum_dw[ice_mask_crp] = np.nan
    # - Save Sub-Glacier Discharges Map
    out_path \
        = os.path.join(output_dir, f'sub_glacial_discharge_map_'
                                   f'{args.routing}.tiff')
    save_raster(accum_dw, res, x_vect.copy(), y_vect.copy(), out_path, crs,
                nodata=np.nan)

    # - Sample the computed Accumulated floe
    src = rasterio.open(out_path)
    sample_pts_in['value'] = [x for x in src.sample(sample_ps_iter)]
    total_discharge = sample_pts_in['value'].sum()[0]
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


if __name__ == '__main__':
    main()
