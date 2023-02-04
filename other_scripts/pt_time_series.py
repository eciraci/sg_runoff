#!/usr/bin/env python
"""
Written by Enrico Ciraci' - (02/2023)

Compute sub-glacial discharge time series at the selected location and compare
with the total discharge time series measured at the grounding line of the
glacier estimated running the script compute_sg_runoff_monthly.py
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

    parser.add_argument('--pt_sample', '-P', type=str,
                        default='subglacial_runoff_sample_main_channel',
                        help='Shapefile containing the location of the point '
                             'of interest.')

    # - Processing Parameters
    args = parser.parse_args()

    # - Path to directory containing monthly sub-glacial discharge maps
    project_dir = args.directory
    bdmch_dir = os.path.join(project_dir, 'BedMachine', 'output.dir')
    monthly_data_dir \
        = os.path.join(bdmch_dir, 'subglacial_runoff_petermann_ts',
                       args.domain, args.routing, 'monthly_est')

    # - Shapefile containing the points used to sample the total
    # - discharge at the grounding line.
    sample_gl_path \
        = os.path.join(project_dir, 'GIS_Data',
                       'Petermann_features_extraction',
                       f'subglacial_runoff_sample_{args.routing}.shp')
    sample_gl_in = gpd.read_file(sample_gl_path)
    sample_gl_iter = []
    for pt in sample_gl_in.geometry:
        sample_gl_iter.append((pt.xy[0][0], pt.xy[1][0]))

    # - Shapefile containing the points used to sample the total
    # - discharge at the selected location
    sample_pts_path \
        = os.path.join(project_dir, 'GIS_Data',
                       'Petermann_features_extraction',
                       f'{args.pt_sample}_{args.routing}.shp')
    sample_pts_in = gpd.read_file(sample_pts_path)

    # - Output Directory
    output_dir \
        = make_dir(os.path.join(bdmch_dir, 'subglacial_runoff_petermann_ts',
                                args.domain, args.routing), 'pt_sample')

    sample_ps_iter = []
    for pt in sample_pts_in.geometry:
        sample_ps_iter.append((pt.xy[0][0], pt.xy[1][0]))

    # - time series vector and time axis
    gl_discharge = []      # - Total discharge at the grounding line
    pt_discharge = []       # - Total discharge at the selected location
    time_ax = []

    # - Compute Monthly Time Series
    for year in range(args.f_year, args.l_year + 1):
        for month in range(1, 13):
            # - Loop over the months
            m_file \
                = os.path.join(monthly_data_dir,
                               f'sub_glacial_discharge_map_{year}-{month:02d}'
                               f'_{args.routing}.tiff')
            # - Load the monthly sub-glacial discharge map
            sdisch_in = load_raster(m_file)
            sdisch_map = sdisch_in['data']
            sdisch_x_vect = sdisch_in['x_coords']
            sdisch_y_vect = sdisch_in['y_coords']
            crs = sdisch_in['crs']
            res = sdisch_in['res']

            # - Sample the computed Accumulated flow
            sampled_pts = sample_in_memory_dataset(sdisch_map, res,
                                                   sdisch_x_vect.copy(),
                                                   sdisch_y_vect.copy(),
                                                   crs, sample_gl_iter,
                                                   nodata=np.nan)
            total_discharge_gl = np.sum(sampled_pts)

            # - Sample the computed Accumulated flow at the points of interest
            sampled_pts = sample_in_memory_dataset(sdisch_map, res,
                                                   sdisch_x_vect.copy(),
                                                   sdisch_y_vect.copy(),
                                                   crs, sample_ps_iter,
                                                   nodata=np.nan)
            total_discharge_pt = np.sum(sampled_pts)

            # - Print the total discharge
            print(f'# - {year} - {month}:')
            print(f'# - Total Discharge at GL [m3/sec]: '
                  f'{total_discharge_gl:.2f}')
            print(f'# - Total Discharge at Pts [m3/sec]: '
                  f'{total_discharge_pt:.2f}')
            print(f'# - Total Discharge Ratio: '
                  f'{total_discharge_pt/total_discharge_gl:.4f}')

            time_ax.append(datetime.datetime(year=year, month=month, day=1))
            gl_discharge.append(total_discharge_gl)
            pt_discharge.append(total_discharge_pt)

    # - save the obtained discharge time series as a netcdf archive
    file_to_save \
        = os.path.join(output_dir, f'{args.pt_sample}_'
                                   f'[m3*sec-1]_routing-{args.routing}'
                                   f'_{args.f_year}'
                                   f'_{args.l_year}.nc')
    time_ax_out = []
    for tm in range(len(time_ax)):
        time_ax_out.append(
            (datetime.date(int(time_ax[tm].year), int(time_ax[tm].month),
                           int(time_ax[tm].day)) - datetime.date(2000, 1,
                                                                 1)).days)

    # - Writing the Output files
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
    var_disch_ts[:] = np.array(pt_discharge)
    var_disch_ts.actual_range = [np.min(pt_discharge), np.max(pt_discharge)]
    var_disch_ts.standard_name = 'Total Sub-glacial Discharge at ' \
                                 'the Selected Location'
    # -
    rootgrp.close()

    # - Save total discharge at the GL time series
    file_to_save \
        = os.path.join(output_dir, f'{args.domain}_sub_glacial_discharge_'
                                   f'[m3*sec-1]_routing-{args.routing}'
                                   f'_{args.f_year}'
                                   f'_{args.l_year}.nc')

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


if __name__ == '__main__':
    main()
