from __future__ import print_function
import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors
import richdem as rd
import rasterio
import rasterio.mask
from rasterio.transform import Affine
from typing import Union, Any
# - Library Dependencies
from utils.make_dir import make_dir
from utils.utility_functions_rasterio import load_raster, clip_raster


def save_raster(raster: np.ndarray, res: Any, x: np.ndarray,
                y: np.ndarray, out_path: str, crs: int,
                nbands: int = 1, nodata: int = -9999) -> None:
    """
    Save the Provided Raster in GeoTiff format
    :param raster: input raster - np.ndarray
    :param res: raster resolution - integer
    :param x: x-axis - np.ndarray
    :param y: y-axis in a figure coordinate system - np.ndarray
    :param crs: - coordinates reference system
    :param out_path: absolute path to output file
    :param nbands: number of raster bands
    :param nodata: output raster no data value
    :return: None
    """
    # - Calculate Affine Transformation of the output raster
    if y[1] > y[0]:
        y_vect = np.flipud(y)
        raster = np.flipud(raster)
        y_vect += res[1]
    print(y_vect[0])
    transform = (Affine.translation(x[0], y_vect[0])
                 * Affine.scale(res[0], -res[1]))

    out_meta = {'driver': 'GTiff',
                'height': raster.shape[0],
                'width': raster.shape[1],
                'nodata': nodata,
                'dtype': str(raster.dtype),
                'compress': 'lzw',
                'count': nbands,
                'crs': crs,
                'transform': transform}

    with rasterio.open(out_path, 'w', **out_meta) as dst:
        dst.write(raster, 1)


def main() -> None:
    # - Conversion Parameters
    rho_water = 1000            # - Water Density [kg / m3]
    rho_ice = 917               # - Ice Density [kg / m3]
    gravity = 9.81              # - Gravity Acceleration [m / s2]

    project_dir = os.path.join('/', 'Volumes', 'Extreme Pro')
    # - Domain Mask
    clip_mask = os.path.join(project_dir, 'GIS_Data', 'Petermann_Domain',
                             f'Petermann_Drainage_Basin_EPSG3413.shp')

    bdmch_dir = os.path.join(project_dir, 'BedMachine', 'output.dir')
    domain = 'Petermann_Drainage_Basin_EPSG3413'
    # - Absolute Path to BedMachine Data
    # - NOTE - these estimates must be computed before running the script.
    # - Bedrock
    bed_path \
        = os.path.join(bdmch_dir, 'BedMachineGreenland-version_05_bed',
                       f'{domain}_bed_EPSG-3413_res150_bilinear.tiff')
    # - Ice Elevation
    elev_path \
        = os.path.join(bdmch_dir, 'BedMachineGreenland-version_05_surface',
                       f'{domain}_surface_EPSG-3413_res150_bilinear.tiff')

    # - Ice Thickness
    thick_path \
        = os.path.join(bdmch_dir, 'BedMachineGreenland-version_05_thickness',
                       f'{domain}_thickness_EPSG-3413_res150'
                       f'_bilinear.tiff')

    output_dir = os.path.join(os.path.expanduser("~"), 'Desktop')

    # - Compute Hydraulic Potential
    print('# - Compute Hydraulic Potential.')
    bed_input = load_raster(bed_path)
    elev_input = load_raster(elev_path)
    thick_input = load_raster(thick_path)
    # - compute ice mask from ice thickness
    thick_mask = thick_input['data']
    ice_mask = np.zeros(np.shape(thick_mask))
    ice_mask[thick_mask > 0] = 1

    surf_elev = elev_input['data']
    bed_elev = bed_input['data']
    x_vect = elev_input['x_coords']
    y_vect = elev_input['y_coords']
    crs = elev_input['crs']
    crs_epsg = crs.to_epsg()
    res = elev_input['res']
    extent = elev_input['extent']
    out_path = os.path.join(output_dir, f'surface_elev_res150.tiff')
    save_raster(surf_elev, res, x_vect.copy(), y_vect.copy(), out_path, crs)

    print(y_vect[0], y_vect[-1])
    # - Hydraulic Potential
    # hydro_pot = ((rho_ice * gravity * surf_elev)
    #              + ((rho_water - rho_ice) * gravity * bed_elev))
    hydro_pot = surf_elev

    # - Consider only ice-covered regions
    # hydro_pot[ice_mask == 0] = np.nan
    # -
    out_path = os.path.join(output_dir, f'hydro_pot_res150.tiff')
    print(y_vect[0], y_vect[-1])
    save_raster(hydro_pot, res, x_vect.copy(), y_vect.copy(), out_path, crs)
    out_path = os.path.join(output_dir, f'hydro_pot_res150_2.tiff')
    print(y_vect[0], y_vect[-1])
    save_raster(hydro_pot, res, x_vect.copy(), y_vect.copy(), out_path, crs)
    plt.imshow(hydro_pot)
    plt.show()
    #
    # # - Load Runoff estimates from RACMO 2.3p2
    # year = 2017
    # month = 2
    # racmo_path \
    #     = os.path.join(project_dir, 'SMB', 'RACMO2.3p2', 'output.dir',
    #                    f'{domain}', 'runoff',
    #                    f'{domain}_runoff_EPSG-{crs_epsg}_res150',
    #                    f'{year}',
    #                    f'runoff.{year}-{month:02}.BN_RACMO2.3p2'
    #                    f'_ERA5_3h_FGRN055.150m_nearest.tiff')
    # runoff_input = load_raster(racmo_path)
    # runoff = runoff_input['data']
    # runoff_x_vect = runoff_input['x_coords']
    # runoff_y_vect = runoff_input['y_coords']
    # # plt.imshow(runoff)
    # # plt.show()
    # # print(runoff_input['src_transform'])
    #
    # # - Ice Velocity Maps
    # velocity_path \
    #     = os.path.join(project_dir, 'Greenland_Ice_Velocity_MEaSUREs',
    #                    f'{domain}', 'interp_vmaps_res150',
    #                    'vel_2017-07-01_2018-06-31',
    #                    'vel_2017-07-01_2018-06-31-rio_EPSG-3413_'
    #                    'res-150_average.tiff')
    # velocity_input = load_raster(velocity_path, nbands=2)
    # velocity = velocity_input['data']
    # velocity_x_vect = velocity_input['x_coords']
    # velocity_y_vect = velocity_input['y_coords']
    # # - Compute velocity magnitude
    # velocity_mag = np.sqrt((velocity[0, :, :] ** 2) + (velocity[1, :, :] ** 2))
    #
    # out_path = os.path.join(output_dir, f'velocity_mag_res150.tiff')
    # save_raster(velocity_mag, res, velocity_x_vect, velocity_y_vect,
    #             out_path, crs)
    # # plt.imshow(velocity_mag)
    # # plt.show()
    #
    # print([runoff_y_vect[0], y_vect[0], velocity_y_vect[0]])
    # print([runoff_y_vect[-1], y_vect[-1], velocity_y_vect[-1]])
    # print(' ')
    # print([runoff_x_vect[0], x_vect[0], velocity_x_vect[0]])
    # print([runoff_x_vect[-1], x_vect[-1], velocity_x_vect[-1]])
    #
    # # - Crop the two datasets over the overlapping domain
    # x_min = np.max([runoff_x_vect[0], x_vect[0], velocity_x_vect[0]])
    # x_max = np.min([runoff_x_vect[-1], x_vect[-1], velocity_x_vect[-1]])
    # y_min = np.max([runoff_y_vect[0], y_vect[0], velocity_y_vect[0]])
    # y_max = np.min([runoff_y_vect[-1], y_vect[-1], velocity_y_vect[-1]])
    # # - Hydraulic Potential
    # ind_hp_x = np.where((x_vect >= x_min) & (x_vect <= x_max))
    # ind_hp_y = np.where((y_vect >= y_min) & (y_vect <= y_max))
    # ind_hp_xx, ind_hp_yy = np.meshgrid(ind_hp_x, ind_hp_y)
    # hydro_pot = hydro_pot[ind_hp_yy, ind_hp_xx]
    # x_vect_crp = x_vect[ind_hp_x]
    # y_vect_crp = y_vect[ind_hp_y]
    #
    # print(y_vect_crp[0], y_vect_crp[-1])
    # print(y_vect[0], y_vect_crp[-1])
    #
    # # - Runoff
    # ind_rf_x = np.where((runoff_x_vect >= x_min) & (runoff_x_vect <= x_max))
    # ind_rf_y = np.where((runoff_y_vect >= y_min) & (runoff_y_vect <= y_max))
    # ind_rf_xx, ind_rf_yy = np.meshgrid(ind_rf_x, ind_rf_y)
    # runoff =  runoff[ind_rf_yy, ind_rf_xx]
    # #
    # # - Ice Velocity Magnitude
    # ind_v_x = np.where((velocity_x_vect >= x_min) & (velocity_x_vect <= x_max))
    # ind_v_y = np.where((velocity_y_vect >= y_min) & (velocity_y_vect <= y_max))
    # ind_v_xx, ind_v_yy = np.meshgrid(ind_v_x, ind_v_y)
    # velocity_mag = velocity_mag[ind_v_yy, ind_v_xx]
    #
    # # - Save the Obtained Rasters
    # # - Hydraulic potential
    # out_hp_crp \
    #     = os.path.join(output_dir, f'hydro_pot_res150_crop.tiff')
    # save_raster(hydro_pot, res, x_vect_crp, y_vect_crp,
    #             out_hp_crp, crs)
    # ice_mask_crp = np.where(np.isnan(np.flipud(hydro_pot)))
    #
    # # - runoff
    # out_rf_crp \
    #     = os.path.join(output_dir, f'runoff_res150_crop.tiff')
    # save_raster(runoff, res, x_vect_crp, y_vect_crp, out_rf_crp, crs)
    #
    # # - velocity magnitude
    # out_vel_crp \
    #     = os.path.join(output_dir, f'velocity_mag_res150_crop.tiff')
    # save_raster(velocity_mag, res, x_vect_crp, y_vect_crp,
    #             out_vel_crp, crs)


if __name__ == '__main__':
    main()



