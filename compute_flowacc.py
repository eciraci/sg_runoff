import os
import pyproj
import numpy as np
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors
import seaborn as sns
from pysheds.grid import Grid


def main() -> None:
    # - path to sample dataset
    data_path = os.path.join('/', 'Volumes', 'Extreme Pro', 'BedMachine',
                             'output.dir', 'BedMachineGreenland_bed',
                             'Petermann_Domain_Velocity_Stereo_bed_EPSG'
                             '-3413_res150_average.tiff')

    # - Read Datasets using pyshed.grid
    grid = Grid.from_raster(data_path, data_name='dem')
    dem = grid.read_raster(data_path)

    # - Show Data
    fig, ax = plt.subplots(figsize=(8, 6))
    fig.patch.set_alpha(0)

    plt.imshow(dem, extent=grid.extent, cmap='terrain',
               zorder=1, vmin=-500, vmax=500)
    plt.colorbar(label='Elevation (m)')
    plt.grid(zorder=0)
    plt.title('Digital elevation map', size=14)
    plt.xlabel('Longitude')
    plt.ylabel('Latitude')
    plt.tight_layout()
    plt.show()

    # Condition DEM
    # ----------------------
    # Fill pits in DEM
    pit_filled_dem = grid.fill_pits(dem)

    # Fill depressions in DEM
    flooded_dem = grid.fill_depressions(pit_filled_dem)

    # Resolve flats in DEM
    inflated_dem = grid.resolve_flats(flooded_dem)


    # - Compute Elevation Flow Direction
    # Determine D8 flow directions from DEM
    # ----------------------
    # Specify directional mapping
    dirmap = (64, 128, 1, 2, 4, 8, 16, 32)

    # Compute flow directions
    # -------------------------------------
    fdir = grid.flowdir(inflated_dem, dirmap=dirmap)


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
    plt.show()

    # Calculate flow accumulation
    # --------------------------
    acc = grid.accumulation(fdir, dirmap=dirmap)

    fig, ax = plt.subplots(figsize=(8, 6))
    fig.patch.set_alpha(0)
    plt.grid('on', zorder=0)
    im = ax.imshow(acc, extent=grid.extent, zorder=2,
                   cmap='cubehelix',
                   norm=colors.LogNorm(1, acc.max()),
                   interpolation='bilinear')
    plt.colorbar(im, ax=ax, label='Upstream Cells')
    plt.title('Flow Accumulation', size=14)
    plt.xlabel('Easting')
    plt.ylabel('Northing')
    plt.tight_layout()
    plt.show()

    # Extract river network
    # ---------------------
    branches = grid.extract_river_network(fdir, acc > 50, dirmap=dirmap)

    sns.set_palette('husl')
    fig, ax = plt.subplots(figsize=(8.5, 6.5))

    plt.xlim(grid.bbox[0], grid.bbox[2])
    plt.ylim(grid.bbox[1], grid.bbox[3])
    ax.set_aspect('equal')

    for branch in branches['features']:
        line = np.asarray(branch['geometry']['coordinates'])
        plt.plot(line[:, 0], line[:, 1])

    _ = plt.title('D8 channels', size=14)
    plt.show()


if __name__ == '__main__':
    main()
