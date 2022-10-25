# Compute Subglacial Runoff 

[![Language][]][1]

### Calculate subglacial discharge for the selected ice covered basin.

### Theory:
Evaluate discharge flow direction based on the Hydraulic Potential gradient evaluated at the interface between ice and bedrock. The definition of Hydraulic Potential and its calculation is presented in:
[Cuffey & Peterson, 2006 - "The Physics of Glaciers" - Chapter #6][].

Hydraulic Potential: phi_x = Pw + rho_water * gravity * Z   
Where:  
- Pw is Water Pressure calculated as: Pw = rho_ice * gravity * [S - Z]
  - units: [kg / m3] * [m / s2] * [m] = [kg / m / s2] 
- Z is Ice Bed Elevation
- S is Ice Surface Elevation

In this case, we used ice thickness and bedrock elevation estimates from [BedMachine version 5][].

Once the Hydraulic Potential has been computed, the accumulated flow is calculated by employing the python package [RichDEM][] Python API.

The List of flow routing algorithms available with RichDEM is available on the project webpage.

**NOTE**: by default, the D-Infinite (Dinf) routing algorithm is used to evaluate the discharge flow direction presented in:
        [Taroboton at al. 1997][].


Compute the runoff generated within each domain pixel that will be routed according the obtained flow routing map. The components of the discharge generated at each location are:

- Runoff from ice melt:
  - Estimates using output from [RACMO 2.3p2][];
  - All the runoff generated is available to generate sub-glacial discharge;

  - Geothermal Heat Melt Production calculated as:
    - *melt = G x area / rho_i / L / 2*
    - Where:
      - *G* - Geothermal heat flux cont. [J/s/m2];
      - *area* - pixel area [m2];
      - *rho_i* - density of ice Ice Density [kg / m3];
      - *L* - latent heat fusion of ice [J/kg]

  - Basal Friction Melt Production calculates as:
    - *melt = tau * Vice x area / rho_i / L / 2*
    - Where:
      - *tau* - basal shear stress equal to 100,000 kPa (1 bar);
      - *Vice* - Ice velocity expressed in meters per second;
      - Annual velocity maps from the NASA MeaSUREs project are used to estimate this component. Note that this means that we are assuming ice velocity at the interface between ice and bedrock equal to the ice surface velocity.

**Note**: In both cases, the factor 2 is included to the formula based on the assumption that only 50% of heat produced by these two processes is available for ice melt production.

### Installation:

1. Setup minimal **conda** installation using [Miniconda][]

2. Create Python Virtual Environment

    > -   Creating an environment with commands ([Link][]);
    > -   Creating an environment from an environment.yml file
    >     ([Link][2])  -> **Recommended**;


### PYTHON DEPENDENCIES:
- [RichDEM: High-Performance Terrain Analysis][]
- [rasterio: access to geospatial raster data][]
- [fiona: Fiona is GDAL’s neat and nimble vector API for Python programmers.][]
- [numpy: The fundamental package for scientific computing with Python.][]
- [xarray: Labelled multi-dimensional arrays in Python.][]
- [matplotlib: Library for creating static, animated, and interactive visualizations in Python.][]
- [tqdm: A Fast, Extensible Progress Bar for Python and CLI.][]
- [necdft4: Provides an object-oriented python interface to the netCDF version 4 library.][]

[Language]: https://img.shields.io/badge/python%20-3.7%2B-brightgreen
[1]: ..%20image::%20https://www.python.org/
[Miniconda]: https://docs.conda.io/en/latest/miniconda.html
[Link]: https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html#creating-an-environment-with-commands
[2]: https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html#creating-an-environment-from-an-environment-yml-file

[RichDEM:  High-Performance Terrain Analysis]:https://richdem.readthedocs.io/en/latest/python_api.html
[xarray: Labelled multi-dimensional arrays in Python.]:https://docs.xarray.dev
[rasterio: access to geospatial raster data]:https://rasterio.readthedocs.io/en/latest/
[matplotlib: Library for creating static, animated, and interactive visualizations in Python.]:https://matplotlib.org
[tqdm: A Fast, Extensible Progress Bar for Python and CLI.]: https://github.com/tqdm/tqdm
[necdft4: Provides an object-oriented python interface to the netCDF version 4 library.]:https://pypi.org/project/netCDF4/
[fiona: Fiona is GDAL’s neat and nimble vector API for Python programmers.] :https://fiona.readthedocs.io/en/latest/
[numpy: The fundamental package for scientific computing with Python.]:https://numpy.org

[Cuffey & Peterson, 2006 - "The Physics of Glaciers" - Chapter #6]:https://www.elsevier.com/books/the-physics-of-glaciers/cuffey/978-0-12-369461-4
[BedMachine version 5]:https://nsidc.org/data/idbmg4/versions/5
[Taroboton at al. 1997]:https://agupubs.onlinelibrary.wiley.com/doi/abs/10.1029/96WR03137
[RACMO 2.3p2]:https://www.science.org/doi/10.1126/sciadv.aaw0123