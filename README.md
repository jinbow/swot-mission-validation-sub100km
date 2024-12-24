# swot-mission-validation-sub100km
This repository contains the code and data necessary to reproduce the analyses and figures presented in Wang et al. (2025).

In Wang et al. (2025), we conducted a validation analysis of sea surface height measurements obtained from the Ka-Band Interferometer (KaRIn) radar aboard the Surface Water and Ocean Topography (SWOT) mission satellite. SWOT is a pioneering mission that utilizes KaRIn for the first time to map the elevation of water surfaces, including both inland and oceanic waters. The mission's primary strength lies in its high precision and its wide-swath coverage of 120 km, with a 20 km gap at nadir. The mission's objectives for oceanography include measuring sea surface height for wavelengths under approximately 100 km in two dimensions. This documentation describes the validation process for the mission.

For validation, we used ground truth data derived from conventional in-situ mooring platforms and confirmed that KaRIn exceeded the mission's science requirements by at least a factor of four.

This code is designed to reproduce all the analyses and figures presented in the paper using published campaign data hosted by PO.DAAC. The analysis is divided into five main steps carried out by five main scripts:

- [0.plot_figure1.py](#0plotfigure1py)
- [1.0.density_all_moorings_gliders_level-2.py](#10densityallmooringsgliderslevel-2py)
- [2.0.calculate_steric_height.py](#20calculatestericheightpy)
- [3.0.colocate.steric.karin.py](#30colocaterickarinpy)
- [4.0.process.colocated.data.py](#40processcolocateddatapy)
- [5.0.wavenumber_spectrum.py](#50wavenumberspectrumpy)

## 0.plot_figure1.py
Create Figure 1. Campaign and mooring locations on a background of SSHA from NeurOST and SWOT KaRIn. 

### Input Data Files


Place the following data files in the `data/` directory:

- **NeurOST SSHA Data**: `NeurOST_SSH-SST_20230403_20240722.nc`
- **SWOT SSHA Data**:
  - `SWOT_L2_LR_SSH_2.0_combined_calval_orbit_pass_013_sub_lat-30-40.nc`
  - `SWOT_L2_LR_SSH_2.0_combined_calval_orbit_pass_026_sub_lat-30-40.nc`
- **Mooring Positions**: `mooring_positions.csv`

These data can be downloaded from the paper page on AGU website. 

### Directory Structure

```markdown
   .
   ├── data/
   │   ├── NeurOST_SSH-SST_20230403_20240722.nc
   │   ├── SWOT_L2_LR_SSH_2.0_combined_calval_orbit_pass_013_sub_lat-30-40.nc
   │   ├── SWOT_L2_LR_SSH_2.0_combined_calval_orbit_pass_026_sub_lat-30-40.nc
   │   └── mooring_positions.csv
   ├── figures/
   │   ├── figure1_locations.png
   │   └── figure1_locations_zoom.png
   ├── utils.py
   └── main_script.py
```

### Usage


1. Ensure the input data files are in the `data/` directory.
2. Run the script to process the data and generate visualizations:

```bash
   python 0.plot_figure1.py
```

3. The output figures are saved in the `figures/` directory:
   - `figure1_locations.png`: Regional map with SSHA data and mooring positions.
   - `figure1_locations_zoom.png`: Zoomed-in map for details of the mooring positions.

4. Use photo editor such as Gimp to create the Figure 1 in Wang et al. (2025)

## 1.0.density_all_moorings_gliders_level-2.py
This script calculates the density from moorings and gliders, combines them into one file, removes outliers, and saves the resulting data.

### Input

* `path_input`: a directory path containing input files (e.g., netCDF files). The relative path should be `../data/`.
* `fn`: the output file name for the combined dataset
* Optional parameters:
	+ `do_moorings`: boolean flag to indicate whether to process mooring data (default: True)
	+ `do_remove_outliers`: boolean flag to indicate whether to remove outliers from the data (default: False)

### Process

1. **Moorings**
	* Open input files and extract relevant variables (temperature, salinity, pressure, density)
	* Calculate seawater density using the `c_rho` function
	* Create a new dataset with the calculated densities
2. **Gliders**
	* Open input files and extract relevant variables (pressure, temperature, salinity, density)
	* Apply outlier removal to each profiler on each mooring (using the `remove_outliers` function)
3. **Removing Outliers**
	* Apply a 5-sigma rule to detect outliers in the density data at each depth with 5-m bin size. 
4. **Combining Data**
	* Combine the processed mooring and glider data into a single dataset using the `xr.Dataset` class
5. **Saving Output**
	* Save the combined dataset to disk as a new netCDF file (`fn_out`)

### Output

* A single netCDF file containing all processed data, including:
	+ Moorings
	+ Gliders (with outlier removal)
* Optional: the original mooring and glider datasets are saved separately.

## 2.0.calculate_steric_height.py
This script calculates the steric height from the density anomaly profiles from moorings and gliders. It uses various functions to process the data, remove spikes, and calculate the steric height.

### Modules Used

* `xarray` for handling netCDF files
* `numpy` for numerical computations
* `pandas` for data manipulation and analysis
* `matplotlib.pyplot` for plotting (not used in this script)
* `gsw` for calculating density from pressure
* `tqdm` for displaying progress bars

### Functions Used

* `despike_profile(x, y, sigma=5)`: removes spikes from a profile
* `remove_spikes(df, sigma=5, varname='steric')`: removes spikes from a time series in a pandas DataFrame
* `calculate_steric_height_from_profile(depth, density, time, integration_depth=[-500, 0], ...)`:
	+ calculates the steric height from a profile
	+ uses the `despike_profile` function to remove spikes
	+ uses the `remove_spikes` function to remove spikes from the data
* `process_glider(glider_name, ...)`: processes glider data and returns a pandas DataFrame with the steric height
* `calculate_mooring_steric(time, data, mooring_name, ...)`: calculates the steric height from a mooring profile
* `remove_spikes(df)` (called by `calculate_mooring_steric`): removes spikes from a time series in a pandas DataFrame

### Main Script

The main script sets up various parameters for calculating the steric height, processes gliders and moorings using these parameters, and saves the results to CSV files.

### Parameters Used

* `integration_depth=[-500, 0]`: depth range for steric height calculation
* `max_spacing=10`: minimum gap threshold for identifying a good profile
* `time_span_minimum=10 * 60`: minimum time span for steric height calculation
* `time_span_maximum=120 * 60`: maximum time span for steric height calculation
* `rho0=1027.5`: density at sea level
* `surface_time_depth=-100`: depth range between surface_time_depth and surface used to estimate mean time
* `bottom_boundary_tolerance=10`: maximum depth gap on the edge of profiles at the bottom
* `top_boundary_tolerance=10`: maximum depth gap on the edge of profiles at the top
* `surface_depth=-8`: the depth of the surface layer, in which density will be replaced by the median of the layer


The paper used the following parameters:

```markdown
integration_depth=[-500, 0]
max_spacing=10
bottom_boundary_tolerance=80
top_boundary_tolerance=10
surface_depth=-8
```

### Output

The script saves the results to CSV files following the format:

* `../data/rutgers/ru32_ru38_steric_heights_depth_{-integration_depth[0]:3d}.csv`
* `../data/mooring.data/all_mooring_steric_heights_depth_{-integration_depth[0]:3d}.csv`

These files contain the steric height and other information for each glider and mooring. The following shows a snippet of the csv data file. 

**Pandas DataFrame**
=====================

| **Column Name** | **Data Type** | **Description** |
| --- | --- | --- |
| `lat` | float64 | Latitude values |
| `lon` | float64 | Longitude values |
| `time_min` | int64 | Minimum time values, seconds since 2023-01-01 |
| `time_max` | int64 | Maximum time values, seconds since 2023-01-01 |
| `surface_time` | int64 | Surface time values, seconds since 2023-01-01 |
| `depth_min` | float64 | Minimum depth values |
| `depth_max` | float64 | Maximum depth values |
| `num_points` | int64 | Number of points in each profile |
| `steric` | float64 | Steric height values (in meters) |
| `Mooring_ID` | object | Mooring ID values |

**Sample Data**
---------------

| lat  | lon   | time_min | time_max | surface_time | depth_min | depth_max | num_points | steric | Mooring_ID |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 36.1798 | -125.1268 | 5354727.0000 | 5355523.5000 | 5355434.0000 | -494.8176 | -4.2868 | 1593.0000 | 52.9755 | S1 |
| 36.1813 | -125.1264 | 5356514.5000 | 5357263.5000 | 5357178.0000 | -496.2219 | -4.7206 | 1494.0000 | 52.1012 | S1 |
| 36.1825 | -125.1261 | 5358204.5000 | 5358941.0000 | 5358857.0000 | -495.5951 | -4.6571 | 1472.0000 | 52.1746 | S1 |
| ...   | ...   | ...       | ...      | ...        | ...     | ...    | ...     | ...    | S1 |
| 36.1841 | -125.1265 | 5363242.0000 | 5363984.0000 | 5363897.0000 | -496.2865 | -4.4438 | 1486.0000 | 50.5210 | S1 |


## 3.0.colocate.steric.karin.py
Co-locate KaRIn and steric height in spatial and temporal domains

## 4.0.process.colocated.data.py
Process colocated data to valid pairs and snapshots

## 5.0.wavenumber_spectrum.py
Wavenumber spectrum analysis 
