# swot-mission-validation-sub100km
This repository contains the code and data necessary to reproduce the analyses and figures presented in Wang et al. (2025).

In Wang et al. (2025), we conducted a validation analysis of sea surface height measurements obtained from the Ka-Band Interferometer (KaRIn) radar aboard the Surface Water and Ocean Topography (SWOT) mission satellite. SWOT is a pioneering mission that utilizes KaRIn for the first time to map the elevation of water surfaces, including both inland and oceanic waters. The mission's primary strength lies in its high precision and its wide-swath coverage of 120 km, with a 20 km gap at nadir. The mission's objectives for oceanography include measuring sea surface height for wavelengths under approximately 100 km in two dimensions. This documentation describes the validation process for the mission.

For validation, we used ground truth data derived from conventional in-situ mooring platforms and confirmed that KaRIn exceeded the mission's science requirements by at least a factor of four.

This code is designed to reproduce all the analyses and figures presented in the paper using published campaign data hosted by PO.DAAC. The analysis is divided into five main steps as outlined below.

## 0.plot_figure1.py
Create Figure 1. Campaign and mooring locations on a background of SSHA from NeurOST and SWOT KaRIn. 

Input Data Files
----------------

Place the following data files in the `data/` directory:

- **NeurOST SSHA Data**: `NeurOST_SSH-SST_20230403_20240722.nc`
- **SWOT SSHA Data**:
  - `SWOT_L2_LR_SSH_2.0_combined_calval_orbit_pass_013_sub_lat-30-40.nc`
  - `SWOT_L2_LR_SSH_2.0_combined_calval_orbit_pass_026_sub_lat-30-40.nc`
- **Mooring Positions**: `mooring_positions.csv`

These data can be downloaded from the paper page on AGU website. 

Usage
-----

1. Ensure the input data files are in the `data/` directory.
2. Run the script to process the data and generate visualizations:

.. code-block:: bash

   python 0.plot_figure1.py

3. The output figures are saved in the `figures/` directory:
   - `figure1_locations.png`: Regional map with SSHA data and mooring positions.
   - `figure1_locations_zoom.png`: Zoomed-in map for detailed analysis.

4. Use photo editor such as Gimp to create the Figure 1 in Wang et al. (2025)

Directory Structure
-------------------

.. code-block:: text

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

## 1.0.density_all_moorings_gliders_level-2.py
Restructure all mooring profiles and calculate density anomaly

## 2.0.calculate_steric_height.py
Calculate steric height from mooring profiles

## 3.0.colocate.steric.karin.py
Co-locate KaRIn and steric height in spatial and temporal domains

## 4.0.process.colocated.data.py
Process colocated data to valid pairs and snapshots

## 5.0.wavenumber_spectrum.py
Wavenumber spectrum analysis 
