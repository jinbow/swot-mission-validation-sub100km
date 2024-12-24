# swot-mission-validation-sub100km
This repository contains the code and data necessary to reproduce the analyses and figures presented in Wang et al. (2025).

In Wang et al. (2025), we conducted a validation analysis of sea surface height measurements obtained from the Ka-Band Interferometer (KaRIn) radar aboard the Surface Water and Ocean Topography (SWOT) mission satellite. SWOT is a pioneering mission that utilizes KaRIn for the first time to map the elevation of water surfaces, including both inland and oceanic waters. The mission's primary strength lies in its high precision and its wide-swath coverage of 120 km, with a 20 km gap at nadir. The mission's objectives for oceanography include measuring sea surface height for wavelengths under approximately 100 km in two dimensions. This documentation describes the validation process for the mission.

For validation, we used ground truth data derived from conventional in-situ mooring platforms and confirmed that KaRIn exceeded the mission's science requirements by at least a factor of four.

This code is designed to reproduce all the analyses and figures presented in the paper using published campaign data hosted by PO.DAAC. The analysis is divided into five main steps as outlined below.

## Create Figure 1. Campaign and mooring locations on a SWOT background

## Restructure all mooring profiles and calculate density anomaly

## Calculate steric height from mooring profiles

## Co-locate KaRIn and steric height in spatial and temporal domains

## Process colocated data to valid pairs and snapshots

## Wavenumber spectrum analysis 
