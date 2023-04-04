# MR_B0B1

## Summary
This module performs analysis of B0/B1 measurements of the large body phantom. Implementation according to department legacy scripts

## Status
Initial version released (20220630). Some cleanup is still needed.

## Dependencies
pip:
- numpy
- pydicom
- matplotlib
- skimage
- cv2 (opencv-contrib-python)

## Acquisition protocol
Use a scan that produces a B0 directly on the scanner. And additionally, the B1 acquisition is 2 different scans with different flip angles (60, 120).

## Selector
A selector for this module should run on dcm_study level as it takes data from the 5 different scans ((gantry) B0 map (with and without shim), B1 60 and B1 120). 

### Input data
dcm_study

### Config
- mr_philips_b0b1.json
Contains filters to select proper acquisitions for each test, based on ImageType and SeriesDescription

### Meta
Limits & analysis are based on https://doi.org/10.3389/fphy.2022.891600

### Rules
StationName, PatientName

## Analysis
- B0 (with and without shim)
	- Select B0 maps
	- Create cilinder (35 cm diameter)
	- Calculate statistics in cilinder
- B1
	- Select scans (flip angle 60 and 120)
	- Calculate B1 map
	- Calculate statistics in 35 cm circle (central slice)

## Results
- B0 figures
	- B0 maps from central slice
- B0 statistics
	- Average, standard deviation, p99-p1 range
- B1 figures
	- Flipangle 60 and 120, B1 map, B1map with mask, Profiles through B1 map
- B1 statistics 
	-  Average, standard deviation, p99-p1 range
