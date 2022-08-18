#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 27 10:26:28 2022

@author: tschakel
"""

from wad_qc.modulelibs import wadwrapper_lib
import numpy as np
import pydicom
import matplotlib.pyplot as plt
from matplotlib.patches import Circle

from MR_B0B1_util import (find_center, create_circular_mask, create_spherical_mask,
                          calc_rms_stats, calc_percentile_stats, b0_create_figure1, 
                          b0_create_figure2, b0_create_figure3, b0_collect_results)

### Helper functions
def getValue(ds, label):
    """Return the value of a pydicom DataElement in Dataset identified by label.

    ds: pydicom Dataset
    label: dicom identifier, in either pydicom Tag object, string or tuple form.
    """
    if isinstance(label, str):
        try:
            # Assume form "0x0008,0x1030"
            tag = pydicom.tag.Tag(label.split(','))
        except ValueError:
            try:
                # Assume form "SeriesDescription"
                tag = ds.data_element(label).tag
            except (AttributeError, KeyError):
                # `label` string doesn't represent an element of the DataSet
                return None
    else:
        # Assume label is of form (0x0008,0x1030) or is a pydicom Tag object.
        tag = pydicom.tag.Tag(label)

    try:
        return str(ds[tag].value)
    except KeyError:
        # Tag doesn't exist in the DataSet
        return None


def isFiltered(ds, filters):
    """Return True if the Dataset `ds` complies to the `filters`,
    otherwise return False.
    """
    for tag, value in filters.items():
        if not str(getValue(ds, tag)) == str(value):
            # Convert both values to string before comparison. Reason is that
            # pydicom can return 'str', 'int' or 'dicom.valuerep' types of data.
            # Similarly, the user (or XML) supplied value can be of any type.
            return False
    return True


def applyFilters(series_filelist, filters):
    """Apply `filters` to the `series_filelist` and return the filtered list.

    First, convert `filters` from an ElementTree Element to a dictionary
    Next, create a new list in the same shape as `series_filelist`, but only
    include filenames for which isFiltered returns True.
    Only include sublists (i.e., series) which are non empty.
    """
    # Turn ElementTree element attributes and text into filters
    #filter_dict = {element.attrib["name"]: element.text for element in filters}
    filter_dict = filters

    filtered_series_filelist = []
    # For each series in the series_filelist (or, study):
    for instance_filelist in series_filelist:
        # Filter filenames within each series
        filtered_instance_filelist = [fn for fn in instance_filelist
                                      if isFiltered(
                pydicom.read_file(fn, stop_before_pixels=True), filter_dict)]
        # Only add the series which are not empty
        if filtered_instance_filelist:
            filtered_series_filelist.append(filtered_instance_filelist)

    return filtered_series_filelist

def acqdatetime(data, results, action):
    """
    Get the date and time of acquisition
    """
    params = action["params"]
    datetime_series = data.getSeriesByDescription(params["datetime_series_description"])
    dt = wadwrapper_lib.acqdatetime_series(datetime_series[0][0])
    results.addDateTime('AcquisitionDateTime', dt) 

def B1_tra_AFI(data, results, action):
    print(">>> B1 TRA AFI <<<")
    params = action["params"]
    filters = action["filters"]
    filename = 'B1 TRA AFI.png'
    
    b1series_filter = {"SeriesDescription":filters.get(item)for item in ["b1_series_description"]}
    b1series = applyFilters(data.series_filelist, b1series_filter)
    
    type_B1_filter = {"ImageType":filters.get(item)for item in ["b1_imageType"]}
    b1map_series = applyFilters(b1series, type_B1_filter)
    dcmInfile,pixeldataIn,dicomMode = wadwrapper_lib.prepareInput(b1map_series[0],headers_only=False)
    b1map = pixeldataIn[int(params['slicenumber'])-1,:,:]
    
    type_M_filter = {"ImageType":filters.get(item)for item in ["M_image_type"]}
    m_series = applyFilters(b1series, type_M_filter)
    dcmInfileM,pixeldataInM,dicomModeM = wadwrapper_lib.prepareInput(m_series[0],headers_only=False)
    m_image = pixeldataInM[int(params['slicenumber'])-1,:,:]
    
    # calc some stats
    x_center_px, y_center_px = find_center(m_image,params)
    circ_rad = 150 / dcmInfile.info.PixelSpacing[0] # 15 cm radius
    circ_mask = create_circular_mask(b1map,(np.int0(x_center_px),np.int0(y_center_px)),circ_rad)
    b1_circ_mean = b1map[circ_mask].mean()
    b1_circ_std = b1map[circ_mask].std()
    
    # create figures
    fig, axs = plt.subplots(2,2,figsize=(12,6))
    fig.suptitle('B1 TRA slice'+params['slicenumber'])
    major_ticks_plotY = np.arange(90,111,5)
    minor_ticks_plotY = np.arange(90,111,1)
    major_ticks_plotX = np.arange(0,b1map.shape[0]+1,int(np.round(b1map.shape[0]/4)))
    
    axs[0,0].imshow(m_image,cmap='gray')
    axs[0,0].set_title('Magnitude')
    axs[0,0].axis('off')
    
    im = axs[0,1].imshow(b1map,vmin=90,vmax=110)
    axs[0,1].set_title('B1 map [%]')
    plt.colorbar(im, ax=axs[0,1], shrink=0.8)
    axs[0,1].axis('off')
    axs[0,1].axhline(x_center_px, c='r')
    axs[0,1].axvline(y_center_px, c='k')
    axs[0,1].add_patch(Circle([x_center_px,y_center_px],circ_rad,fc='magenta',lw=2,ec='magenta',alpha=0.3))
    
    axs[1,0].plot(b1map[np.int0(x_center_px),:], c='r')
    axs[1,0].set_title('Midline')
    axs[1,0].set_ylabel('B1 [%]')
    axs[1,0].set_ylim([90,110])
    axs[1,0].set_xticks(major_ticks_plotX)
    axs[1,0].set_yticks(major_ticks_plotY)
    axs[1,0].set_yticks(minor_ticks_plotY, minor=True)
    axs[1,0].grid(which='minor', alpha=0.2,axis='y')
    axs[1,0].grid(which='major', alpha=0.8,axis='both')
    
    axs[1,1].plot(b1map[:,np.int0(y_center_px)], c='k')
    axs[1,1].set_title('Midline')
    axs[1,1].set_ylim([90,110])
    axs[1,1].set_xticks(major_ticks_plotX)
    axs[1,1].set_yticks(major_ticks_plotY)
    axs[1,1].set_yticks(minor_ticks_plotY, minor=True)
    axs[1,1].grid(which='minor', alpha=0.2,axis='y')
    axs[1,1].grid(which='major', alpha=0.8,axis='both')
    
    fig.savefig(filename,dpi=150)
    results.addFloat("B1 mean", b1_circ_mean)
    results.addFloat("B1 std", b1_circ_std)
    results.addObject("B1 figure", filename)
    
def B1_tra(data, results, action):
    print(">>> B1 TRA <<<")
    params = action["params"]
    filters = action["filters"]
    filename = 'B1 TRA.png'
    
    b1_60_series_filter = {"SeriesDescription":filters.get(item)for item in ["b1_60_series_description"]}
    b1_60_series = applyFilters(data.series_filelist, b1_60_series_filter)
    
    if len(b1_60_series) < 1:
        print("ERROR: B1 60 series not present")
        return
    
    dcmInfile,pixeldataIn,dicomMode = wadwrapper_lib.prepareInput(b1_60_series[0],headers_only=False)
    if dcmInfile.info.FlipAngle == 60.0:
        b1_60_image_data = pixeldataIn[int(params['slicenumber'])-1,:,:]
    else:
        raise Exception('Series for 60 does not have flip angle 60')
    
    b1_120_series_filter = {"SeriesDescription":filters.get(item)for item in ["b1_120_series_description"]}
    b1_120_series = applyFilters(data.series_filelist, b1_120_series_filter)
    
    if len(b1_120_series) < 1:
        print("ERROR: B1 120 series not present")
        return
    
    dcmInfile,pixeldataIn,dicomMode = wadwrapper_lib.prepareInput(b1_120_series[0],headers_only=False)
    if dcmInfile.info.FlipAngle == 120.0:
        b1_120_image_data = pixeldataIn[int(params['slicenumber'])-1,:,:]
    else:
        raise Exception('Series for 120 does not have flip angle 120')
    
    # filter background
    b1_60_image_data[ b1_60_image_data < 100 ] = 0
    b1_120_image_data[ b1_120_image_data < 100 ] = 0
    
    # calc b1map
    tmp = np.divide(b1_120_image_data, (2 * b1_60_image_data), out=np.zeros_like(b1_120_image_data), where=b1_60_image_data!=0)
    #tmp = b1_120_image_data / (2 * b1_60_image_data)
    b1map = 100 * np.degrees(np.arccos(tmp, out=np.zeros_like(b1_120_image_data), where=tmp!=0)) / 60.0
    #b1map = 100 * np.degrees(np.arccos(tmp)) / 60.000
    b1map[np.isnan(b1map)] = 0
    
    # calc some stats
    x_center_px, y_center_px = find_center(b1_60_image_data,params)
    circ_rad = 150 / dcmInfile.info.PixelSpacing[0] # 15 cm radius
    circ_mask = create_circular_mask(b1map,(np.int0(x_center_px),np.int0(y_center_px)),circ_rad)
    b1_circ_mean = b1map[circ_mask].mean()
    b1_circ_std = b1map[circ_mask].std()
    
    # create figures
    fig, axs = plt.subplots(2,3,figsize=(12,6))
    fig.suptitle('B1 TRA')
    major_ticks_plotY = np.arange(90,111,5)
    minor_ticks_plotY = np.arange(90,111,1)
    major_ticks_plotX = np.arange(0,b1map.shape[0]+1,int(np.round(b1map.shape[0]/4)))
    
    axs[0,0].imshow(b1_60_image_data,cmap='gray')
    axs[0,0].set_title('FlipAngle 60')
    axs[0,0].axis('off')
    
    axs[0,1].imshow(b1_120_image_data,cmap='gray')
    axs[0,1].set_title('FlipAngle 120')
    axs[0,1].axis('off')
    
    im = axs[0,2].imshow(b1map,vmin=90,vmax=110)
    axs[0,2].set_title('B1 map [%]')
    plt.colorbar(im, ax=axs[0,2], shrink=0.8)
    axs[0,2].axis('off')
    
    axs[1,0].imshow(b1map,vmin=90,vmax=110)
    axs[1,0].set_title('B1 map [%]')
    axs[1,0].axis('off')
    axs[1,0].axhline(x_center_px, c='r')
    axs[1,0].axvline(y_center_px, c='k')
    axs[1,0].add_patch(Circle([x_center_px,y_center_px],circ_rad,fc='magenta',lw=2,ec='magenta',alpha=0.3))
    
    axs[1,1].plot(b1map[np.int0(x_center_px),:], c='r')
    axs[1,1].set_title('Midline')
    axs[1,1].set_ylabel('B1 [%]')
    axs[1,1].set_ylim([90,110])
    axs[1,1].set_xticks(major_ticks_plotX)
    axs[1,1].set_yticks(major_ticks_plotY)
    axs[1,1].set_yticks(minor_ticks_plotY, minor=True)
    axs[1,1].grid(which='minor', alpha=0.2,axis='y')
    axs[1,1].grid(which='major', alpha=0.8,axis='both')
    
    axs[1,2].plot(b1map[:,np.int0(y_center_px)], c='k')
    axs[1,2].set_title('Midline')
    axs[1,2].set_ylim([90,110])
    axs[1,2].set_xticks(major_ticks_plotX)
    axs[1,2].set_yticks(major_ticks_plotY)
    axs[1,2].set_yticks(minor_ticks_plotY, minor=True)
    axs[1,2].grid(which='minor', alpha=0.2,axis='y')
    axs[1,2].grid(which='major', alpha=0.8,axis='both')
    
    fig.savefig(filename,dpi=300)
    results.addFloat("B1 mean", b1_circ_mean)
    results.addFloat("B1 std", b1_circ_std)
    results.addObject("B1 figure", filename)
    
def B0_tra(data, results, action):
    print(">>> B0 TRA <<<")
    params = action["params"]
    filters = action["filters"]
    basename = 'B0 '
    
    b0_series_filter = {"SeriesDescription":filters.get(item)for item in ["b0_series_description"]}
    type_b0_filter = {"ImageType":filters.get(item)for item in ["b0_image_type"]}
    type_M_filter = {"ImageType":filters.get(item)for item in ["M_image_type"]}
    data_series = applyFilters(data.series_filelist, b0_series_filter)
    data_series_b0 = applyFilters(data_series, type_b0_filter)
    data_series_M = applyFilters(data_series, type_M_filter)
    
    dcmInfile,b0map,dicomMode = wadwrapper_lib.prepareInput(data_series_b0[0],headers_only=False)
    dcmInfileM,pixeldataIn,dicomMode = wadwrapper_lib.prepareInput(data_series_M[0],headers_only=False)
    
    b0map = np.transpose(b0map,(1,2,0)) #transpose to x,y,z
    pixeldataIn = np.transpose(pixeldataIn,(1,2,0))
    
    phantom = pixeldataIn > 200
    imaging_frequency = dcmInfile.info.ImagingFrequency
    b0map_ppm = b0map / imaging_frequency
    
    # get the orientation
    seriesname = dcmInfile.info.SeriesDescription
    if seriesname.startswith('t'):
        scanori = 'TRANSVERSAL'
    elif seriesname.startswith('c'):
        scanori = 'CORONAL'
    elif seriesname.startswith('s'):
        scanori = 'SAGITTAL'
    else:
        scanori = 'UNKNOWN'

    # slicenumber
    slicenumber = int(params['slicenumber'])
    
    # define volumes for statistics
    # Load dimensions based on the number of rows, columns, and slices (along the Z axis)
    pixelDims = (int(dcmInfile.info.Rows), int(dcmInfile.info.Columns), dcmInfile.info[0x20011018].value)
    
    # Load spacing values (in mm)
    gap = float(dcmInfile.info.SpacingBetweenSlices)
    pixelSpacing = (float(dcmInfile.info.PixelSpacing[0]), float(dcmInfile.info.PixelSpacing[1]),float(dcmInfile.info.SliceThickness)+gap)
    
    # find pixel coords of geometric center (0,0)
    # ImagePositionPatient gives the x,y,z coordinates of the center of the pixel in the upper lefthand corner
    x0 = np.round((0 - dcmInfile.info.ImagePositionPatient[0] + (pixelSpacing[0] / 2) ) / pixelSpacing[0])
    y0 = np.round((0 - dcmInfile.info.ImagePositionPatient[1] + (pixelSpacing[1] / 2) ) / pixelSpacing[1])
    z0 = np.round((0 - dcmInfile.info.ImagePositionPatient[2] + (pixelSpacing[2] / 2) ) / pixelSpacing[2])

    # Create DSV masks
    dsv100 = create_spherical_mask(b0map_ppm.shape,voxel_spacing=(pixelSpacing[0],pixelSpacing[1],pixelSpacing[2]), center=(x0,y0,z0), radius=50)
    dsv200 = create_spherical_mask(b0map_ppm.shape,voxel_spacing=(pixelSpacing[0],pixelSpacing[1],pixelSpacing[2]), center=(x0,y0,z0), radius=100)
    dsv300 = create_spherical_mask(b0map_ppm.shape,voxel_spacing=(pixelSpacing[0],pixelSpacing[1],pixelSpacing[2]), center=(x0,y0,z0), radius=150)
    dsv350 = create_spherical_mask(b0map_ppm.shape,voxel_spacing=(pixelSpacing[0],pixelSpacing[1],pixelSpacing[2]), center=(x0,y0,z0), radius=175)
    
    # calculate rms statistics within the different masks
    rms_dsv100 = calc_rms_stats(b0map_ppm,dsv100)
    rms_dsv200 = calc_rms_stats(b0map_ppm,dsv200)
    rms_dsv300 = calc_rms_stats(b0map_ppm,dsv300)
    rms_dsv350 = calc_rms_stats(b0map_ppm,dsv350)
    rms_phantom = calc_rms_stats(b0map_ppm, phantom)
    
    # calculate different percentile statistics for each slice
    perc_stats = calc_percentile_stats(b0map_ppm,phantom)
    
    acqdate = dcmInfile.info.StudyDate
    acqtime = dcmInfile.info.StudyTime
    fname_fig1 = b0_create_figure1(b0map, phantom, acqdate, acqtime, scanori,basename)
    fname_fig2 = b0_create_figure2(b0map_ppm,b0map,phantom,dsv100,dsv200,dsv300,
                                   dsv350,acqdate,acqtime,imaging_frequency, scanori,basename,slicenumber)
    fname_fig3 = b0_create_figure3(b0map_ppm, phantom, acqdate, acqtime, scanori,basename)
    
    # collect results
    results = b0_collect_results(results,scanori,seriesname,
                                 fname_fig1,fname_fig2,fname_fig3,
                                 rms_dsv100,rms_dsv200,rms_dsv300,rms_dsv350,
                                 perc_stats,rms_phantom,pixelDims,basename)
    
def B0_tra_noshim(data, results, action):
    print(">>> B0 TRA NOSHIM <<<")
    params = action["params"]
    filters = action["filters"]
    basename = 'B0 NOSHIM '
    
    b0_series_filter = {"SeriesDescription":filters.get(item)for item in ["b0_series_description"]}
    type_b0_filter = {"ImageType":filters.get(item)for item in ["b0_image_type"]}
    type_M_filter = {"ImageType":filters.get(item)for item in ["M_image_type"]}
    data_series = applyFilters(data.series_filelist, b0_series_filter)
    
    # sometimes, the no_shim is not acquired: applyFilters will return an empyt list to data_series
    # catch this and stop further analysis
    if len(data_series) < 1:
        print("ERROR: NO SHIM series not present")
        return
    
    data_series_b0 = applyFilters(data_series, type_b0_filter)
    data_series_M = applyFilters(data_series, type_M_filter)
    
    dcmInfile,b0map,dicomMode = wadwrapper_lib.prepareInput(data_series_b0[0],headers_only=False)
    dcmInfileM,pixeldataIn,dicomMode = wadwrapper_lib.prepareInput(data_series_M[0],headers_only=False)
    
    b0map = np.transpose(b0map,(1,2,0)) #transpose to x,y,z
    pixeldataIn = np.transpose(pixeldataIn,(1,2,0))
    
    phantom = pixeldataIn > 200
    imaging_frequency = dcmInfile.info.ImagingFrequency
    b0map_ppm = b0map / imaging_frequency
    
    # get the orientation
    seriesname = dcmInfile.info.SeriesDescription
    if seriesname.startswith('t'):
        scanori = 'TRANSVERSAL'
    elif seriesname.startswith('c'):
        scanori = 'CORONAL'
    elif seriesname.startswith('s'):
        scanori = 'SAGITTAL'
    else:
        scanori = 'UNKNOWN'

    # slicenumber
    slicenumber = int(params['slicenumber'])
    
    # define volumes for statistics
    # Load dimensions based on the number of rows, columns, and slices (along the Z axis)
    pixelDims = (int(dcmInfile.info.Rows), int(dcmInfile.info.Columns), dcmInfile.info[0x20011018].value)
    
    # Load spacing values (in mm)
    gap = float(dcmInfile.info.SpacingBetweenSlices)
    pixelSpacing = (float(dcmInfile.info.PixelSpacing[0]), float(dcmInfile.info.PixelSpacing[1]),float(dcmInfile.info.SliceThickness)+gap)
    
    # find pixel coords of geometric center (0,0)
    # ImagePositionPatient gives the x,y,z coordinates of the center of the pixel in the upper lefthand corner
    x0 = np.round((0 - dcmInfile.info.ImagePositionPatient[0] + (pixelSpacing[0] / 2) ) / pixelSpacing[0])
    y0 = np.round((0 - dcmInfile.info.ImagePositionPatient[1] + (pixelSpacing[1] / 2) ) / pixelSpacing[1])
    z0 = np.round((0 - dcmInfile.info.ImagePositionPatient[2] + (pixelSpacing[2] / 2) ) / pixelSpacing[2])

    # Create DSV masks
    dsv100 = create_spherical_mask(b0map_ppm.shape,voxel_spacing=(pixelSpacing[0],pixelSpacing[1],pixelSpacing[2]), center=(x0,y0,z0), radius=50)
    dsv200 = create_spherical_mask(b0map_ppm.shape,voxel_spacing=(pixelSpacing[0],pixelSpacing[1],pixelSpacing[2]), center=(x0,y0,z0), radius=100)
    dsv300 = create_spherical_mask(b0map_ppm.shape,voxel_spacing=(pixelSpacing[0],pixelSpacing[1],pixelSpacing[2]), center=(x0,y0,z0), radius=150)
    dsv350 = create_spherical_mask(b0map_ppm.shape,voxel_spacing=(pixelSpacing[0],pixelSpacing[1],pixelSpacing[2]), center=(x0,y0,z0), radius=175)
    
    # calculate rms statistics within the different masks
    rms_dsv100 = calc_rms_stats(b0map_ppm,dsv100)
    rms_dsv200 = calc_rms_stats(b0map_ppm,dsv200)
    rms_dsv300 = calc_rms_stats(b0map_ppm,dsv300)
    rms_dsv350 = calc_rms_stats(b0map_ppm,dsv350)
    rms_phantom = calc_rms_stats(b0map_ppm, phantom)
    
    # calculate different percentile statistics for each slice
    perc_stats = calc_percentile_stats(b0map_ppm,phantom)
    
    acqdate = dcmInfile.info.StudyDate
    acqtime = dcmInfile.info.StudyTime
    fname_fig1 = b0_create_figure1(b0map, phantom, acqdate, acqtime, scanori,basename)
    fname_fig2 = b0_create_figure2(b0map_ppm,b0map,phantom,dsv100,dsv200,dsv300,
                                   dsv350,acqdate,acqtime,imaging_frequency, scanori,basename,slicenumber)
    fname_fig3 = b0_create_figure3(b0map_ppm, phantom, acqdate, acqtime, scanori,basename)
    
    # collect results
    results = b0_collect_results(results,scanori,seriesname,
                                 fname_fig1,fname_fig2,fname_fig3,
                                 rms_dsv100,rms_dsv200,rms_dsv300,rms_dsv350,
                                 perc_stats,rms_phantom,pixelDims,basename)
    