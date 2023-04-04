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
from skimage.restoration import unwrap_phase
import MR_B0B1_dcm_input
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

def isFilteredContains(ds, filters):
    """Return True if the Dataset `ds` complies to the `filters`,
    otherwise return False.
    """
    for tag, value in filters.items():
        if str(value) not in str(getValue(ds, tag)):
        #if not str(getValue(ds, tag)).startswith(str(value)):
        #if not str(getValue(ds, tag)) == str(value):
            # Convert both values to string before comparison. Reason is that
            # pydicom can return 'str', 'int' or 'dicom.valuerep' types of data.
            # Similarly, the user (or XML) supplied value can be of any type.
            return False
    return True

def flatten(nestedlist):
    """Flatten a list of lists"""
    return [item for sublist in nestedlist for item in sublist]

def getInstanceByTagsPartial(data, filters):
    """Return a list of dicom instances which satisfy the specified filters.

    filters: dictionary where each key, value pair is a pydicom DataElement
        name or tag number and the corresponding filter

    Example:
    myfilter = {
        "ExposureTime": "100",
        "PatientName": "WAD",
        (0x0018,0x1405): "13"
    }
    """
    func = lambda fn: pydicom.read_file(fn)
    
    instances = (func(fn) for fn in flatten(data.series_filelist))
    return [ds for ds in instances if isFilteredContains(ds, filters)]

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
    filters = action["filters"]
    datetime_series = data.getSeriesByDescription(params["datetime_series_description"])
    
    # Add an exception for when the user does not follow the instructions
    # Try some partial matching on the scan names
    if len(datetime_series) < 1:
        datetime_series = getInstanceByTagsPartial(data,filters["datetime_filter_partial"])
        dt = wadwrapper_lib.acqdatetime_series(datetime_series[0])
    else:
       dt = wadwrapper_lib.acqdatetime_series(datetime_series[0][0])
       
    results.addDateTime('AcquisitionDateTime', dt) 

def B1_AFI(data, results, action):
    # Analyse B1 acquired with AFI
    # The scanner calculates B1 maps, here just read and get stats
    print(">>> B1 TRA AFI <<<")
    params = action["params"]
    filters = action["filters"]
    
    b1series_filter = {"SeriesDescription":filters.get(item)for item in ["b1_series_description"]}
    b1series = applyFilters(data.series_filelist, b1series_filter)  
    if len(b1series) < 1:
        print(">>> B1_AFI_50_150 TRA not found <<<")
        return
    
    type_B1_filter = {"ImageType":filters.get(item)for item in ["b1_imageType"]}
    b1map_series = applyFilters(b1series, type_B1_filter)
    if len(b1map_series) < 1:
        print("ERROR: t B1_AFI_50_150 series does not contain B1 map")
        return
    
    #dcmInfile,pixeldataIn,dicomMode = wadwrapper_lib.prepareInput(b1map_series[0],headers_only=False)
    dcmInfile,pixeldataIn,dicomMode = MR_B0B1_dcm_input.prepareInput(b1map_series[0],headers_only=False)
    b1map = pixeldataIn[int(params['slicenumber'])-1,:,:]
    
    type_M_filter = {"ImageType":filters.get(item)for item in ["M_image_type"]}
    m_series = applyFilters(b1series, type_M_filter)
    if len(b1map_series) < 1:
        print("ERROR: t B1_AFI_50_150 series does not contain magnitude image")
        return
    
    #dcmInfileM,pixeldataInM,dicomModeM = wadwrapper_lib.prepareInput(m_series[0],headers_only=False)
    dcmInfileM,pixeldataInM,dicomModeM = MR_B0B1_dcm_input.prepareInput(m_series[0],headers_only=False)
    m_image = pixeldataInM[int(params['slicenumber'])-1,:,:]
    
    # Define mask (DSV) for analysis
    # Load dimensions based on the number of rows, columns
    pixelDims = (int(dcmInfile.info.Rows), int(dcmInfile.info.Columns))
    
    # Load spacing values (in mm)
    pixelSpacing = (float(dcmInfile.info.PixelSpacing[0]), float(dcmInfile.info.PixelSpacing[1]))
    
    # Create mask
    radius = 175
    dsv350 = create_circular_mask(b1map.shape, voxel_spacing=(pixelSpacing[0],pixelSpacing[1]), 
                                  center=(pixelDims[0]/2,pixelDims[1]/2), radius=radius)
    
    b1map_masked = b1map * dsv350
    b1map_masked = b1map_masked[np.nonzero(b1map_masked)]
    
    # Get main results
    mean_dsv350 = np.mean(b1map_masked)
    std_dsv350 = np.std(b1map_masked)
    range_dsv350 = np.percentile(b1map_masked, 99) - np.percentile(b1map_masked, 1)
    
    # create figures
    figtitle = 'B1 AFI [%] TRA '+dcmInfile.info.StudyDate+' '+dcmInfile.info.StudyTime
    fig, axs = plt.subplots(2,2,figsize=(12,6))
    fig.suptitle(figtitle)
    major_ticks_plotY = np.arange(80,111,5)
    minor_ticks_plotY = np.arange(80,111,1)
    major_ticks_plotX = np.arange(0,b1map.shape[0]+1,int(np.round(b1map.shape[0]/4)))
    
    axs[0,0].imshow(m_image,cmap='gray')
    axs[0,0].set_title('Magnitude')
    axs[0,0].axis('off')
    
    im = axs[0,1].imshow(b1map,vmin=80,vmax=110)
    axs[0,1].set_title('B1 map [%]')
    plt.colorbar(im, ax=axs[0,1], shrink=0.8)
    axs[0,1].axis('off')
    axs[0,1].axhline(pixelDims[0]/2, c='r')
    axs[0,1].axvline(pixelDims[1]/2, c='k')
    axs[0,1].add_patch(Circle([pixelDims[0]/2,pixelDims[1]/2],radius/dcmInfile.info.PixelSpacing[0],fc='magenta',lw=2,ec='magenta',alpha=0.3))
    
    axs[1,0].plot(b1map[np.int0(pixelDims[0]/2),:], c='r')
    axs[1,0].set_title('Midline')
    axs[1,0].set_ylabel('B1 [%]')
    axs[1,0].set_ylim([80,110])
    axs[1,0].set_xticks(major_ticks_plotX)
    axs[1,0].set_yticks(major_ticks_plotY)
    axs[1,0].set_yticks(minor_ticks_plotY, minor=True)
    axs[1,0].grid(which='minor', alpha=0.2,axis='y')
    axs[1,0].grid(which='major', alpha=0.8,axis='both')
    
    axs[1,1].plot(b1map[:,np.int0(pixelDims[1]/2)], c='k')
    axs[1,1].set_title('Midline')
    axs[1,1].set_ylim([80,110])
    axs[1,1].set_xticks(major_ticks_plotX)
    axs[1,1].set_yticks(major_ticks_plotY)
    axs[1,1].set_yticks(minor_ticks_plotY, minor=True)
    axs[1,1].grid(which='minor', alpha=0.2,axis='y')
    axs[1,1].grid(which='major', alpha=0.8,axis='both')
    
    filename = 'B1_AFI_TRA.png'
    fig.savefig(filename,dpi=150)
    
    # Save results
    results.addFloat("B1_AFI_TRA_mean", mean_dsv350)
    results.addFloat("B1_AFI_TRA_std", std_dsv350)
    results.addFloat("B1_AFI_TRA_p99-p1", range_dsv350)
    results.addObject("B1_AFI_TRA_figure", filename)
    

def process_B1(data_series_b1_60,data_series_b1_120,scanori,params,results):
    # Perform B1 analysis 
    slicenumber = int(params['slicenumber'])
    
    # Read the data and undo the scaling
    dcmInfile,pixeldataIn,dicomMode = MR_B0B1_dcm_input.prepareInput(data_series_b1_60[0],headers_only=False)
    if not "RescaleSlope" in dcmInfile.info:
        # if slope is not present, private tags
        # look at the scaling in the reader function. there are differences between enhanced and normal
        # due to changes in the dicomformat sent from the scanners, these exceptions for missing dicomtags are required
        # should only be needed for the historic data.
        print("Getting Slope/Intercept from private tags")
        dcmInfile.info.RescaleSlope = dcmInfile.info[0x2005,0x100e].value
        dcmInfile.info.RescaleIntercept = dcmInfile.info[0x2005,0x100d].value

    slope60 = dcmInfile.info.RescaleSlope
    inter60 = dcmInfile.info.RescaleIntercept
    pixeldataIn = np.transpose(pixeldataIn,(1,2,0)) 
    b1_60_image_data = (pixeldataIn[:,:,slicenumber] - inter60) / slope60
            
    dcmInfile,pixeldataIn,dicomMode = MR_B0B1_dcm_input.prepareInput(data_series_b1_120[0],headers_only=False)
    if not "RescaleSlope" in dcmInfile.info:
        # if slope is not present, private tags
        print("Getting Slope/Intercept from private tags")
        dcmInfile.info.RescaleSlope = dcmInfile.info[0x2005,0x100e].value
        dcmInfile.info.RescaleIntercept = dcmInfile.info[0x2005,0x100d].value
        
    slope120 = dcmInfile.info.RescaleSlope
    inter120 = dcmInfile.info.RescaleIntercept
    pixeldataIn = np.transpose(pixeldataIn,(1,2,0)) 
    b1_120_image_data = (pixeldataIn[:,:,slicenumber] - inter120) / slope120
    
    
    # calc b1map
    # return result of calculation; eps is added for numerical stability
    eps = 1e-6
    b1map = 100 * np.degrees(np.arccos(b1_120_image_data/(2 * b1_60_image_data+eps), out=np.zeros_like(b1_120_image_data), where=b1_120_image_data!=0)) / 60.0
    b1map[np.isnan(b1map)] = 0
        
    # Define mask (DSV) for analysis
    # Load dimensions based on the number of rows, columns
    pixelDims = (int(dcmInfile.info.Rows), int(dcmInfile.info.Columns))
    
    # Load spacing values (in mm)
    pixelSpacing = (float(dcmInfile.info.PixelSpacing[0]), float(dcmInfile.info.PixelSpacing[1]))
    
    # Create mask
    radius = 175
    dsv350 = create_circular_mask(b1map.shape, voxel_spacing=(pixelSpacing[0],pixelSpacing[1]), 
                                  center=(pixelDims[0]/2,pixelDims[1]/2), radius=radius)
    
    b1map_masked = b1map * dsv350
    b1map_masked = b1map_masked[np.nonzero(b1map_masked)]
    
    # Get main results
    mean_dsv350 = np.mean(b1map_masked)
    std_dsv350 = np.std(b1map_masked)
    range_dsv350 = np.percentile(b1map_masked, 99) - np.percentile(b1map_masked, 1)
    
    # create figures
    figtitle = 'B1 [%] '+scanori+' '+dcmInfile.info.StudyDate+' '+dcmInfile.info.StudyTime
    fig, axs = plt.subplots(2,3,figsize=(12,6))
    fig.suptitle(figtitle)
    major_ticks_plotY = np.arange(80,111,5)
    minor_ticks_plotY = np.arange(80,111,1)
    major_ticks_plotX = np.arange(0,b1map.shape[0]+1,int(np.round(b1map.shape[0]/4)))
    
    axs[0,0].imshow(b1_60_image_data,cmap='gray')
    axs[0,0].set_title('FlipAngle 60')
    axs[0,0].axis('off')
    
    axs[0,1].imshow(b1_120_image_data,cmap='gray')
    axs[0,1].set_title('FlipAngle 120')
    axs[0,1].axis('off')
    
    im = axs[0,2].imshow(b1map,vmin=80,vmax=110)
    axs[0,2].set_title('B1 map [%]')
    plt.colorbar(im, ax=axs[0,2], shrink=0.8)
    axs[0,2].axis('off')
    
    axs[1,0].imshow(b1map,vmin=80,vmax=110)
    axs[1,0].set_title('B1 map [%]')
    axs[1,0].axis('off')
    axs[1,0].axhline(pixelDims[0]/2, c='r')
    axs[1,0].axvline(pixelDims[1]/2, c='k')
    axs[1,0].add_patch(Circle([pixelDims[0]/2,pixelDims[1]/2],radius/dcmInfile.info.PixelSpacing[0],fc='magenta',lw=2,ec='magenta',alpha=0.3))
    
    axs[1,1].plot(b1map[np.int0(pixelDims[0]/2),:], c='r')
    axs[1,1].set_title('Midline')
    axs[1,1].set_ylabel('B1 [%]')
    axs[1,1].set_ylim([80,110])
    axs[1,1].set_xticks(major_ticks_plotX)
    axs[1,1].set_yticks(major_ticks_plotY)
    axs[1,1].set_yticks(minor_ticks_plotY, minor=True)
    axs[1,1].grid(which='minor', alpha=0.2,axis='y')
    axs[1,1].grid(which='major', alpha=0.8,axis='both')
    
    axs[1,2].plot(b1map[:,np.int0(pixelDims[1]/2)], c='k')
    axs[1,2].set_title('Midline')
    axs[1,2].set_ylim([80,110])
    axs[1,2].set_xticks(major_ticks_plotX)
    axs[1,2].set_yticks(major_ticks_plotY)
    axs[1,2].set_yticks(minor_ticks_plotY, minor=True)
    axs[1,2].grid(which='minor', alpha=0.2,axis='y')
    axs[1,2].grid(which='major', alpha=0.8,axis='both')
    
    filename = 'B1_'+scanori+'.png'
    fig.savefig(filename,dpi=150)
    
    # Save results
    results.addFloat("B1_"+scanori+"_mean", mean_dsv350)
    results.addFloat("B1_"+scanori+"_std", std_dsv350)
    results.addFloat("B1_"+scanori+"_p99-p1", range_dsv350)
    results.addObject("B1_"+scanori+"_figure", filename)
        
    return results
  
def B1_60_120(data, results, action):
    # Run the B1 analysis for the 60-120 method for all different orientations
    print(">>> B1 60/120<<<")
    params = action["params"]
    filters = action["filters"]

    # Filter for transverse
    b1_tra_60_filter = {"SeriesDescription":filters.get(item)for item in ["b1_tra_60_series_description"]}
    data_series_b1tra60 = applyFilters(data.series_filelist,b1_tra_60_filter)
    b1_tra_120_filter = {"SeriesDescription":filters.get(item)for item in ["b1_tra_120_series_description"]}
    data_series_b1tra120 = applyFilters(data.series_filelist,b1_tra_120_filter)
    if len(data_series_b1tra60) < 1 or len(data_series_b1tra120) < 1:
        print(">>> B1 TRA 60 or 120 not found <<<")
    else:
        scanori = 'TRA'
        results = process_B1(data_series_b1tra60,data_series_b1tra120,scanori,params,results)
        
    # Filter for coronal
    b1_cor_60_filter = {"SeriesDescription":filters.get(item)for item in ["b1_cor_60_series_description"]}
    data_series_b1cor60 = applyFilters(data.series_filelist,b1_cor_60_filter)
    b1_cor_120_filter = {"SeriesDescription":filters.get(item)for item in ["b1_cor_120_series_description"]}
    data_series_b1cor120 = applyFilters(data.series_filelist,b1_cor_120_filter)
    if len(data_series_b1cor60) < 1 or len(data_series_b1cor120) < 1:
        print(">>> B1 COR 60 or 120 not found <<<")
    else:
        scanori = 'COR'
        results = process_B1(data_series_b1cor60,data_series_b1cor120,scanori,params,results)
        
    # Filter for sagittal
    b1_sag_60_filter = {"SeriesDescription":filters.get(item)for item in ["b1_sag_60_series_description"]}
    data_series_b1sag60 = applyFilters(data.series_filelist,b1_sag_60_filter)
    b1_sag_120_filter = {"SeriesDescription":filters.get(item)for item in ["b1_sag_120_series_description"]}
    data_series_b1sag120 = applyFilters(data.series_filelist,b1_sag_120_filter)
    if len(data_series_b1sag60) < 1 or len(data_series_b1sag120) < 1:
        print(">>> B1 SAG 60 or 120 not found <<<")
    else:
        scanori = 'SAG'
        results = process_B1(data_series_b1sag60,data_series_b1sag120,scanori,params,results)
        
    
def filter_B0(filters,data,b0_series_filter):
    type_b0_filter = {"ImageType":filters.get(item)for item in ["b0_image_type"]}
    data_series = applyFilters(data.series_filelist, b0_series_filter)
    data_series_b0 = applyFilters(data_series, type_b0_filter)
    
    return data_series_b0 

def process_B0(data_series_b0,scanori,shim,params,results):
    
    # Read the data
    dcmInfile,b0map,dicomMode = MR_B0B1_dcm_input.prepareInput(data_series_b0[0],headers_only=False)
    b0map = np.transpose(b0map,(1,2,0)) 

    # Extract the correct slice
    slicenumber = int(params['slicenumber'])
    
    # Convert to PPM
    imaging_frequency = dcmInfile.info.ImagingFrequency
    b0map_ppm = b0map / imaging_frequency
    
    # Define mask (DSV) for analysis
    # Load dimensions based on the number of rows, columns
    pixelDims = (int(dcmInfile.info.Rows), int(dcmInfile.info.Columns))
    
    # Load spacing values (in mm)
    pixelSpacing = (float(dcmInfile.info.PixelSpacing[0]), float(dcmInfile.info.PixelSpacing[1]))
    
    # Create mask
    radius = 175
    dsv350 = create_circular_mask(b0map_ppm.shape[:2], voxel_spacing=(pixelSpacing[0],pixelSpacing[1]), 
                                  center=(pixelDims[0]/2,pixelDims[1]/2), radius=radius)
    
    # apply mask
    for slice in np.arange(b0map_ppm.shape[2]):
        b0map_ppm[~dsv350,slice] = np.nan
    
    # Get main results
    mean_dsv350 = np.nanmean(b0map_ppm)
    std_dsv350 = np.nanstd(b0map_ppm)
    range_dsv350 = np.nanpercentile(b0map_ppm, 99) - np.nanpercentile(b0map_ppm, 1)
    
    # Create figure
    figtitle = 'B0 [ppm] '+shim+' '+scanori+' '+dcmInfile.info.StudyDate+' '+dcmInfile.info.StudyTime
    fig, axs = plt.subplots(1,1,figsize=(7,6))
    pos = axs.imshow(b0map_ppm[:,:,slicenumber],vmin=-5,vmax=5)
    axs.add_patch(Circle([pixelDims[0]/2,pixelDims[1]/2],175/dcmInfile.info.PixelSpacing[0],fc='None',lw=2,ec='red'))  
    fig.suptitle(figtitle,fontsize=20)
    fig.colorbar(pos, ax=axs, label='[ppm]')
    result_str = 'Mean = %.2f\nStDev = %.2f\np99-p1 = %.2f'%(mean_dsv350, std_dsv350,range_dsv350)
    axs.text(0.01, 0.99, result_str, transform=axs.transAxes, fontsize=14,verticalalignment='top')
    filename = 'B0_'+shim+'_'+scanori+'.png'
    fig.savefig(filename,dpi=150)
    
    # Save results
    results.addFloat("B0_"+shim+"_"+scanori+"_mean", mean_dsv350)
    results.addFloat("B0_"+shim+"_"+scanori+"_std", std_dsv350)
    results.addFloat("B0_"+shim+"_"+scanori+"_p99-p1", range_dsv350)
    results.addObject("B0_"+shim+"_"+scanori+"_figure", filename)
        
    return results
    
def B0_shim(data,results,action):
    # Run the B0 analysis for all different orientations with shimming
    print(">>> B0 Shim <<<")
    params = action["params"]
    filters = action["filters"]
    shim = "SHIM"
    
    # Filter for transverse
    b0_tra_filter = {"SeriesDescription":filters.get(item)for item in ["b0_tra_series_description"]}
    data_series_b0tra = filter_B0(filters,data,b0_tra_filter)
    if len(data_series_b0tra) < 1:
        print(">>> B0 Shim TRA not found <<<")
    else:
        scanori = 'TRA'
        results = process_B0(data_series_b0tra,scanori,shim,params,results)
        
    # Filter for coronal
    b0_cor_filter = {"SeriesDescription":filters.get(item)for item in ["b0_cor_series_description"]}
    data_series_b0cor = filter_B0(filters,data,b0_cor_filter)
    if len(data_series_b0cor) < 1:
        print(">>> B0 Shim COR not found <<<")
    else:
        scanori = 'COR'
        results = process_B0(data_series_b0cor,scanori,shim,params,results)
        
    # Filter for sagittal
    b0_sag_filter = {"SeriesDescription":filters.get(item)for item in ["b0_sag_series_description"]}
    data_series_b0sag = filter_B0(filters,data,b0_sag_filter)
    if len(data_series_b0sag) < 1:
        print(">>> B0 Shim SAG not found <<<")
    else:
        scanori = 'SAG'
        results = process_B0(data_series_b0sag,scanori,shim,params,results)
        
def B0_noshim(data,results,action):
    # Run the B0 analysis for all different orientations without shimming
    # (could combine with the shim function)
    print(">>> B0 NO Shim <<<")
    params = action["params"]
    filters = action["filters"]
    shim = "NOSHIM"
    
    # Filter for transverse
    b0_tra_filter = {"SeriesDescription":filters.get(item)for item in ["b0_tra_series_description"]}
    data_series_b0tra = filter_B0(filters,data,b0_tra_filter)
    if len(data_series_b0tra) < 1:
        print(">>> B0 NO Shim TRA not found <<<")
    else:
        scanori = 'TRA'
        results = process_B0(data_series_b0tra,scanori,shim,params,results)
        
    # Filter for coronal
    b0_cor_filter = {"SeriesDescription":filters.get(item)for item in ["b0_cor_series_description"]}
    data_series_b0cor = filter_B0(filters,data,b0_cor_filter)
    if len(data_series_b0cor) < 1:
        print(">>> B0 NO Shim COR not found <<<")
    else:
        scanori = 'COR'
        results = process_B0(data_series_b0cor,scanori,shim,params,results)
        
    # Filter for sagittal
    b0_sag_filter = {"SeriesDescription":filters.get(item)for item in ["b0_sag_series_description"]}
    data_series_b0sag = filter_B0(filters,data,b0_sag_filter)
    if len(data_series_b0sag) < 1:
        print(">>> B0 NO Shim SAG not found <<<")
    else:
        scanori = 'SAG'
        results = process_B0(data_series_b0sag,scanori,shim,params,results)   
        
def B0_gantry(data,results,action):
    print(">>> B0 Gantry <<<")
    params = action["params"]
    filters = action["filters"]
    slicenumber = int(params['slicenumber'])
    tolerance = float(params['p99-p1_tolerance_ppm'])
    
    # Filter for Gantry B0
    b0_gantry_filter = {"SeriesDescription":filters.get(item)for item in ["b0_gantry_series_description"]}
    b0_gantry_data_series = applyFilters(data.series_filelist, b0_gantry_filter)
    
    if len(b0_gantry_data_series) < 1:
        print(">>> B0 Gantry: Series not found <<<")
        return
    
    # Filter for phase images
    p_b0_gantry_filter = {"ImageType":filters.get(item)for item in ["P_image_type"]}
    data_series = applyFilters(b0_gantry_data_series, p_b0_gantry_filter)
    
    dcmInfile,phase_maps,dicomMode = MR_B0B1_dcm_input.prepareInput(data_series[0],headers_only=False)
    phase_maps = np.reshape(phase_maps,(-1,3,256,256),'F') #reshape to [angles,slices,row,col]
    
    gantry_angles=[-180, -150, -120, -90, -60, -30, 0, 30, 60, 90, 120, 150, 180]
    
    # Check number of different datasets found
    if(phase_maps.shape[0] != 13):
        print(">>> B0 Gantry: Incorrect number of datasets found <<<")
        return
    
    # Set the reference phase map
    phase_map_G0 = phase_maps[6]
    
    # Define mask (DSV) for analysis
    # Load dimensions based on the number of rows, columns
    pixelDims = (int(dcmInfile.info.Rows), int(dcmInfile.info.Columns))
    
    # Load spacing values (in mm)
    pixelSpacing = (float(dcmInfile.info.PixelSpacing[0]), float(dcmInfile.info.PixelSpacing[1]))
    
    # Create mask
    radius = 175
    dsv350 = create_circular_mask(phase_map_G0[slicenumber].shape, voxel_spacing=(pixelSpacing[0],pixelSpacing[1]), 
                                  center=(pixelDims[0]/2,pixelDims[1]/2), radius=radius)
    
    range_dsv350 = []
    for phase_map_Ga,angle in zip(phase_maps,gantry_angles):
        
        if angle == 0:
            range_dsv350.append(0.0)
        else:
            #unwrap phase division by 1000 to bring phase (-pi,pi)
            phase_Ga_unwrapped = unwrap_phase(phase_map_Ga[slicenumber]/1000.)
            phase_Gref_unwrapped = unwrap_phase(phase_map_G0[slicenumber]/1000.)
    
            phase_diff=(phase_Ga_unwrapped-phase_Gref_unwrapped)
        
            #correct centre of unwrapping for factors of +/- 2pi
            if abs(phase_diff[np.int0(pixelDims[0]/2),np.int0(pixelDims[1]/2)]) > np.pi:
                phase_diff -= 2*np.pi*np.sign(phase_diff[np.int0(pixelDims[0]/2),np.int0(pixelDims[1]/2)])
       
            #apply mask
            phase_masked = phase_diff * dsv350
            phase_masked = phase_masked[np.nonzero(phase_masked)]
            
            # convert to ppm and calc range within mask
            b0ppm = phase_masked/(2*np.pi*15.66*1e-3*63.88177)
            range_dsv350.append(np.percentile(b0ppm, 99) - np.percentile(b0ppm, 1))
            
    # Check if output is within tolerance
    range_dsv350 = np.array(range_dsv350)
    angle_passed = range_dsv350 < tolerance
    results.addBool('B0_Gantry_Passed',all(angle_passed))
    
    # Create figure
    figtitle = 'B0 Gantry, difference with G0 [ppm] '+dcmInfile.info.StudyDate+' '+dcmInfile.info.StudyTime
    fig, axs = plt.subplots(1,1,figsize=(7,6))
    fig.suptitle(figtitle)
    major_ticks_plotY = np.arange(0,0.51,0.1)
    minor_ticks_plotY = np.arange(0,0.51,0.05)
    major_ticks_plotX = np.arange(-180,181,30)
    axs.plot(gantry_angles, range_dsv350,linewidth=2,marker='o',markersize=5)
    axs.set_xticks(major_ticks_plotX)
    axs.set_yticks(major_ticks_plotY)
    axs.set_yticks(minor_ticks_plotY, minor=True)
    axs.grid(which='minor', alpha=0.2,axis='y')
    axs.grid(which='major', alpha=0.8,axis='both')
    
    filename = 'B0_Gantry.png'
    fig.savefig(filename,dpi=150)
    results.addObject("B0_Gantry_figure", filename)
    