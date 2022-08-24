#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 27 15:05:16 2022

@author: tschakel
"""
from skimage import feature
import cv2
import numpy as np
import matplotlib.pyplot as plt

def detect_edges(image, sigma=0.3, low_threshold=750, high_threshold=None
):
    """
    Detect edges on a 2d array
    :param high_threshold: high threshold for the hysteresis thresholding
    :param low_threshold: low threshold for the hysteresis thresholding
    :param sigma: width of the Gaussian
    :param image: 2d numpy array
    :return binary array of same dimensions as numpy_array representing the detected edges
    """

    # canny requires floats
    edges = feature.canny(
        image.astype("float32"),
        sigma=sigma,
        low_threshold=low_threshold,
        high_threshold=high_threshold,
    )

    return edges

def find_center(image_data,params):
    """
    alternative for retrieve_ellipse_parameters
    """
    low_thresh = np.mean(np.nonzero(image_data))*0.1
    edges = detect_edges(image_data, int(params['canny_sigma']), low_thresh)
    contours, hierarchy = cv2.findContours(edges.astype('uint8'), cv2.RETR_TREE, cv2.CHAIN_APPROX_SIMPLE)
    c = max(contours, key = cv2.contourArea) # select the biggest contour
    (center_x_pix,center_y_pix),radius_pix = cv2.minEnclosingCircle(c)
    
    return center_x_pix, center_y_pix

def create_circular_mask(image, center=None, radius=None):
    (w,h) = image.shape
    
    if center is None: # use the middle of the image
        center = (int(w/2), int(h/2))
    if radius is None: # use the smallest distance between the center and image walls
        radius = min(center[0], center[1], w-center[0], h-center[1])

    Y, X = np.ogrid[:h, :w]
    dist_from_center = np.sqrt((X - center[0])**2 + (Y-center[1])**2)

    mask = dist_from_center <= radius
    return mask

def create_spherical_mask(shape, voxel_spacing=(1.0,1.0,1.0), center=None, radius=None):
    # Create a spherical mask for height h, width w and depth d centered at center with given radius

    x,y,z = shape

    if center is None:  # use the middle of the image
        center = np.array(np.array(shape)/2,dtype=int)
    if radius is None:  # use the smallest distance between the center and image walls
        radius = min(center, shape-center)

    X,Y,Z=(np.ogrid[:x, :y, :z])
    dist_from_center = np.sqrt(((X - center[0]) * voxel_spacing[0]) ** 2  +\
                       ((Y - center[1]) * voxel_spacing[1]) ** 2  +\
                       ((Z - center[2]) * voxel_spacing[2]) ** 2)

    mask = dist_from_center <= radius
    return mask

def calc_rms_stats(pixeldataInB0ppm,mask):
    # calc stats per slice
    rms_stats = np.empty(pixeldataInB0ppm.shape[2]+1)
    for slice in range(pixeldataInB0ppm.shape[2]):
        slicedata = pixeldataInB0ppm[:, :, slice] * mask[:, :, slice]
        rms_stats[slice] = np.sqrt(np.mean(slicedata**2))
        
    rms_stats[slice+1] = np.sqrt(np.mean(pixeldataInB0ppm[mask]**2))
    
    return rms_stats
    
def calc_percentile_stats(pixeldataInB0ppm,mask):
    # calc percentile stats
    percentiles = [5,25,50,75,95]
    perc_stats = np.empty([pixeldataInB0ppm.shape[2]+1,5])
    for slice in range(pixeldataInB0ppm.shape[2]):
        slicedata = pixeldataInB0ppm[:, :, slice] * mask[:, :, slice]
        if np.max(slicedata) == 0:
            perc_stats[slice,:] = [0,0,0,0,0]
        else:
            perc_stats[slice,:] = np.percentile(slicedata[slicedata != 0],percentiles)
        
    perc_stats[slice+1,:] = np.percentile(pixeldataInB0ppm[mask],percentiles)
    
    return perc_stats

def b0_create_figure1(b0map,phantom,acqdate,acqtime,scanori,basename):
    # figure 1
    # display slices of masked fieldmap (in Hz), use p5 and p95 as clip range

    fig1data = b0map * phantom #masked image
    fig1min = np.percentile(fig1data[fig1data != 0],5)
    fig1max = np.percentile(fig1data[fig1data != 0],95)
    filename = basename+scanori+'_figure1.jpg'
    
    nslices = b0map.shape[2]
    plotcols = 5
    if nslices <=5:
        plotrange = np.arange(plotcols)
    else:
        plotstep = np.ceil(nslices/5)
        plotrange = np.int0(np.arange(0,nslices,plotstep))
        
    fig, axs = plt.subplots(ncols=plotcols, nrows=1, figsize=(10,3))
    title = basename + scanori +" "+ acqdate +" "+ acqtime
    fig.suptitle(title,fontsize=20)
    
    n = 0
    for slice in plotrange:
        ax = axs[n]
        im = ax.imshow(fig1data[:,:,slice],cmap='hot',vmin=fig1min,vmax=fig1max)
        ax.set_title('Slice '+str(slice+1))
        ax.axis('off')
        n = n+1
    
    fig.colorbar(im, shrink=0.6)
    fig.savefig(filename,dpi=150)
    return filename

def b0_create_figure2(b0map_ppm,b0map,phantom,dsv100,dsv200,dsv300,
                   dsv350,acqdate,acqtime,imaging_frequency,scanori,basename,slicenumber):
    # figure 2
    # Display the analysis on the center slice
    # mask of phantom, fieldmap in Hz, PPM map
    # overlay the r=5,10,15,17.5 masks
    fig2data = b0map * phantom #masked image
    fig2min = np.percentile(fig2data[fig2data != 0],5)
    fig2max = np.percentile(fig2data[fig2data != 0],95)

    slice=slicenumber # center slice
    filename = basename+scanori+'_figure2.jpg'
    #fig, axs = plt.subplots(ncols=4, nrows=2, sharey=True, sharex=True,
    #                        constrained_layout=True, figsize=(10,5))
    fig, axs = plt.subplots(ncols=4, nrows=2,figsize=(10,5))
    title = basename+scanori +" "+ acqdate +" "+ acqtime + " results slice "+str(slicenumber+1)
    fig.suptitle(title,fontsize=20)
    
    ax = axs[0,0]
    im = ax.imshow(phantom[:,:,slice],cmap='gray')
    ax.set_title('Mask slice '+str(slice+1))
    ax.axis('off')
    
    ax = axs[0,1]
    im = ax.imshow(b0map[:,:,slice] * phantom[:, :, slice],cmap='hot',vmin=fig2min,vmax=fig2max)
    ax.set_title('Fieldmap [Hz] slice  '+str(slice+1))
    fig.colorbar(im, ax=ax)
    ax.axis('off')
    
    ax = axs[0,2]
    im = ax.imshow(b0map_ppm[:,:,slice] * phantom[:, :, slice],cmap='hot',vmin=fig2min/imaging_frequency,vmax=fig2max/imaging_frequency)
    ax.set_title('Fieldmap [ppm] slice '+str(slice+1))
    fig.colorbar(im, ax=ax)
    ax.axis('off')
    
    ax = axs[0,3]
    ax.remove()
    
    ax = axs[1,0]
    im = ax.imshow(b0map_ppm[:,:,slice] * dsv100[:,:,slice],
                   cmap='hot',vmin=fig2min/imaging_frequency,vmax=fig2max/imaging_frequency)
    ax.set_title('r=5cm')
    ax.axis('off')
    
    ax = axs[1,1]
    im = ax.imshow(b0map_ppm[:,:,slice] * dsv200[:,:,slice],
                   cmap='hot',vmin=fig2min/imaging_frequency,vmax=fig2max/imaging_frequency)
    ax.set_title('r=10cm')
    ax.axis('off')
    
    ax = axs[1,2]
    im = ax.imshow(b0map_ppm[:,:,slice] * dsv300[:,:,slice],
                   cmap='hot',vmin=fig2min/imaging_frequency,vmax=fig2max/imaging_frequency)
    ax.set_title('r=15cm')
    ax.axis('off')
    
    ax = axs[1,3]
    im = ax.imshow(b0map_ppm[:,:,slice] * dsv350[:,:,slice],
                   cmap='hot',vmin=fig2min/imaging_frequency,vmax=fig2max/imaging_frequency)
    ax.set_title('r=17.5cm')
    ax.axis('off')
    plt.tight_layout()
    
    fig.savefig(filename,dpi=150)
    return filename

def b0_create_figure3(b0map_ppm,phantom,acqdate,acqtime,scanori,basename):
    # figure 3
    # histogram of ppm of all slices (50 bins)
    # cumulative ppm histogram
    # every slice: mean ppm with min and max as errorbars
    fig3data = b0map_ppm * phantom #masked image
    nslices = b0map_ppm.shape[2]
    
    filename = basename+scanori+'_figure3.jpg'
    fig, axs = plt.subplots(ncols=3, nrows=1,
                            constrained_layout=True, figsize=(12,3))
    title = basename+scanori +" "+ acqdate +" "+ acqtime
    fig.suptitle(title,fontsize=24)
    
    ax1 = axs[0]
    ax1.hist(fig3data[fig3data != 0].flatten(), bins=50)
    ax1.set(title='Histogram',xlabel='ppm',ylabel='number of pixels')
    
    ax2 = axs[1]
    ax1.get_shared_x_axes().join(ax1, ax2)
    ax2.hist(fig3data[fig3data != 0].flatten(), bins=50, cumulative=True,histtype='step')
    ax2.set(title='Cumulative Histogram',xlabel='ppm',ylabel='number of pixels')
    
    ax3 = axs[2]
    x = np.arange(1,nslices+1)
    slicemedian=np.empty(nslices)
    sliceuperr=np.empty(nslices)
    slicelowerr=np.empty(nslices)
    
    for slice in range(nslices):
        slicedata = fig3data[:,:,slice]
        slicemedian[slice] = np.median(slicedata[slicedata != 0])
        slicelowerr[slice] = np.percentile(slicedata[slicedata != 0],5)
        sliceuperr[slice] = np.percentile(slicedata[slicedata != 0],95)
        
    asymmetric_error = [slicelowerr-slicemedian,slicemedian-sliceuperr]
    ax3.errorbar(x,slicemedian,yerr=asymmetric_error,fmt='o-')
    ax3.set(xlabel='slice number', ylabel='ppm',title='Median ppm +/- p5/95')
    
    fig.savefig(filename,dpi=160)    
    return filename

def b0_collect_results(results,scanori,seriesname,
                       fname_fig1,fname_fig2,fname_fig3,
                       rms_dsv100,rms_dsv200,rms_dsv300,rms_dsv350,
                       perc_stats,rms_phantom,pixelDims,basename):
    results.addString(basename + "Orientation", scanori)
    results.addString(basename + "SeriesDescription", seriesname)
    results.addObject(basename + "Figure1",fname_fig1)
    results.addObject(basename + "Figure2",fname_fig2)
    results.addObject(basename + "Figure3",fname_fig3)
    
    reportkeyvals = []
    for slice in range(pixelDims[2]):
        # results for the different slices
        idname = "_slice"+str(slice+1)
        
        reportkeyvals.append( (basename + "RMS d10"+idname,rms_dsv100[slice]) )
        reportkeyvals.append( (basename + "RMS d20"+idname,rms_dsv200[slice]) )
        reportkeyvals.append( (basename + "RMS d30"+idname,rms_dsv300[slice]) )
        reportkeyvals.append( (basename + "RMS d35"+idname,rms_dsv350[slice]) )
        
        reportkeyvals.append( (basename + "p5 "+idname,perc_stats[slice,0]) )
        reportkeyvals.append( (basename + "p25 "+idname,perc_stats[slice,1]) )
        reportkeyvals.append( (basename + "p50 "+idname,perc_stats[slice,2]) )
        reportkeyvals.append( (basename + "p75 "+idname,perc_stats[slice,3]) )
        reportkeyvals.append( (basename + "p95 "+idname,perc_stats[slice,4]) )
    
    #whole phantom results
    reportkeyvals.append( (basename + "RMS phantom",rms_phantom[pixelDims[2]]) )
    reportkeyvals.append( (basename + "p5 phantom",perc_stats[pixelDims[2],0]) )
    reportkeyvals.append( (basename + "p25 phantom",perc_stats[pixelDims[2],1]) )
    reportkeyvals.append( (basename + "p50 phantom",perc_stats[pixelDims[2],2]) )
    reportkeyvals.append( (basename + "p75 phantom",perc_stats[pixelDims[2],3]) )
    reportkeyvals.append( (basename + "p95 phantom",perc_stats[pixelDims[2],4]) )
    
    for key,val in reportkeyvals:
        results.addFloat(key, val)
        
    return results