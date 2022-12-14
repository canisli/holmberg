#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 22 12:14:37 2022

Calculate the seeing for all the calibrated Holmberg images
Save results to csv file
Plot histogram of seeings
@author: canis


Parameters:
    scope : 13 or 61
"""

import os
import time
import glob
import math
import numpy as np
import xarray as xr

from astropy.io import fits
from astropy.table import Table, vstack

from astro import astroutils
from astro import mphot
from astro import mimage
from utils import mutils
from utils import mplot as plt
from utils import mlogging

from icecream import ic

astrometry_image_templ = '*.a.*'
log = mlogging.getLogger(name='cseeing', console=True, filename='cseeing_log.txt')

def calc_seeing_image(image, quick=False, verbose=False):
    """
    Retrieves average seeing in arcseconds calculated from centroid finding algorithm mimage.Image.stars()

    Args: 
        image : mimage.Image
    Returns: 
        Tuple (Seeing in arcseconds as a float, number of stars)
    """
    if isinstance(image, str):
        image = mimage.Image(image)
    stars_data = image.stars(frame=1, snr=10, nstd=3.0, fraction=0.05, nfwhmsdev=1.8,
               diskradius=4,
               autoradius='auto fwhm', Nthresh=1.8, bufferannulus=7, skyannulus=10,
               halfwidth=1.5, emax=0.4, fwhmpixmin=1.5, fwhmpixmax=15, fit2d='moffat',
               partialpixels=True, quick=quick, detection_x=None, detection_y=None,
               detection_r1=None, detection_r2=None, verbose=verbose, maxcpus=None,
               timeout=None, tmax=300, Nmax=None )

    hdr, data = mimage.image2datadict(image)
    
    pixscale = np.nan
    if 'pixscal1' in hdr:
        pixscale = float(hdr['pixscal1'][0])
    elif 'pixscale' in hdr:
        pixscale = float(hdr['pixscale'][0])

    num_stars = len(stars_data[0])
    
    fwmean_pix = stars_data[2]
    fwmean_arc = fwmean_pix * pixscale
    
    
    log.info(f'{num_stars} stars were found')
    log.info(f'The average seeing is {fwmean_pix} pixels = {fwmean_arc} arcsec')
    
    return (fwmean_arc, num_stars)

def calc_seeing_night(path, outpath=None,
    quick=False, overwrite=False):
    """
    Retrieve average seeing in arcsecond for each image in a night. Option to save results to seeings.csv

    Args: path : Path of directory of calibrated images for one night
           overwrite : Boolean for whether to ignore existing csv file and recalculate seeings
    Returns: 
        Astropy Table of image names and seeings
    
    """
    ext = "_quick" if quick else ""
    bad_files = []
    if not overwrite:
        seeing_file = glob.glob(path + f'/seeings{ext}.ascii')
        if len(seeing_file) > 0:
            log.info(f'Night already has seeings{ext}.ascii. Fetching ascii file')
            return (Table.read(seeing_file[0], format='ascii.fixed_width_two_line'), 
                    bad_files)
    files = sorted(glob.glob(f'{path}/*.a.*'))
    seeings = Table({'Image': [], 'Seeing (arcsec)' : [], 'Num Stars': []}, 
                    dtype=('str', 'float', 'int32'))
    if len(files) == 0:
        log.info('Folder is empty')
        return (seeings, bad_files)
    for fn in files:
        fn_name = fn.split('/')[-1]
        log.info(fn_name)
        image = mimage.Image(fn)
        seeing = -1
        num_stars = -1
        #try:
        (seeing, num_stars) = calc_seeing_image(image, quick=quick)
        # except Exception as e:
        #     log.error(e)
        #     bad_files.append(fn)
        seeings.add_row([fn_name, seeing, num_stars])
    if outpath is not None:
        seeings.write(f'{outpath}/seeings{ext}.ascii', format='ascii.fixed_width_two_line', overwrite=True)
        log.info(f'Wrote results to {path}/seeings{ext}.ascii')
    return (seeings, bad_files)

def get_avg_fwhm_year(scope, year, quick=False, overwrite=False):
    """
    Run get_avg_fwhm_night on all the night dirs in a year
    """
    
    log.info(f'Using quick? {quick}')
    log.info(f'Overwrite? {overwrite}')
    nights = sorted(glob.glob(f'{images_dir[scope]}/{year}/???/calibrated'))
    total_bad_files = []
    for night in nights:
        log.info(f'Getting avg fwhms for {night}')
        #try:
            
        seeings, bad_files = get_avg_fwhm_night(night, quick=quick, save=True, overwrite=overwrite)
        for x in bad_files:
            total_bad_files.append(x)
        #except Exception as e:
        #    log.error(e)
    log.info(f'Done getting fwhms for {images_dir[scope]}/{year}')
    log.info(total_bad_files)

def combine_seeings(scope, year, quick=False, save=True):
    """
    Take all the seeing.csv files generated by get_avg_fwhm_year and combine them into one master file.
    
    Args: 
        yeardir: Path of directory of images for one year.
        save : Boolean for whether to write results to csv file 
    Returns: 
        Nothing.
    """
    ext = "_quick" if quick else ""
    seeings_files = sorted(glob.glob(f'{images_dir[scope]}/{year}/???/calibrated/seeings{ext}.ascii'))
    log.info(f'Found {len(seeings_files)} seeings{ext}.ascii files')
    log.info('Combining files')
    datas = []
    
    for fn in seeings_files:
        night_number = fn.split('/')[-3]
        data = Table.read(fn, format='ascii')
        if len(data) > 0:
            log.info(f'Appending {night_number}')
            datas.append(data)
        else:
            log.info(f'Skipped night because empty')
    
    total_data =vstack(datas)
        
    if save:
        total_data.write(f'{images_dir[scope]}/{year}/total_seeings{ext}.ascii', format='ascii.fixed_width_two_line', overwrite=True)
        log.info(f'Wrote results to {images_dir[scope]}/{year}/total_seeings{ext}.ascii')
        
        # sort by filename and seeing

    return total_data

def construct_histogram(scope, year, quick=False):
    """
    Read all the seeings for a year and draw a histogram of the seeings rounded to the nearest decimal point
    
    Args: 
        total_data: Either table returned by combine_csvs or path to total_seeings.csv file
    Returns: 
        Nothing.
    """
    ext = '_quick' if quick else ''
    
    total_data = Table.read(f'{images_dir[scope]}/{year}/total_seeings{ext}.ascii', 
                            format='ascii.fixed_width_two_line')
    seeings = np.round(np.array(total_data['Seeing (arcsec)']), decimals=1)
    plt.figure()
    plt.histoplot(seeings, colors=['blue'], xlab='Seeing (arcsec)', 
                  ylab='Frequency (# images)', 
                  xrange=(-1, max(seeings)),
                  dogrid=True, plottitle=f'Seeings{ext} {scope} {year} for Holmberg II X-1 Obserations in arcseconds')

def filter_seeings_camera(scope,cameras=['nd12','tek2k'], max_seeings=[2.0,2.5,3.0,3.5,4.0,6.0]):
    for max_seeing in max_seeings:
        log.info(max_seeing)
        for camera in cameras:
            outpath = f'{data_dir}/cameras/{camera}/seeing_below_{max_seeing}'
            log.info(f'Output folder: {outpath}') 
            if not os.path.exists(outpath):
                os.makedirs(outpath)
        files = sorted(glob.glob(f'{images_dir[scope]}/????/seeing_below_{max_seeing}/{astrometry_image_templ}'))
        
        for fn in files:
            hdr, data = mimage.image2datadict(mimage.Image(fn))
            camera = str(hdr['instrume'][0])
            os.system(f'cp {fn} {data_dir}/cameras/{camera}/seeing_below_{max_seeing}')

def total_seeings(max_seeing):
    files = sorted(np.concatenate(
        [glob.glob(f'/home/canis/data/13/images/acam/seeing_below_{max_seeing}/{astrometry_image_templ}'),
         glob.glob(f'/home/canis/data/13/images/Marana/seeing_below_{max_seeing}/{astrometry_image_templ}'),
         glob.glob(f'/home/canis/data/61/images/tek2k/seeing_below_{max_seeing}/{astrometry_image_templ}'),
         glob.glob(f'/home/canis/data/61/images/nd12/seeing_below_{max_seeing}/{astrometry_image_templ}'),]
        ))
    
    outpath = f'/home/canis/data/all_seeing_below_{max_seeing}'
    log.info(f'Output folder: {outpath}') 
    if not os.path.exists(outpath):
        os.makedirs(outpath)
    for fn in files:
        os.system(f'cp {fn} {outpath}')
        
        
        

def filter_seeings(scope, years=[2019,2020,2021], max_seeings=[2.0, 2.5, 3.0, 3.5, 4.0, 6.0], overwrite=False, verbose=False):
    """
    Check the seeing of each calibrated image and if it is below `max_seeing` copy it to the output folder
    
    Args: 
        scope: 13 or 61
        year
        max_seeing: Maximum seeing allowed to filer images.
    Returns: 
        Nothing.
    """
    if not isinstance(years, list):
        years = [years]
    for year in years:
        for max_seeing in max_seeings:
            outpath = f'{images_dir[scope]}/{year}/seeing_below_{max_seeing}'
            log.info(f'Output folder: {outpath}') 
            if not os.path.exists(outpath):
                os.makedirs(outpath)
            else:
                if not overwrite:
                    log.info('Already exists.')
                    continue
                
            total_data = Table.read(f'{images_dir[scope]}/{year}/total_seeings.ascii', 
                                    format='ascii.fixed_width_two_line')
            
            files = sorted(glob.glob(f'{images_dir[scope]}/{year}/???/calibrated/{astrometry_image_templ}'))
            # Both files and total_data are sorted by alphabetical order of the name of the image
            log.info(f'Found {len(files)} calibrated files with astrometry')
            log.info(f'Filtering the ones with seeing below {max_seeing}')
            
            num_files = len(total_data)
            log.info(f'Total rows in table {num_files}')
            
            good_files = []
            index = -1
            bad_seeing_count = 0
            # wrong_camera_count = 0
            # good_camera = 'tek2k'
            for fn in files: # slow
                #image = mimage.Image(fn)
                #hdr, data = mimage.image2datadict(image)
                name = fn.split('/')[-1]
        # =============================================================================
        #         camera = str(hdr['instrume'][0])
        #         if camera != good_camera:
        #             wrong_camera_count += 1
        #             log.debug(f'{name} is not using {good_camera}. Skipping')
        #             continue
        # =============================================================================
                while total_data[index]['Image'][:-5] != name[:-5]:
                    index += 1
                row = total_data[index]
                seeing = float(row['Seeing (arcsec)'])
                if 0 < seeing < max_seeing:
                    if verbose:
                        log.debug(f'{name} has a seeing {seeing} between 0 and {max_seeing}. Adding and copying')
                    os.system(f'cp {fn} {outpath}/')
                    good_files.append(fn)
                else:
                    bad_seeing_count += 1
                    if verbose:
                        log.debug(f'{name} has a seeing {seeing} not between 0 and {max_seeing}. Skipping')
                    
            log.info('SUMMARY')
            log.info(f'Found {len(good_files)} out of {num_files} images with seeing below {max_seeing}')
            log.info(f'Skipped {bad_seeing_count} files with seeing above {max_seeing}')
            #log.info(f'Skipped {wrong_camera_count} files taken without {good_camera}')
            log.info(f'Copied them over to {outpath}')
            
            log.info('REMEMBER TO COMBINE SEEINGS!')

def count_files(fdir, exp=astrometry_image_templ):
    return len(glob.glob(f'{fdir}/{exp}'))