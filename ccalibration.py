#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 21 14:21:49 2022

Calibrate all the Holmberg II X-1 images in a pathectory
@author: canis
"""

import os
import time
import glob
import math
import numpy as np
import xarray as xr

from astropy.io import fits

from astro import astroutils
from astro import mphot
from astro import mimage
from utils import mutils
from utils import mplot as plt
from utils import mlogging

from icecream import ic

from astro.mconstants import *

log = mlogging.getLogger( name='ccalibration', console=True, filename='ccalibration_log.txt')

def generate_master_calibration(path):
    """
    Get master calibration data. If does not exist it generates it using calibration images.
    
    Args: 
        path: Path of pathectory of images for one night.
    Returns: 
        Dictionary containg master calibration data.
    """
    masters = mphot.master_calibration_files(path)
    if len(masters['mbias']) == 0:
        log.info('Generating master calibration images')
        mphot.masterbias(path)
        mphot.masterdark(path)
        mphot.masterflat(path)
        masters = mphot.master_calibration_files(path)
    else:
        log.info('Found existing master calibration images')
    return masters

def calibrate_night(path, format='nc',
                ftemplate='g??d???.[0-9][0-9][0-9]*', 
                astrometry=False, anetcfg=None, 
                verbose=False, overwrite=False):
    """
    Calibrate all the science images in the folder assuming calibration files are created.
    Creates master calibration files and then runs calibration on each image.
        and saves the result to `path/calibration`.
        
    Args:
        path: Path of pathectory of images for one night.
    Returns: 
        Nothing.
    """
    outpath = f'{path}/calibrated'
    
    
    if not os.path.exists(outpath):
        os.system(f'mkdir {outpath}')
    elif len(glob.glob(f'{outpath}/*')) > 0:
        if not overwrite:
            log.info('Night is already calibrated')
            return

    fnames = mutils.fixpath(path)
    files = sorted(glob.glob(f'{fnames}/{ftemplate}'))
    
    masters = generate_master_calibration(path)
    mbias_data = mphot.get_mbias_data(masters['mbias'])
    mdarks_data = mphot.get_mdarks_data(masters['mdarks'])
    mflats_data = mphot.get_mflats_data(masters['mflats'])
    
    num = len(files)
    num_skip = 0
    log.info(f'Found {num} images matching pattern {ftemplate}')
    log.info('Searching for science images')
    sciencefiles = []
    for fn in files:
        hdr, data = mimage.image2datadict(fn)
        if data is None:
            continue
        expinfo = mphot.get_exposure_info(hdr)
        name = str(expinfo['object'])
        utc  = expinfo['utc']
        filt = expinfo['filter']
        expt = expinfo['expt']
        b = expinfo['binning']
        if not "Holmberg II X-1" in name:
            if verbose:
                log.info(name + ' is not Holmberg II X-1 in' + fn)
            num_skip += 1
            continue
        if 'calibrat' in hdr['obstype'][0].lower():
            continue
        if ('bias' in hdr['obstype'][0].lower()
            or 'bias' in hdr['object'][0].lower()):
            continue
        if ('dark' in hdr['obstype'][0].lower()
            or 'dark' in hdr['object'][0].lower()
            or 'initialized' in hdr['object'][0].lower()
            or ('shutstat' in hdr and 'close' in hdr['shutstat'][0].lower())):
            continue
        if 'flat' in hdr['obstype'][0].lower():
            continue
        sciencefiles.append(fn)
    remaining = num - num_skip
    log.info(f'Remaining: {remaining} science images')

    for path in sciencefiles:
        log.info(f'Calibrating {path}')
        fn = path.split('/')[-1]
        calibrated = mphot.calibrate_image(ident=path, biasdata=mbias_data, darksdata = mdarks_data, flatsdata=mflats_data)
        fout = f'{outpath}/{fn}.c.{format}'
        #calibrated.save(fname=fout, format=format)
        if astrometry:
            acalibrated = calibrate_astrometry(calibrated, anetcfg=anetcfg)
            calibrated.hdr = acalibrated.hdr
            fout = f'{outpath}/{fn}.a.{format}'
        calibrated.data -= 32768
        calibrated.save(fname=fout, format=format)
    log.info(f'Done calibrating {path}')

def calibrate_year(path, astrometry=False):
    """
    Calibrate all the nights in a year using calibrate_night.
    
    Args: 
        yearpath: Path of pathectory of images for one year.
    Returns: 
        Nothing.
    """
    nights = sorted(glob.glob(f'{path}/???/'))
    for night in nights:
        log.info(f'Calibrating {night}')
        calibrate_night(night, astrometry=astrometry)
    log.info(f'Done calibrating {path}')
    
    
def calibrate_astrometry(input, outpath=None, anetcfg=None, verbose=False):
    #FIXME anetcfg
    """
    Recalibrate a calibrated file with astrometry
    
    Args: 
        path: Path to calibrated image
    Returns: 
        Nothing.
    """
    if isinstance(input, str):
        log.info(f'Astrometry {input.split("/")[-1]}')
    acalibrated = mphot.calibrate_image(ident=input, flags='a',
                                       anetdir='/usr/local/cellar/astrometry-net/0.91/bin',
                                       anetcfg=anetcfg, # '/opt/nofs/astrometry/etc/astrometry.cfg'
                                       sigmalimit=False, verbose=verbose)
    if np.ndim(acalibrated.data) > 2:
        acalibrated.data = acalibrated.data[0]
    if outpath is None:
        return acalibrated
    else:
        acalibrated.save(outpath)
    
def calibrate_year_astrometry(scope, year):
    """
    Recalibrate all calibrated files with astrometry
    
    Args: 
        yearpath: Path of pathectory of images for one year.
    Returns: 
        Nothing.
    """
    files = sorted(glob.glob(f'{images_path[scope]}/{year}/???/calibrated/{calibrated_image_templ}'))
    log.info(f'Found {len(files)} calibrated images')
    for fn in files:
        name = fn.split('/')[-1]
        log.info(f'{name}')

        if len(glob.glob(f'{images_path[scope]}/{year}/???/calibrated/{name[:-4]}a.nc')) > 0:
            log.info('skipping')
            continue
        calibrate_astrometry(fn)

    log.info('Done')
    
def remove_bad_astrometry(path):
    #files = sorted(glob.glob(f'{images_path[scope]}/{year}/???/calibrated/{astrometry_image_templ}', recursive=True))
    files = sorted(glob.glob(f'{path}/*a*', recursive=True))
    log.info(f'Found {len(files)} calibrated images with astrometry')
    bad_files = []
    for fn in files:
        image = mimage.Image(fn)
        hdr, data = mimage.image2datadict(image)
        if not 'a_order' in hdr:
            bad_files.append(fn)
            log.error(f'{fn} does not have astrometry. Removing')
            os.system(f'rm {fn}')
    log.info(len(bad_files))
    return bad_files
    
def fix_calibration_bug(fn):
    img_a = mimage.Image(fn)
    img_c = mimage.Image(fn[:-5]+'.c.nc')
    img_c.hdr = img_a.hdr
    img_c.save(fn)
      
# =============================================================================
# def subtract_32768_file(fn):
#     image = mimage.Image(fn)
#     if image.data[500][500] < 32678:
#         log.info(f'{fn} already got subtracted. Skipping')
#         return
#     log.info(f'Subtracting 32768 from {fn}')
#     image.data -= 32768
#     image.save(fn)
# 
# def subtract_32768(path):
#     files = sorted(glob.glob(f'{path}/**/{astrometry_image_templ}', recursive=True))
#     log.info(f'Found {len(files)} calibrated images with astrometry')
#     count = 0
#     for fn in files:
#         image = mimage.Image(fn)
#         if np.ndim(image.data) > 2:
#             if image.data[0][500][500] < 32678:
#                 log.info(f'{fn} already got subtracted. Skipping')
#                 continue
#         else:
#             if image.data[500][500] < 32678:
#                 log.info(f'{fn} already got subtracted. Skipping')
#                 continue
#         log.info(f'Subtracting 32768 from {fn}')
#         count += 1
#         image.data -= 32768
#         image.save(fn)
#     log.info(f'Subtracted 32768 from {count} out of {len(files)} files')
# 
#     
# def move_all_non_Holmbergs(scope, year):
#     """
#     Move all calibrated images that aren't Holmberg to a subfolder.
#     I updated calibrate_night so that this function shouldn't be needed.
    
#     Args: 
#         yearpath: Path of pathectory of images for one year.
#     Returns: 
#         Nothing.
#     """
#     nights = sorted(glob.glob(f'{images_path[scope]}/{year}/???/'))
#     for night in nights:
#         log.info(f'Moving non holmberg images for {night}')
#         non_holmberg = []
#         images = sorted(glob.glob(f'{night}calibrated/{calibrated_image_templ}'))
#         for fn in images:
#             hdr, data = mimage.image2datadict(fn)
#             if data is None:
#                 continue
#             expinfo = mphot.get_exposure_info(hdr)
#             name = expinfo['object']
#             utc  = expinfo['utc']
#             filt = expinfo['filter']
#             expt = expinfo['expt']
#             b = expinfo['binning']
#             if not "olmberg" in name:
#                 log.info(f'{name} is not Holmberg II X-1 in {fn}')
#                 non_holmberg.append(fn)
#         if len(non_holmberg) is not 0:
#             outpath = f'{night}calibrated/nonholmberg/'
#             if not os.path.exists(outpath):
#                 os.system(f'mkdir {outpath}')
#             for fn in non_holmberg:
#                 os.system('mv {fn} {outpath}')              
#     log.info(f'Done moving all non Holmberg images for {scope, year}')
# =============================================================================