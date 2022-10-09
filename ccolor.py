#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul  6 15:50:17 2022

@author: canis
"""

import os
import time
import glob
import importlib
import tempfile
import numpy as np
import astropy
import astropy.io.fits as fits
from astropy.table import Table, hstack
from astropy.wcs import WCS
from astro import astroutils
from astro import mastrometry
from astro import mimage
from astro import mphot
from utils import mutils
from utils import mplot as plt
from utils import mlogging

from cconstants import holmberg_ra, holmberg_dec, comp_stars, astrometry_image_templ

log = mlogging.getLogger('ccolor')

def plot_color(starcat_file, image, VR_max= 999):
    fn = image.split('/')[-1]
    image = mimage.Image(image)
    t = Table.read(starcat_file, format='ascii.fixed_width_two_line')
    log.info(f'{len(t)} rows in table')
    i = 0
    
    xs = []
    ys = []
    
# =============================================================================
#     while i < len(t): 
#         try:
#             if float(t[i]['V']) - float(t[i]['R']) > VR_max:
#                 raise ValueError
#             xs.append(float(t[i]['V']) - float(t[i]['R']))
#             ys.append(float(t[i]['R']) - float(t[i]['I']))
#             i += 1
#         except ValueError:
#             t.remove_row(i)
#     log.info(f'{len(t)} good rows in table')
#     plt.scatter(xs, ys)
#     plt.xlabel('V-R')
#     plt.ylabel('R-I')
#     plt.title('Color plot for stars within 5\'\' of Holmberg')
#     out = f'{starcat_file.split(".")[-2]}_with_V-r_below_{VR_max}.ascii'
# =============================================================================
    out='gaia_coordinates_with_pix.ascii'
    t.sort('d')
    
    log.info('Converting Ra and dec to pixels')
    t2 = Table({'file': [], 'x': [], 'y': []}, dtype=(str, float, float))
    for row in t:
        (x, y) = mphot.equatorial2pix((row['RA'], row['DEC']), image)
        t2.add_row([fn, x, y])
    
    t = hstack([t2, t])
    t.write(out, format='ascii.fixed_width_two_line', overwrite=True)
    log.info(f'Wrote to {out}')
    return t

                                    
                                    
        

    
    