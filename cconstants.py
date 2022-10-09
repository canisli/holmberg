#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 24 16:01:39 2022
Constants

@author: canis
"""

images_dir = {13: '/home/canis/data/13/images', 61:'/home/canis/data/61/images'}
data_dir = '/home/canis/data'
calibrated_image_templ = '*.c.nc'
astrometry_image_templ = '*.a.nc'

holmberg_ra = '08 19 28.99'
holmberg_dec = '70 42 19.4'
      
# =============================================================================
# comp_stars = [('08 19 22.11' , '70 42 22.765'),  #1 bright star SIMBAD
#               ('08 19 58.466', '70 43 18.408'),  #2 bright star SIMBAD
#               ('08 19 35.292', '70 42 46.714'),  #3 slightly bright star SIMBAD   
#               ('08 19 31.496', '70 44 21.869'),  #4 blue star   GAIA
#               ('08 19 06.009', '70 41 48.596'),  #5 blue star   GAIA
#               ('08 18 50.035', '70 42 51.406'),  #6 blue star   GAIA
#               ('08 19 42.614', '70 40 23.9842')  # V - R is not < 0.2 GAIA
#              ]
# 
# =============================================================================

holmberg= ('08 19 28.99', '70 42 19.4')
comp_stars = [#('08 19 22.08335', '70 42 22.44100'),  #1
             #('08 19 58.48721', '70 43 18.39930'),  #2 saturated, commented
              ('08 19 35.27861', '70 42 46.68284'),  #3
              ('08 19 31.496', '70 44 21.869'),      #4 dim blue star   
              ('08 19 06.009', '70 41 48.596'),      #5 dim blue star   
              ('08 18 50.035', '70 42 51.406'),      #6 dim blue star   
              ('08 19 42.614', '70 40 23.9842'),     #7 NEW
             # ('08 19 13.55795', '70 42 27.63916'),  #8 NEW similar to #2
              ('08 19 57.87424', '70 42 53.07892'),  #9 NEW nearby #1
              ('08 19 55.25860', '70 42 08.06211'),   #10
            ]


"""
from astropy.coordinates import SkyCoord
import astropy.units as u

holmberg = SkyCoord(ra=ra, dec=dec, unit=(u.h, u.deg))

wcs.world_to_pixel(holmberg)
"""