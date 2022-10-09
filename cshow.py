#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 26 13:48:09 2022

@author: canis
"""
import glob
import pickle
import numpy as np

from utils import mplot as plt

def show(seeing=3.0, plots=[1,2,3,4,5,6,7,8,9,10], interval=None, version=5, ext='_optimize', mean=False):
    print(f'_v{version}{ext}')
    if mean:
        fdir=f'/home/canis/data/photresults_mean_v{version}{ext}'
    else:
        fdir=f'/home/canis/data/photresults_v{version}{ext}'
    if isinstance(plots,int):
        plots = [plots]
    if plots == 'astropy':
        plots = [1, 5,6,7]
    if interval is None or (interval[0] == -np.inf and interval[1] == np.inf):
        print(interval)
        files = sorted(glob.glob(f'{fdir}/{seeing}/*.plt'))
    else:
        print(interval)
        files = sorted(glob.glob(f'{fdir}/{seeing}/{interval[0]}-{interval[1]}/*.plt'))
    i = 0
    for fn in files:
        i += 1
        if i in plots:
            fig = pickle.load(open(f'{fn}', 'rb'))
            fig.show()

def showc(seeing=3.0, plots=[1,2,3,4,5,6,7,8,9,10], c=1, interval=None, version=5, ext='_optimize'):
    print(f'_v{version}{ext}')
    fdir=f'/home/canis/data/photresults_v{version}{ext}'
    if c=='all':
        c = [1,2,3,4,5,6,7,8,9,10]
    if isinstance(plots,int):
        plots = [plots]
    if plots == 'astropy':
        plots = [1, 5,6,7]
    if interval is None or (interval[0] == -np.inf and interval[1] == np.inf):
        print(interval)
        files = sorted(glob.glob(f'{fdir}/{seeing}/comp{c}/*.plt'))
    else:
        print(interval)
        files = sorted(glob.glob(f'{fdir}/{seeing}/{interval[0]}-{interval[1]}/comp{c}/*.plt'))
    i = 0
    for fn in files:
        i += 1
        if i in plots:
            fig = pickle.load(open(f'{fn}', 'rb'))            
            fig.show()
            fig.savefig(f'/home/canis/fig{i}.png')

def close():
    plt.close('all')