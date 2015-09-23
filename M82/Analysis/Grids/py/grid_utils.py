# Utilities for MonteCarlo grids

from __future__ import print_function, absolute_import, division, unicode_literals

import numpy as np
import gzip
import subprocess, glob 
import json, io

from astropy import units as u
from astropy.io import ascii, fits
from astropy.table import QTable, Table, Column

from xastropy.xutils import fits as xxf

from xastropy.xutils import xdebug as xdb

def monte_to_fits(spec_fil, fits_fil):
    '''Convert MonteCarlo ASCII output to fits
    '''
    # Read (fast)
    f = open(spec_fil, 'r')
    lines = f.readlines()
    f.close()

    # Parse header
    nlam, nx, ny = [int(ii) for ii in lines[0].split(' ')]
    lam0, dlam, dx = [float(ii) for ii in lines[1].split(' ')]

    # Wavelength
    wave  = lam0 + np.arange(nlam)*dlam

    # x_arr
    x_arr = np.arange(nx)*dx

    # Generate farr for quick generation
    farr = np.array( [float(line[0:16]) for line in lines[2:]] )  
    #carr = np.array( [float(line[20:].strip()) for line in lines[2:]] )   # "Clicks"
    #xdb.set_trace()

    # Data
    data = np.zeros((nx,ny,nlam))
    flux = np.zeros(nlam)
    cnt=0
    for ii in range(nx):
        for jj in range(ny):
            data[ii,jj,:] = farr[cnt:cnt+nlam]
            cnt += nlam

    # Normalize (not really following this..)
    L_tot = np.sum(data)
    xdb.set_trace()
    nrm = 1.
    # dlambda
    if nlam > 1:
        L_tot *= (wave[1] - wave[0])
        nrm *= (wave[1] - wave[0])*nlam
                
    # spatial
    if nx > 1:
        L_tot *= (x_arr[1] - x_arr[0])**2
        nrm *= (x_arr[1] - x_arr[0])**2
    #
    print('L_tot = {:g}, nrm = {:g}'.format(L_tot,nrm))

    # Write
    xxf.write_quick_fits([data*nrm, wave, x_arr], fits_fil)

#### ########################## #########################
#### ########################## #########################
#### ########################## #########################


# Command line execution
if __name__ == '__main__':
    flg_run = 0
    flg_run += 2**0 # First try

    # First try with M82
    if (flg % 2**1) >= 2**0:
        # Param file
        param = init_param()
        param['flag_density'] = 0
        param['flag_wind'] = 0
        param['title'] = 'M82_first'
        # Generate arrays and output files
        main(param, 'M82_first')


