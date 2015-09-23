# Following grid_biconical.pro from P+11 Wind paper

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


def init_arrays(ngrid):
    ''' Set up arrays for grid calculations
    '''
    xval = np.arange(ngrid) - ngrid/2 + 0.5  #; Dimensionless
    cell_center = np.zeros( (ngrid,ngrid,ngrid,3) )

    # Get started
    xx = np.outer(xval, np.ones(ngrid)) 
    yy = np.outer(np.ones(ngrid), xval) # These may be flipped relative to IDL
    for k in range(ngrid):
        cell_center[:,:,k,0] = xx
        cell_center[k,:,:,1] = xx # Seems this is asymmetric..
        cell_center[k,:,:,2] = yy

    #;;;;;;;
    #;; Corners
    x_corners = np.zeros((ngrid,ngrid,ngrid,8))
    y_corners = np.zeros((ngrid,ngrid,ngrid,8))
    z_corners = np.zeros((ngrid,ngrid,ngrid,8))

    #;; x
    xc = 0.5*np.array([-1,1,-1,1,-1,1,-1,1])
    tmp = (cell_center[:,:,:,0])[:] 
    for qq in range(8):
        x_corners[:,:,:,qq] = tmp + xc[qq]
    #;; y
    yc = 0.5*np.array([1,1,-1,-1,1,1,-1,-1])
    tmp = (cell_center[:,:,:,1])[:] 
    for qq in range(8):
        y_corners[:,:,:,qq] = tmp + yc[qq]
    #;; z
    zc = 0.5*np.array([-1,-1,-1,-1,1,1,1,1])
    tmp = (cell_center[:,:,:,2])[:] 
    for qq in range(8):
        z_corners[:,:,:,qq] = tmp + zc[qq]

    # Return
    return cell_center, (x_corners, y_corners, z_corners)

def set_radii(param, cell_center, corners):
    '''Simple calculation of radii
    '''
    r_center = param['dl'] * np.reshape(np.sqrt( cell_center[:,:,:,0]**2 + 
        cell_center[:,:,:,1]**2 + cell_center[:,:,:,2]**2), 
        (param['ngrid'],param['ngrid'],param['ngrid']))
    # Corners
    r_corners = np.zeros((param['ngrid'],param['ngrid'],param['ngrid'],8))
    for xyz_corners in corners:
        r_corners += xyz_corners**2 
    r_corners = param['dl'] * np.sqrt(r_corners)
    # Return
    return r_center, r_corners

def set_density(param, grid_stuff=None, flag=1):
    '''Generate the density grid
    Parameters:
    ----------
    param: dict
      Dict summarizing grid Parameters
    flag: int, optional
      Specifies the type of density profile
      0: Constant density
      1: Power-law

    Returns:
    ----------
    rho_grid: ndarray
      Density grid
    '''
    if grid_stuff is None:
        # Init cells
        cell_center, corners = init_arrays(param['ngrid'])

        # Radius
        r_center, r_corners = set_radii(param, cell_center, corners)
        del  corners
    else: # Unpack
        cell_center, r_center, r_corners = grid_stuff

    # Good portion of the wind (radial)
    msk_corners = np.zeros(r_corners.shape,dtype='int')
    tst1 = r_corners >= param['rw_inner']
    tst2 = r_corners <= param['rw_outer']
    gd_corn = tst1 & tst2
    msk_corners[gd_corn] = 1

    # Density Law
    if flag==0:
        density = param['n0'] 
    elif flag==1:
        density = param['n0'] * (np.maximum(r_center,1e-5)/param['rw_inner'])**param['density_pow']

    # Deal with partial grid cells
    frac_cell = np.sum(msk_corners,3)
    rho_grid = density * frac_cell / 8.

    # Biconical
    theta = np.reshape( np.arctan( np.abs(cell_center[:,:,:,2])/ 
        np.sqrt(cell_center[:,:,:,0]**2 + cell_center[:,:,:,1]**2) ), rho_grid.shape) * 180. / np.pi
    #xdb.set_trace()
    zero = theta < param['theta_min']
    rho_grid[zero] = 0.

    # Return
    return rho_grid

def set_velocity(param, flg_wind=1, grid_stuff=None):
    '''Generate the density grid
    Parameters:
    ----------
    param: dict
      Dict summarizing grid Parameters
    flg_wind: int, optional
      Specifies the type of density profile
      0: Constant velocity
      1: r^-1 law

    Returns:
    ----------
    vx_grid,vy_grid,vz_grid: tuple
      Tuple of the velocity grids
    '''
    # Init cells
    if grid_stuff is None:
        if cell_center is None:
            cell_center, corners = init_arrays(param['ngrid'])

        # Radius
        if r_center is None:
            r_center, r_corners = set_radii(param, cell_center, corners)
    else: # Unpack
        cell_center, r_center = grid_stuff

    # Wind profile (always radial speed)
    if flg_wind==0:      
        speed_cell = param['v_wind']
    elif flg_wind==1: 
        speed_cell = param['v_wind'] * r_center/param['rw_inner'] 
    elif flg_wind==2: 
        speed_cell = param['v_wind'] + dv_wind*(r_cell-rw_inner)/(rw_outer-rw_inner)
    else:
        raise ValueError('Not ready for this wind profile!')

    # Grid
    ngrid=param['ngrid']
    dl = param['dl']
    vx_grid = speed_cell * (dl*cell_center[:,:,:,0]) / r_center 
    vy_grid = speed_cell * (dl*cell_center[:,:,:,1]) / r_center 
    vz_grid = speed_cell * (dl*cell_center[:,:,:,2]) / r_center 

    return (vx_grid, vy_grid, vz_grid)

def set_emiss(param, r_corners=None):
    '''Generate the emissivity grid
    Parameters:
    ----------
    param: dict
      Dict summarizing grid Parameters

    Returns:
    ----------
    emiss: ndarray
      Density grid
    '''
    # Radius
    if r_corners is None:
        cell_center, corners = init_arrays(param['ngrid'])
        r_center, r_corners = set_radii(param, cell_center, corners)
    # 
    msk_corners = np.zeros(r_corners.shape,dtype='int')
    gd_corn = r_corners < param['rg_outer']
    msk_corners[gd_corn] = 1

    frac_cell = np.sum(msk_corners,3)
    emiss = frac_cell / 8.
    # Return
    return emiss

def write_monte(param, out_arrays, outfil):
    '''Turn grid files into an ASCII file for the MonteCarlo code
    '''
    # Unpack
    rho_grid, vx_grid, vy_grid, vz_grid, emiss, dopp = out_arrays
    # Open file
    f = open(outfil,'w')
    f.write('{:d} {:g}\n'.format(rho_grid.shape[0],param['dl']))

    # Loop to loop
    nx, ny, nz = rho_grid.shape
    for ii in range(nx):
        for jj in range(ny):
            for kk in range(nz):
                f.write('{:g} {:.1f} {:.1f} {:.1f} {:.2f} {:.1f}\n'.format(
                    rho_grid[ii,jj,kk],
                    vx_grid[ii,jj,kk],
                    vy_grid[ii,jj,kk],
                    vz_grid[ii,jj,kk],
                    emiss[ii,jj,kk],
                    dopp[ii,jj,kk]))
    f.close()
    print('Wrote {:s}'.format(outfil))


def init_param(ngrid=300):

    # Some defaults
    param = dict(ngrid=ngrid)
    param['rw_inner'] = 1.            #; Inner boundary of the wind (kpc)
    param['rw_outer'] = 3.            #; Outer boundary of the wind (kpc)
    param['box_size'] = param['rw_outer'] * 2
    param['dl'] = param['box_size'] / param['ngrid'] # Cell size (kpc)
    param['n0'] =  0.0005        # Density at inner radius (Hydrogen; cm^-3)
    param['theta_min'] =  45.    # Opening angle of the cone (deg)
    param['v_wind'] = 450.       # Parameterization of the wind speed (radial; km/s)
    param['rg_outer'] = 0.1      # Outer boundary of galaxy (kpc)
    param['flag_density'] = 1    # Flag for density profile
    param['flag_wind'] = 0       # Flag for wind speed profile
    param['doppler'] = 15.       # Doppler parameter (km/s)
    #param['density_pow'] = -2.

    return param


def main(param, out_root=None):
    '''Put all the pieces together
    '''
    # Init cells
    cell_center, corners = init_arrays(param['ngrid'])
    print('Done making cells')
    r_center, r_corners = set_radii(param, cell_center, corners)
    print('Done with radius')

    # Density
    rho_grid = set_density(param, flag=param['flag_density'], grid_stuff=(cell_center,r_center,r_corners))
    print('Done with density')

    # Velocity
    vxyz_grid = set_velocity(param, flg_wind=param['flag_wind'], grid_stuff=(cell_center, r_center))
    print('Done with velocity')

    # Doppler (and more)
    gd_cell = rho_grid > 0. 
    for vgrid in vxyz_grid:
        vgrid[~gd_cell] = 0.
    dopp = np.zeros(rho_grid.shape)
    dopp[gd_cell] = param['doppler']

    # Emissivity
    emiss = set_emiss(param,r_corners=r_corners)

    # Pack it all up 
    out_arrays = [rho_grid] + [vgrid for vgrid in vxyz_grid] + [emiss] + [dopp]

    # Write files
    if out_root is not None:
        # JSON for parameter fil
        with io.open(out_root+'_param.json', 'w', encoding='utf-8') as f:
            f.write(unicode(json.dumps(param, sort_keys=True, indent=4, 
                separators=(',', ': '))))
        # FITS file 
        xxf.write_quick_fits(out_arrays, out_root+'.fits')

        # MonteCarlo type file
        write_monte(param,out_arrays,out_root+'.asc')

    return out_arrays


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


