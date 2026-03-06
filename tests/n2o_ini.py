#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb  7 09:27:23 2023

@author: jahaerri
"""

import xarray as xr
import netCDF4 as nc4
import numpy as np
import cartopy.crs as ccrs
import matplotlib.pyplot as plt

#import site
import sys
#sys.path.insert(0, site.USER_SITE)
sys.path.append('/home/loher/python/ROMSpy/romspy-0.9.0')

from romspy import PreProcessor, clim_adjustments

target = '/net/kryo/work/loher/ROMS/ROMSpy_n2o/humpac15_grd.nc'
out = '/net/kryo/work/loher/ROMS/ROMSpy_n2o/humpac_n2o_ini.nc'
sources = [
    {
        'variables': [
            {'out': 'N2O', 'in': 'N2O', 'vertical': True}
        ],
        'files': ['/net/kryo/work/jahaerri/roms/inputs/humpac15_Ncycle/n2o_gridded_landfill.nc'],
        'interpolation_method': 'bil',
    }
]

roms_clim = PreProcessor(target, out, sources, theta_s = 10.0, theta_b = 4.0, layers = 64, hc = 250.0,
                         tcline = 250.0, sigma_type = 3, use_ROMS_grdfile=True, verbose=True)
roms_clim.adjustments = clim_adjustments
roms_clim.obc = [0, 0, 0, 1] # open at [S, E, N, W]

roms_clim.make()

print('at end of n2o_ini.py')

"""
#combine n2o and no2 initial conditions with inital conditions from Eike
n2o = xr.open_dataset("/net/kryo/work/jahaerri/roms/inputs/humpac15_Ncycle/humpac_n2o_ini_0_0.nc").variables["N2O"].values
#o2 = xr.open_dataset("/net/kryo/work/jahaerri/roms/inputs/humpac15_Ncycle/ini/spinup_r103_humpac15_rst_bec2phys_timereset.nc").variables["O2"].values
fn_mask = '/nfs/kryo/work/koehne/roms/output/humpac15/hindcast_1979_2019/hindcast_r105_humpac15/daily/avg/humpac15_1979_avg.nc'
mask_rho = xr.open_dataset(fn_mask).variables['mask_rho'][:].values
lon = xr.open_dataset(fn_mask).variables['lon_rho'][:].values
lat = xr.open_dataset(fn_mask).variables['lat_rho'][:].values
xr.open_dataset(fn_mask).attrs

n2o_nan = np.where((mask_rho ==1) & (np.isnan(n2o[0,0,:,:])==True), 1, np.nan)

fig = plt.figure(figsize = (10,10), dpi = 300)
ax = plt.axes(projection = ccrs.Robinson(-110))
ax.set_global()
ax.pcolormesh(lon,lat, n2o_nan, vmin = 0, vmax = 6, transform = ccrs.PlateCarree())

ds = nc4.Dataset("/net/kryo/work/jahaerri/roms/inputs/humpac15_Ncycle/ini/spinup_r103_humpac15_rst_bec2phys_timereset.nc", 'a', format='NETCDF4')

n2o_nc = ds.createVariable("N2O", 'f8', ('time', 's_rho', 'eta_rho', 'xi_rho'))
n2o_nc[:,:,:] = n2o[:,:,:,:]
n2o_nc.units = "mMol/m3"
n2o_nc.long_name = "Nitrous Oxide"
      
ds.close()

#add boundary conditions
n2o_bry = np.array([n2o[0,:,:,0]]*12)

ds_bry = nc4.Dataset("/net/kryo/work/jahaerri/roms/inputs/humpac15_Ncycle/bry/humpac15_bry_replaced.nc", 'a', format='NETCDF4')

n2obry_nc = ds_bry.createVariable("N2O_west", 'f8', ('bry_time', 's_rho', 'eta_rho'))
n2obry_nc[:,:,:] = n2o_bry[:,:,:]
n2obry_nc.units = "mMol/m3"
n2obry_nc.long_name = "Nitrous Oxide"
  
ds_bry.close()
"""
