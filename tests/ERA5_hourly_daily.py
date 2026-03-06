#! /usr/bin/env python

import netCDF4
import cartopy.crs as ccrs
import matplotlib.colors as colors
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy as np
import os
#import pdb

# Computes daily averages from hourly ERA5 data for some user specified
# variables and compares them to daily ERA5 data. Currently this
# is done only for Jan 1, 1979.

# Folder where plots will be stored:
OUTDIR = '/net/kryo/work/loher/plots/ERA5_test'
# ERA5 hourly and daily data files:
ERADIR = '/net/kryo/work/updata/ecmwf-reanalysis/era5_netcdf'
DATA_1h = '{}/hourly/1979/ERA5_1979_01.nc'.format(ERADIR)
DATA_1d = '{}/daily/1979/ERA5_1979_01_daily.nc'.format(ERADIR)
# ERA5 variables to be checked:
ERA_VARS = ['sst','ewss','str']
# Flag indicating whether a given variable is an accumulated
# quantity in ERA5:
ERA_IS_ACC = {'ewss': True, 'str': True, 'sst': False}
# Units of ERA5 time averages to be used in the plots:
ERA_AVG_UNITS = {'ewss': 'N m**-2', 'str': 'W m**-2', 'sst': 'K'}
# ERA5 unit conversion factors:
#ERA_UNIT_CONV = {'ewss': 0.001, 'str': 1}
ERA_UNIT_CONV = {'ewss': 1, 'str': 1, 'sst': 1}
# Long names of ERA5 variables:
ERA_LONGNAME = {'ewss': 'Eastward turbulent surface stress',
                'str':  'Surface net thermal radiation',
                'sst':  'Sea surface temperature' }
# Flag indicating whether a white centered colorbar is to be
# used for a given variable:
WHITE_CENTERED_CBAR = {'ewss': True, 'str': True, 'sst': False}
# Matplotlib colormaps to use for "normal" and white centered maps:
CMAP_NORMAL = 'nipy_spectral'
CMAP_WCENTERED = 'seismic'

# Function for producing map plots
def map_plot_global(data1,data2,lat,lon,title1,title2,white_centered,fig_title,unit,out_file):
    fig, ax = plt.subplots(2,2, subplot_kw=dict(projection=ccrs.Mercator()),
                                sharex=True, sharey=True, figsize=(12,9))
    fig.tight_layout()
    #fig.suptitle(fig_title, fontsize=12)
    # Compute max and min of the 2 data arrays:
    max1 = data1.max()
    max2 = data2.max()
    vmax = max(max1,max2)
    min1 = data1.min()
    min2 = data2.min()
    vmin = min(min1,min2)
    if white_centered:
        norm = colors.TwoSlopeNorm(vmin=vmin, vcenter=0, vmax=vmax)
        cmap = CMAP_WCENTERED
    else:
        norm = colors.Normalize(vmin=vmin, vmax=vmax)
        cmap = CMAP_NORMAL
    # Create plots:
    proj = ccrs.PlateCarree()
    # For data1:
    cs1 = ax[0,0].pcolormesh(lon,lat,data1,norm=norm,transform=proj,cmap=cmap)
    ax[0,0].set_title(title1, fontsize=9)
    ax[0,0].coastlines()
    divider = make_axes_locatable(ax[0,0])
    ax_cb = divider.new_horizontal(size="5%", pad=0.1, axes_class=plt.Axes)
    fig.add_axes(ax_cb)
    plt.colorbar(cs1, cax=ax_cb)
    ax_cb.tick_params(labelsize=9)
    #ax_cb.set_ylabel(unit, fontsize=8)
    ax_cb.minorticks_on()
    # data2:
    cs2 = ax[0,1].pcolormesh(lon,lat,data2,norm=norm,transform=proj,cmap=cmap)
    ax[0,1].set_title(title2, fontsize=9)
    ax[0,1].coastlines()
    divider = make_axes_locatable(ax[0,1])
    ax_cb = divider.new_horizontal(size="5%", pad=0.1, axes_class=plt.Axes)
    fig.add_axes(ax_cb)
    plt.colorbar(cs2, cax=ax_cb)
    ax_cb.tick_params(labelsize=9)
    ax_cb.set_ylabel(unit, fontsize=9)
    ax_cb.minorticks_on()
    # Difference:
    diff = data1-data2
    min_diff = diff.min()
    max_diff = diff.max()
    if min_diff < 0 and max_diff > 0:
        norm = colors.TwoSlopeNorm(vmin=min_diff, vcenter=0, vmax=max_diff)
        cmap = CMAP_WCENTERED
    else:
        norm = colors.Normalize(vmin=vmin, vmax=vmax)
        cmap = CMAP_NORMAL
    cs3 = ax[1,0].pcolormesh(lon,lat,diff,norm=norm,transform=proj,cmap=cmap)
    ax[1,0].set_title('{} - {}'.format(title1,title2), fontsize=9)
    ax[1,0].coastlines()
    divider = make_axes_locatable(ax[1,0])
    ax_cb = divider.new_horizontal(size="5%", pad=0.05, axes_class=plt.Axes)
    fig.add_axes(ax_cb)
    plt.colorbar(cs3, cax=ax_cb)
    ax_cb.tick_params(labelsize=9)
    ax_cb.set_ylabel(unit, fontsize=9)
    ax_cb.minorticks_on()
    # Plot title is shown in the lower right axes:
    ax[1,1].set_axis_off()
    ax[1,1].text(0.5, 0.5, fig_title, horizontalalignment='center',
        verticalalignment='center', fontweight='bold',
        fontsize=12, transform=ax[1,1].transAxes)
    # Save plot:
    print('save file: {}'.format(out_file))
    plt.savefig(out_file, bbox_inches='tight')
    plt.close()


if not os.path.exists(OUTDIR):
    os.makedirs(OUTDIR)

# The following computation of stimei (to get time averages from accumulated
# quantities) is taken from the Romstools:
# Number of time steps per day in ERA5 hourly files:
stepday_1h = 24
# Summing time of ERA data of data [h]:
stime_1h = 24/stepday_1h
# Conversion factor to convert to time averages, for hourly accumulated data:
stimei_1h = 1/(stime_1h*3600)  # [sec^-1]
# For daily ERA5 data:
stepday_1d = 1
stime_1d = 24/stepday_1d
stimei_1d = 1/(stime_1d*3600)  # [sec^-1]

nc_h = netCDF4.Dataset(DATA_1h, 'r')
nc_d = netCDF4.Dataset(DATA_1d, 'r')
lat = nc_h.variables['latitude'][:]
lon = nc_h.variables['longitude'][:]
for var in ERA_VARS:
    print('var = '+var)
    if ERA_IS_ACC[var]:
        # Here we need stimei to get time averages:
        fac_h = stimei_1h
        fac_d = stimei_1d
    else:
        fac_h = 1
        fac_d = 1
    # Average 24 time records for Jan 1, to get the daily average:
    dat_h = fac_h*ERA_UNIT_CONV[var]*np.mean(nc_h.variables[var][:24,:], axis=0)
    # Read daily data for Jan 1 and compute daily average:
    dat_d = fac_d*ERA_UNIT_CONV[var]*nc_d.variables[var][0,:]
    # Plot results:
    fname = '{}/{}_avg_1979_01_01.png'.format(OUTDIR,var)
    map_plot_global(dat_h,dat_d,lat,lon,var+' hourly',var+' daily',WHITE_CENTERED_CBAR[var],
                    ERA_LONGNAME[var]+'\nDaily average for Jan 1, 1979',ERA_AVG_UNITS[var],fname)
nc_h.close()
nc_d.close()
