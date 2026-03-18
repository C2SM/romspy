#! /usr/bin/env python

# Run within conda env "romspy".

import sys
#sys.path.append('/home/loher/.local/lib/python3.9/site-packages')
sys.path.append('/home/loher/python/romspy')

from romspy import PreProcessorFrc, forcing_adjustments
from romspy import PreProcessorClm, clim_adjustments
from romspy import PreProcessorIni, ini_adjustments
from romspy import PreProcessorBry, bry_adjustments
from romspy import clim_adjustments
from romspy import UP_data_paths

# Name of ROMS setup:
ROMSsetup = 'pactcs30'
# Path to ROMS grid file:
ROMS_grid = "/net/sea/work/loher/ROMSOC/grd/pactcs30_grd.nc"
# Folder where output files will be stored:
outdir = '/net/sea/work/loher/ROMSpy_output/clm02'

# ROMS atm forcing:
# The following list of dictionaries contains the information about what data is to
# be used for which forcing variable.
sources_frc = [
    # ERA5 variables:
    {
        'ROMS_setup': ROMSsetup,
        'data_source': 'ERA5',
        'variables': [
            {'out': 'sustr', 'in': 'ewss'},
            {'out': 'svstr', 'in': 'nsss'},
            {'out': 'swrad', 'in': 'ssr'},
            {'out': 'shflux', 'in': ["sshf","slhf","str"], 'expr': 'shflux=sshf+slhf+str;'},
            {'out': 'evap', 'in': 'e'},
            {'out': 'precip', 'in': 'tp'},
            {'out': 'SST', 'in': 'sst'}
        ],
        'base_folder': UP_data_paths.UP_data_dir+'/gridded/atmosphere/2d/reanalysis/era5/',
        'interpolation_method': 'bil',
        'time_resolution': '1d_1h',   # '1d': daily forcing, '4h': 4h resolution
        'start_year': 1979,
        'end_year': 1979,
        'start_year_run': 1979,
        'end_year_run': 2020,
        'use_cyclic_time_axes': True, # whether the time axes in the forcing are cyclic or not
        'auxiliary_folder': '/net/sea/work/loher/ROMSpy_output/auxiliary_files',
    },
    # COADS variables: SSS and those needed for dQdSST
    {
        'data_source': 'COADS05',
        'variables': [
            {'out': 'SSS', 'in': 'salinity', 'files': [UP_data_paths.UP_data_dir+'/gridded/ocean/2d/observation/sss/coads05/sss_landfill.cdf'] },
            {'out': 'SST', 'in': 'sst', 'files': [UP_data_paths.UP_data_dir+'/gridded/ocean/2d/observation/sst/coads05/sst_landfill.cdf'] },
            {'out': 'humidity', 'in': 'qsea', 'files': [UP_data_paths.UP_data_dir+'/gridded/atmosphere/2d/observation/humidity/coads05/qsea_landfill.cdf'] },
            {'out': 't_air', 'in': 'sat', 'files': [UP_data_paths.UP_data_dir+'/gridded/atmosphere/2d/observation/surface_temp/coads05/sat_landfill.cdf'] },
            {'out': 'rho_air', 'in': 'airdens', 'files': [UP_data_paths.UP_data_dir+'/gridded/atmosphere/2d/observation/rho/coads05/airdens_landfill.cdf'] },
            {'out': 'u_air', 'in': 'w3', 'files': [UP_data_paths.UP_data_dir+'/gridded/atmosphere/2d/observation/wind/coads05/w3_landfill.cdf'] },
        ],
        'interpolation_method': 'bil',
    },
    # {
    #     'variables': [('dust', 'DSTSF')],
    #     'files': ['/net/kryo/work/updata/roms_tools_data/dust/dst79gnx_gx3v5.nc_inputCCSM_fill.nc'],
    #     'interpolation_method': 'bil'
    # },
    # {
    #     'variables': [
    #         {'out': 'sustr', 'in': 'tauy'},
    #         {'out': 'svstr', 'in': 'taux'}
    #     ],
    #     'files': ['/work/nicomuen/ERAinterim_tau_landfill.nc'], # This file got deleted, find a new one
    #     'interpolation_method': 'bil',
    # }
]

sources_clm = [
    {
        'variables': [
            {'out': 'temp', 'in': 'temp', 'vertical': True},
            {'out': 'salt', 'in': 'salt', 'vertical': True},
            {'out': 'u', 'in': 'u', 'vertical': True},
            {'out': 'v', 'in': 'v', 'vertical': True},
            {'out': 'zeta', 'in': 'ssh', 'vertical': False},
        ],
        'files': [UP_data_paths.UP_data_dir+'/gridded/ocean/3d/reanalysis/soda/SODA2.1.6/for_romstools/SODA_2.1.6_1979-2008_clm_landfill.nc'],
        'interpolation_method': 'bil',
    },
    # Biology 2d variables:
    {
        'var_group': 'Bio_BEC_2d',
        'variables': [
            # Chl [mg/m^3]:
            {'out': 'TOT_CHL_SURF',
             'in': 'chlorophyll',
             'files': [UP_data_paths.UP_data_dir+'/gridded/ocean/2d/observation/chl/seawifs_romstools/chla_mon_landfill.cdf']
            },
            
        ],
        'vertical': False,
        'interpolation_method': 'bil',
    },
    # Biology 3d variables:
    {
        'var_group': 'Bio_BEC_3d',
        'variables': [
            {'out': 'Alk',
             'in': 'AT_NNGv2',
             'files': [UP_data_paths.UP_data_dir+'/gridded/ocean/3d/observation/alk/broullon_2019/AT_NNGv2_climatology_time_first.nc']
            },
            {'out': 'DIC',
             'in': 'TCO2_NNGv2LDEO',
             'files': [UP_data_paths.UP_data_dir+'/gridded/ocean/3d/observation/dic/broullon_2020/TCO2_NNGv2LDEO_climatology_time_first.nc'],
             # Year represented by the climatology given by 'files':
             'year': 1995,
            },
        ],
        'vertical': True,
        'interpolation_method': 'bil',
    },
]

print("Output folder: "+outdir)
# Forcing files:
preproc_frc = PreProcessorFrc(ROMSsetup, outdir, ROMS_grid, sources_frc, layers=64, hc=250.0,
                           tcline=250.0, theta_s=10.0, theta_b=4,
                           fillmiss_after_hor=True, verbose=1)
preproc_frc.adjustments = forcing_adjustments
preproc_frc.mark_as_vectors('u','v')
preproc_frc.obc = [1, 0, 0, 0] # open at [S, E, N, W]
preproc_frc.make()

# Create clim files, with 64 vertical levels:
preproc_clm = PreProcessorClm(ROMSsetup, outdir, ROMS_grid, sources_clm, layers=64, hc=250.0,
                           tcline=250.0, theta_s=10.0, theta_b=4,
                           fillmiss_after_hor=True, verbose=1)
preproc_clm.adjustments = clim_adjustments
preproc_clm.mark_as_vectors('u','v')
preproc_clm.obc = [1, 0, 0, 0] # open at [S, E, N, W]
preproc_clm.make()

#sys.exit(0)

# Initial condition:
import glob
ini_year = 1979
clm_files = glob.glob(f"{outdir}/pactcs30_clm_*.nc")
preproc_ini = PreProcessorIni(ROMSsetup, outdir, ROMS_grid, ini_year, clm_files)
preproc_ini.adjustments = ini_adjustments
preproc_ini.make()

# Boundary files:
preproc_bry = PreProcessorBry(clm_files, ROMSsetup, ROMS_grid,
                              obc=[1,0,0,0], layers=64, hc=250.0, tcline=250.0,
                              theta_s=10.0, theta_b=4.0, verbose=True)
preproc_bry.adjustments = bry_adjustments
preproc_bry.make()
