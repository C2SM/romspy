#! /usr/bin/env python

import sys
#sys.path.append('/home/loher/.local/lib/python3.9/site-packages')
sys.path.append('/home/loher/python/romspy')

from romspy import PreProcessor, forcing_adjustments

# Path to ROMS grid file:
ROMS_grid = "/net/kryo/work/loher/ROMSOC/grd/pactcs30_grd.nc"
# Path where output files are generated:
outfile = "/net/kryo/work/loher/ROMSpy_output/frc_cyclic_test/pactcs30_frc.nc"

sources = [
    # ERA5 variables:
    {
        'ROMS_setup': 'pactcs30',
        'data_source': 'ERA5',
        'variables': [
            {'out': 'sustr', 'in': 'ewss'},
            {'out': 'svstr', 'in': 'nsss'},
            {'out': 'swrad', 'in': 'ssr'},
            {'out': 'shflux', 'expr': 'shflux=sshf+slhf+str;'},
            {'out': 'evap', 'in': 'e'},
            {'out': 'precip', 'in': 'tp'},
            {'out': 'SST', 'in': 'sst'}
        ],
        'base_folder': '/net/kryo/work/updata/ecmwf-reanalysis/era5_netcdf',
        'interpolation_method': 'bil',
        'time_resolution': '1d_1h',   # '1d': daily forcing, '4h': 4h resolution
        'start_year': 1979,
        'end_year': 1979,
        'start_year_run': 1979,
        'end_year_run': 2020,
        'use_cyclic_time_axes': True, # whether the time axes in the forcing are cyclic or not
        'auxiliary_folder': '/net/kryo/work/loher/ROMSpy_output/auxiliary_files',
    },
    # COADS variables: SSS and those needed for dQdSST
    {
        'data_source': 'COADS05',
        'variables': [
            {'out': 'SSS', 'in': 'salinity', 'files': ['/net/kryo/work/updata/Roms_tools_data/COADS05/sss_landfill.cdf'] },
            {'out': 'SST', 'in': 'sst', 'files': ['/net/kryo/work/updata/Roms_tools_data/COADS05/sst_landfill.cdf'] },
            {'out': 'humidity', 'in': 'qsea', 'files': ['/net/kryo/work/updata/Roms_tools_data/COADS05/qsea_landfill.cdf'] },
            {'out': 't_air', 'in': 'sat', 'files': ['/net/kryo/work/updata/Roms_tools_data/COADS05/sat_landfill.cdf'] },
            {'out': 'rho_air', 'in': 'airdens', 'files': ['/net/kryo/work/updata/Roms_tools_data/COADS05/airdens_landfill.cdf'] },
            {'out': 'u_air', 'in': 'w3', 'files': ['/net/kryo/work/updata/Roms_tools_data/COADS05/w3_landfill.cdf'] },
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

my_preprocessor = PreProcessor(ROMS_grid, outfile, sources, hc=250.0, tcline=250.0, theta_s=10.0,
                               theta_b=4.0, fill_missing=True, verbose=True)
my_preprocessor.adjustments = forcing_adjustments
my_preprocessor.mark_as_vectors('sustr','svstr')
my_preprocessor.make()
