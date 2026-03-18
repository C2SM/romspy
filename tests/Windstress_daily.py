#! /usr/bin/env python

# Create 1h wind stress forcing files for Humpac15 setup.

import sys
#sys.path.append('/home/loher/.local/lib/python3.9/site-packages')
sys.path.append('/home/loher/python/romspy')

from romspy import PreProcessorFrc, forcing_adjustments
from romspy import PreProcessorClm, clim_adjustments
from romspy import PreProcessorIni, ini_adjustments
from romspy import PreProcessorBry, bry_adjustments
from romspy import clim_adjustments

# Name of ROMS setup:
ROMSsetup = 'humpac15'
# Path to ROMS grid file:
ROMS_grid = "/nfs/sea/work/jahaerri/roms/inputs/humpac15_Ncycle/grd/humpac15_grd.nc"
# Folder where output files will be stored:
outdir = '/net/sea/work/loher/ROMSpy_output/Windstress_daily'

# The following list of dictionaries contains the information about what data is to
# be used for which variable.
sources = [
    # 1) Setting related to atm forcing:
    # ERA5 variables:
    {
        'ROMS_setup': ROMSsetup,
        'data_source': 'ERA5',
        'variables': [
            {'out': 'sustr', 'in': 'ewss'},
            {'out': 'svstr', 'in': 'nsss'},
        ],
        'base_folder': '/net/sea/work/datasets/gridded/atmosphere/2d/reanalysis/era5/',
        'interpolation_method': 'bil',
        'time_resolution': '1d',   # '1d': daily forcing, '4h': 4h resolution
        'start_year': 1979,
        'end_year': 1980,
        'start_year_run': 1979,
        'end_year_run': 2020,
        'use_cyclic_time_axes': True, # whether the time axes in the forcing are cyclic or not
        'auxiliary_folder': '/net/sea/work/loher/ROMSpy_output/auxiliary_files',
    },
]

print("Output folder: "+outdir)
# Forcing files:
preproc_frc = PreProcessorFrc(ROMSsetup, outdir, ROMS_grid, sources, layers=64, hc=250.0,
                           tcline=250.0, theta_s=10.0, theta_b=4,
                           fillmiss_after_hor=True, verbose=1)
preproc_frc.adjustments = forcing_adjustments
preproc_frc.mark_as_vectors('sustr','svstr')
preproc_frc.obc = [1, 0, 0, 0] # open at [S, E, N, W]
preproc_frc.make()
