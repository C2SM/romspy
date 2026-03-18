#! /usr/bin/env python

import sys
#sys.path.append('/home/loher/.local/lib/python3.9/site-packages')
sys.path.append('/home/loher/python/romspy')

import glob
from romspy import PreProcessor, clim_adjustments
from romspy import PreProcessorIni, ini_adjustments
from romspy import PreProcessorBry, bry_adjustments

# UP data path:
updata = '/net/kryo/work/updata'
# Path to ROMS grid file:
target_grid = "/net/kryo/work/loher/ROMSOC/grd/pactcs30_grd.nc"
# Folder where output files will be stored:
outdir = '/net/kryo/work/loher/ROMSpy_output/clm02'

# Climatology:
outfile_clm = "{}/pactcs30_clm.nc".format(outdir)
sources = [
    {
        'variables': [
            {'out': 'temp', 'in': 'temp', 'vertical': True},
            {'out': 'salt', 'in': 'salt', 'vertical': True},
            {'out': 'u', 'in': 'u', 'vertical': True},
            {'out': 'v', 'in': 'v', 'vertical': True},
            {'out': 'zeta', 'in': 'ssh', 'vertical': False},
        ],
        'files': [updata+'/soda_clim/soda_2.1.6/SODA_2.1.6_1979-2008_clm_landfill.nc'],
        'interpolation_method': 'bil',
    },
]

# Clim file, with 64 vertical levels:
preproc_clm = PreProcessor(target_grid, outfile_clm, sources, layers=64, hc=250.0,
                           tcline=250.0, theta_s=10.0, theta_b=4,
                           fillmiss_after_hor=True, verbose=True)
preproc_clm.adjustments = clim_adjustments
preproc_clm.mark_as_vectors('u','v')
preproc_clm.obc = [1, 0, 0, 0] # open at [S, E, N, W]
preproc_clm.make()

# Initial condition:
outfile_ini = "{}/pactcs30_ini.nc".format(outdir)
outfiles_clm = glob.glob("{}/pactcs30_clm_*.nc".format(outdir))
preproc_ini = PreProcessorIni(outfile_ini, outfiles_clm, layers=64, hc=250.0,
                              tcline=250.0, theta_s=10.0, theta_b=4.0,
                              verbose=True)
preproc_ini.adjustments = ini_adjustments
preproc_ini.make()

# Boundary files:
outfile_bry = "{}/pactcs30_bry.nc".format(outdir)
preproc_bry = PreProcessorBry(outfile_bry, outfiles_clm, 'Pactcs30', target_grid,
                              obc=[1,0,0,0], layers=64, hc=250.0, tcline=250.0,
                              theta_s=10.0, theta_b=4.0, verbose=True)
preproc_bry.adjustments = bry_adjustments
preproc_bry.make()
