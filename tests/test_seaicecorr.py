#! /usr/bin/env python

# Create 1h wind stress forcing files for Humpac15 setup.

import sys
#sys.path.append('/home/loher/.local/lib/python3.9/site-packages')
sys.path.append('/home/loher/python/romspy/c++/build')

import roms_forcing_utils
import time

grd_file = "/nfs/sea/work/jahaerri/roms/inputs/humpac15_Ncycle/grd/humpac15_grd.nc"
out_dir = "/net/sea/work/loher/ROMSpy_output/ERA5_hourly_2"
roms_setup = "humpac15"
roms_frc = roms_forcing_utils.ROMS_frc(roms_setup, grd_file, out_dir, True)
flist = [f"{out_dir}/humpac15_frc_0_001.nc", f"{out_dir}/humpac15_frc_0_002.nc"]
year = 1979
seaice_file = "/net/sea/work/koehne/roms/inputs/humpac15/hindcast_1979_2019/frc/corr_fields_seaice/365days/humpac15_seaice_frc.nc"
snowice_file = "/net/sea/work/koehne/roms/inputs/humpac15/hindcast_1979_2019/frc/corr_fields_seaice/365days/humpac15_snowice_frc.nc"
time_res = 1
t1 = time.time()
roms_frc.make_seaice_correction(flist, year, seaice_file, snowice_file, time_res, True)
t2 = time.time()
print(f"time for make_seaice:correction: {t2-t1:.2f} s")
