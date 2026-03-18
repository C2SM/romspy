#! /bin/bash

indir=/net/meso/work/loher/CESM_output/GCB_2024/GCB_RunA_2024_gfort/ocn/hist
outdir=/net/sea/work/datasets/gridded/ocean/3d/model/romspy_input
yr_shift=1750
# List of variables, separated by comma:
vars=DIC

cd $outdir
for myr in `seq 1 273`; do
    yr=$(( $myr + $yr_shift ))
    printf -v myrf '%03d' $myr
    fin=$indir/GCB_RunA_2024_gfort.pop.h.0$myrf.nc
    ncra -h -v $vars $fin CESM_DIC_$yr.nc
    echo annual DIC done for $yr
done

echo concatenating annual means of DIC 
ncrcat -h CESM_DIC_*.nc CESM_DIC_ann_1751-2023.nc

echo done!