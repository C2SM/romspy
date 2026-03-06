"""
Author: Damian Loher
License: GNU GPL2+
"""

import netCDF4
from romspy.adjustments import biology_params

# Corrects DIC based on CESM output given on a 1x1 lat-lon grid:
def dic_adjustment(file: str, **kwargs):
    # Years for which to compute DIC difference
    # More precisely, DIC[year2] - DIC[year1] is computed on the ROMS grid:
    import pdb
    pdb.set_trace()
    sources = kwargs['sources']
    group_index = kwargs['group_index']
    if not 'var_group' in sources[group_index]:
        msg = "did not find 'var_group' in sources[group_index]"
        raise ValueError(msg)
    vgrp = sources[group_index]['var_group']
    if vgrp != "Bio_BEC_3d":
        msg = f"unkown 'var_group' found in sources[group_index]: {vgrp}"
        raise ValueError(msg)
    # Look for DIC entry and find year of DIC clim used:
    clm_year = None
    for vdata in sources[group_index]['variables']:
        if vdata["out"] == "DIC":
            clm_year = vdata["year"]
    if clm_year is None:
        msg = "no information found about DIC in 'sources'"
        raise ValueError(msg)
    year1 = clm_year
    year2 = kwargs['years_second']
    if year1 < biology_params.CESM_DIC_ANN_FIRST_YR:
        msg = f'year_first must be >= {biology_params.CESM_DIC_ANN_FIRST_YR}'
        raise ValueError(msg)
    if year2 > biology_params.CESM_DIC_ANN_LAST_YR:
        msg = f'year_second must be <= {biology_params.CESM_DIC_ANN_LAST_YR}'
        raise ValueError(msg)
    # Get difference of DIC annual means:
    nc = netCDF4.Dataset(biology_params.CESM_DIC_ANN_MEANS,'r')
    t1 = year1 - biology_params.CESM_DIC_ANN_FIRST_YR
    t2 = year2 - biology_params.CESM_DIC_ANN_FIRST_YR
    delta_dic = (nc.variables['DIC'][t2,:] - nc.variables['DIC'][t1,:])
    nc.close()
    return delta_dic


bry_adjustments = [
    {'out_var_names': set(), 'in_var_names': {'DIC'},
     'func': dic_adjustment},
]