import cdo
import netCDF4
import os
import numpy as np
from itertools import count
import sys

from romspy.interpolation.interpolator import ShiftPairCollection
from romspy.verification import test_cdo, has_vertical, verify_sources
from .grid_routines import scrip_grid_from_nc
from .interpolation.vertical import sigma_stretch_sc, sigma_stretch_cs, get_z_levels
from .interpolation import Interpolator

"""
Author: Nicolas Munnich
License: GNU GPL2+
"""

# Class for creating clim files:
class PreProcessorClm:
    def __init__(self, ROMSsetup: str, outdir: str, target_grid: str, sources: list, **kwargs):
        """
        Contains universal methods

        :param ROMSsetup: name of ROMS setup, e.g. 'Pactcs30'
        :param outdir: folder where output files are to be stored
        :param target_grid: ROMS grid file
        :param sources: Information about which variables should be taken from which files using
                which interpolation method onto which type of grid_routines and whether to interpolate them vertically
        :param scrip_grid: Scrip file on rho grid_routines to interpolate with. Will be created if not provided
        :param verbose: whether text should be printed as the program is running
        Can have any of the following optional arguments:
            theta_s - S-coordinate surface control parameter - default 7.0
            theta_b - S-coordinate bottom control parameter - default 0.0
            layers - Number of S-coordinate layers - default 32
            hc - S-coordinate critical depth - default 150
            tcline - S-coordinate surface/bottom layer width - default 150
            sigma_type - default 3
            zeta_source - dict with keys 'name' 'file'. Sets zeta to the first timestep of 'name' in 'file'.
            file_type - output filetype -  default nc4c
            processes - number of processes cdo should use - default 8
            scrip_grid - SCRIP version of target_grid - optional, will be created if not passed
            verbose - level of verbosity - default 1
            time_underscored - whether 'time' variables should be replaced with '#_time' - default False
            keep_weights - whether to keep calculated weights - default False
            keep_z_clim - whether to keep zclim - default False
            use_ROMS_grdfile - whether the ROMS grid file is to be used instead of a SCRIP format grid file - default False
        """
        print()
        print('\033[1;32m===================\033[0m')
        print('\033[1;32mClimatology files\033[0m')
        print('\033[1;32m===================\033[0m')
        print()
        # Level of verbosity:
        verbose = kwargs.get('verbose', 1)
        # cdo
        cdo_debug = (verbose>2)
        self.cdo = cdo.Cdo(debug=cdo_debug)
        test_cdo(self.cdo)
        # Sources
        #verify_sources(sources, kwargs.get('verbose', False))
        self.sources, self.target_grid = sources, target_grid

        # Vertical interpolation information
        # Replace anything not passed in with default values
        self.has_vertical = has_vertical(self.sources)
        self.theta_s, self.theta_b, self.layers, self.hc, self.tcline, self.sigma_type = (
            kwargs.get('theta_s', 7.0), kwargs.get('theta_b', 0.0),
            kwargs.get("layers", 32), kwargs.get("hc", 150),
            kwargs.get('tcline', 150), kwargs.get("sigma_type", 3)
        )
        self.sc = sigma_stretch_sc(self.layers, True)
        self.cs = sigma_stretch_cs(self.theta_s, self.theta_b, self.sc, self.sigma_type)

        # Get z_levels
        with netCDF4.Dataset(target_grid, mode='r') as my_grid:
            self.h = my_grid.variables['h'][:]
        self.zeta = self.get_zeta(kwargs['zeta_source']) if 'zeta_source' in kwargs else np.zeros_like(self.h)
        self.z_level_rho, self.z_level_u, self.z_level_v = get_z_levels(self.h, self.sc, self.cs, self.hc,
                                                                        self.zeta)

        # CDO options
        self.file_type, self.processes = kwargs.get('file_type', 'nc4c'), kwargs.get('processes', 8)
        # Other Options
        self.verbose, self.time_underscored, self.keep_weights, self.keep_z_clim = (
            kwargs.get('verbose', 1), kwargs.get('time_underscored', False),
            kwargs.get('keep_weights', False), kwargs.get('keep_z_clim', False)
        )
        self._adjustments = None
        self.ROMS_setup = ROMSsetup
        self.outfile = f'{outdir}/{ROMSsetup}_clm.nc'
        self.use_ROMS_grdfile = kwargs.get('use_ROMS_grdfile', False)
        # Fill in missing values after horizontal interplolation:
        self.fillmiss_after_hor = kwargs.get('fillmiss_after_hor', False)
        # Boolean or list of variables for which to extrapolate to land at the end:
        self.fill_missing = kwargs.get('fill_missing', True)
        # Interpolator
        self.scrip_grid = kwargs.get('scrip_grid', scrip_grid_from_nc(target_grid,outdir))
        self.shift_pairs = ShiftPairCollection()

        if self.verbose > 0:
            print("Finished PreProcessorClm setup")

    @property
    def options(self):
        return " -b F32 -f " + self.file_type + " -P " + str(self.processes)

    @property
    def adjustments(self):
        return self._adjustments

    @adjustments.setter
    def adjustments(self, value: list):
        if not isinstance(value, list):
            print("ERROR: Your adjustments are in an incorrect format. Adjustments must be a list of dictionaries!")
            return
        for adj in value:
            if not isinstance(adj, dict):
                print("ERROR: Your adjustments are in an incorrect format. Adjustments must be a list of dictionaries!")
                print(adj)
                return
            if not isinstance(adj.get('in_var_names'), set):
                print("ERROR: An adjustment has an incorrect 'in_var_name' key. "
                      "All dicts must have the key 'in_var_names' pointing to a set of input variable strings")
                print(adj)
                return
            if not isinstance(adj.get('out_var_names'), set):
                print("ERROR: An adjustment has an incorrect 'out_var_name' key. "
                      "All dicts must have the key 'out_var_names' pointing to a set of output variable strings")
                print(adj)
                return
            if adj.get('func') is None:
                print("ERROR: All dicts must have the key 'func' pointing to a function!")
                print(adj)
                return
        self._adjustments = value

    def make(self):
        if self.adjustments is None:
            print("Please set your adjustments first!")
            return

        print("PreProcessorClm --> make called")
        # dict of all variables produced
        all_vars = {var['out'] for sublist in [x['variables'] for x in self.sources] for var in sublist}
        # For each group of variables
        for group, group_index in zip(self.sources, count()):
            if "var_group" in group:
                print("var group = "+group["var_group"])
            elif "data_source" in group:
                print("data source = "+group["data_source"])
            #else:
            #    msg = "group found without neither 'var_group' nor 'data_source' att"
            #    raise ValueError(msg)
            if not 'data_source' in group:
                group['data_source'] = None
            variables = group['variables']
            # set of all variables present in out_file after interpolation
            out_variables = {i["out"] for i in variables}
            fidx_iter = count(start=1)
            if 'files' in group:
                lsources = [self.sources[group_index]]
                # Set up interpolator:
                interp = Interpolator(self.cdo, os.path.split(self.outfile)[0], lsources, self.target_grid,
                            self.scrip_grid,
                            (self.z_level_rho, self.z_level_u, self.z_level_v) if self.has_vertical else None,
                            self.shift_pairs,
                            self.options, '', self.keep_weights, self.keep_z_clim, self.use_ROMS_grdfile,
                            self.fillmiss_after_hor, self.verbose)
                # For each file associated with the group of variables
                for in_file, file_index in zip(group['files'], fidx_iter):
                    # Get the unique filename for each file
                    out_file = '{}_{}_{:03}.nc'.format(self.outfile[:-3],group_index,file_index)
                    # Continue if output file exists already:
                    #print("infile = "+in_file)
                    if os.path.exists(out_file):
                        continue
                    # Interpolate everything necessary
                    interp.interpolate(in_file, out_file, group, group_index, variables, group['files'])
                    # Make any adjustments to variables necessary
                    self.make_adjustments(out_file, out_variables, group_index, group['files'], all_vars)
                    # Rename time to starting with an underscore if necessary
                    if self.time_underscored:
                        self.__rename_time(out_file, group_index)
                    # add 1D information
                    if has_vertical(self.sources):
                        self.add_1d_attrs(out_file)
            else:
                # Input files are given per variable:
                outfiles = []
                lsources = [self.sources[group_index]]
                for x in variables:
                    for in_file, file_index in zip(x['files'], count(start=1)):
                        # Set up interpolator:
                        interp = Interpolator(self.cdo, os.path.split(self.outfile)[0], lsources, self.target_grid,
                            self.scrip_grid,
                            (self.z_level_rho, self.z_level_u, self.z_level_v) if self.has_vertical else None,
                            self.shift_pairs,
                            self.options, in_file, self.keep_weights, self.keep_z_clim, self.use_ROMS_grdfile,
                            self.verbose)
                        out_file = '{}_{}_{}_{:03}.nc'.format(self.outfile[:-3],x['out'],
                                                              group_index,file_index)
                        outfiles.append(out_file)
                        # Continue if output file exists already:
                        if os.path.exists(out_file):
                            continue
                        #out_file = self.outfile[:-3] + '_' + x['out'] + '_' + \
                        #    str(group_index) + '_' + '{:03}'.format(str(file_index)) + '.nc'
                        var_files = ','.join(x['files'])
                        lvars = [x]

                        # Interpolate everything necessary
                        print(f"before interp.interpolate: group = {group}")
                        out_file = interp.interpolate(in_file, out_file, group, group_index, lvars, var_files)
                        # Make any adjustments to variables necessary
                        #self.make_adjustments(out_file, out_variables, group_index, var_files, all_vars)
                        # Rename time to starting with an underscore if necessary
                        if self.time_underscored:
                            self.__rename_time(out_file, group_index)
                        # add 1D information
                        if has_vertical(self.sources):
                            self.add_1d_attrs(out_file)
                # Merge per variable cdo output files:
                file_index = next(fidx_iter)
                out_file = self.outfile[:-3] + '_' + str(group_index) + '_' + str(file_index) + '.nc'
                if not os.path.exists(out_file):
                    print("create file: "+out_file)
                    print(f"    from files: {outfiles}")
                    self.cdo.merge(input=(' '.join(outfiles)), output=out_file, options=self.options)
                # Make any adjustments to variables necessary
                print(f"adjustements: out_file = {out_file}")
                print(f"adjustements: out_variables = {out_variables}")
                self.make_adjustments(out_file, out_variables, group_index, group['files'], all_vars)

    @staticmethod
    def __rename_time(file, group_index):
        with netCDF4.Dataset(file, mode='r+') as my_file:
            my_file.renameDimension('time', str(group_index) + '_time')
            my_file.renameVariable('time', str(group_index) + '_time')

    def add_1d_attrs(self, file_name):
        vertical = [
            {'name': 'theta_s', 'long_name': 'S-coordinate surface control parameter', 'datatype': 'f',
             'dimensions': 'one', 'units': '-', 'data': self.theta_s},
            {'name': 'theta_b', 'long_name': 'S-coordinate bottom control parameter', 'datatype': 'f',
             'dimensions': 'one', 'units': '-', 'data': self.theta_b},
            {'name': 'Tcline', 'long_name': 'S-coordinate surface/bottom layer width', 'datatype': 'f',
             'dimensions': 'one', 'units': 'meter', 'data': self.tcline},
            {'name': 'hc', 'long_name': 'S-coordinate critical depth', 'datatype': 'f',
             'dimensions': 'one', 'units': 'meter', 'data': self.hc},
            {'name': 'sc_r', 'long_name': 'S-coordinate at RHO-points', 'datatype': 'f',
             'dimensions': 's_rho', 'units': '-', 'data': self.sc},
            {'name': 'Cs_r', 'long_name': 'S-coordinate stretching curve at RHO-points', 'datatype': 'f',
             'dimensions': 's_rho', 'units': '-', 'data': self.cs}
        ]
        with netCDF4.Dataset(file_name, mode="r+") as my_file:  # with automatically opens and closes
            for var in vertical:
                #my_file.setncattr(var['name'], str(var['long_name'] + " := " + str(var['data'])))
                my_file.setncattr(var['name'], var['data'])

    def add_time_underscores(self):
        self.time_underscored = True

    def make_adjustments(self, file: str, out_variables: set, group_index: int, group_files: str, all_vars: set):
        if self.verbose > 0:
            print("Making adjustments to file contents as per adjustments.")

        for adjustment in self.adjustments:
            # If a variable in out_var_names is not preprovided or it doesn't produce variables
            if len(adjustment['out_var_names'] - all_vars) > 0 or len(adjustment['out_var_names']) == 0:
                # If the file has all the necessary inputs
                if len(adjustment['in_var_names'] & out_variables) == len(adjustment['in_var_names']):
                    # If the output isn't already pre-calculated and all the inputs are in the same file
                    if self.verbose > 0:
                        print("Calling " + str(adjustment['func'].__name__))
                        sys.stdout.flush()
                        # import pdb
                        # pdb.set_trace()
                        adjustment['func'](file, group_files=group_files, group_index=group_index, options=self.options,
                                            adjustments=self.adjustments, **vars(self))
                        # Update out_variables, if this adjustment has added a variable:
                        out_variables |= adjustment['out_var_names']
                elif len(adjustment['in_var_names'] & out_variables) >= 1:
                    print("Warning: The variables needed to calculate " + str(adjustment['out_var_names']) +
                          " were not all present in the same file! Variables needed: " +
                          str(adjustment['in_var_names']))

    def get_zeta(self, zeta_source):
        temp = self.cdo.remapbil(self.scrip_grid, input="-selname," + zeta_source['name'] + " " + zeta_source['file'],
                                 options=self.options)
        with netCDF4.Dataset(temp) as zetafile:
            zeta = np.array(zetafile.variables[zeta_source['name']][0])
        return zeta

    def mark_as_vectors(self, vector_u_name, vector_v_name):
        self.shift_pairs.add_shift_pair(vector_u_name, vector_v_name)
