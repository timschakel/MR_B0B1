#!/usr/bin/env python
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
# This code is an analysis module for WAD-QC 2.0: a server for automated 
# analysis of medical images for quality control.
#
# The WAD-QC Software can be found on 
# https://bitbucket.org/MedPhysNL/wadqc/wiki/Home
# 
#
# Changelog:
#   20220629: Initial version
#
# runfile('/nfs/arch11/researchData/USER/tschakel/projects/wadqc/QAtests/MRI_B0B1/MR_B0B1/MR_B0B1_wadwrapper.py', args='-r results.json -c config/module_config.json -d /nfs/arch11/researchData/USER/tschakel/projects/wadqc/QAtests/MRI_B0B1/MR_B0B1/testdata/data1', wdir='/nfs/arch11/researchData/USER/tschakel/projects/wadqc/QAtests/MRI_B0B1/MR_B0B1')

# this will fail unless wad_qc is already installed
from wad_qc.module import pyWADinput
import MR_B0B1_lib

# import matplotlib
# matplotlib.use('Agg') # Force matplotlib to not use any Xwindows backend.

if __name__ == "__main__":
    data, results, config = pyWADinput()
    
    # Log which series are found
    data_series = data.getAllSeries()
    print("The following series are found:")
    for item in data_series:
        print(item[0]["SeriesDescription"].value+" with "+str(len(item))+" instances")
    
    """
    Perform the analysis of the B0 and B1 measurements.
    """
    
    # read runtime parameters for module
    # TO DO: read data once and pass to all functions instead of repeating for every test
    for name,action in config['actions'].items():
        if name == 'acqdatetime':
            MR_B0B1_lib.acqdatetime(data, results, action)
        elif name == 'B1_tra':
            MR_B0B1_lib.B1_tra(data, results, action)
        elif name == 'B1_tra_AFI':
            MR_B0B1_lib.B1_tra_AFI(data, results, action)
        elif name == 'B0_tra':
            MR_B0B1_lib.B0_tra(data, results, action)
        elif name == 'B0_tra_noshim':
            MR_B0B1_lib.B0_tra_noshim(data, results, action)

    results.write()
