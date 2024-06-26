#!/usr/bin/env python

""" @file filter_mag_and_size_info.py

    Created 8 Oct 2015

    Filter out only magnitude and size info of galaxies from the CFHTLenS
    catalogue.

    ---------------------------------------------------------------------

    Copyright (C) 2015 Bryan R. Gillis

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""

import sys
import argparse
from os.path import join
from os import remove
import numpy as np

from astropy.io import ascii

from find_data_dir import determine_data_dir
from get_fields import get_fields

def main(argv):
    """ @TODO main docstring
    """
    parser = argparse.ArgumentParser()
    
    # Data directory
    parser.add_argument("--data_dir",type=str, default=None,
                        help="The directory where CFHTLenS data can be found.")
    
    # Execute command-line parsing
    args = parser.parse_args()
    
    desired_data_dir = args.data_dir
    
    data_dir = determine_data_dir(desired_data_dir)
            
    print("Using " + data_dir + " as data directory.")
    
    field_names = get_fields(data_dir)
    
    # Loop through tables and add each to the output table
    for field_name in field_names:
        
        field_filename = join(data_dir,"full_tables",field_name + ".dat")
        
        field_table = ascii.read(field_filename)
        
        # Discard bad rows
        bad = np.logical_or(np.logical_or(field_table['Z_B']<0.2,field_table['fitclass']>=1),
                            field_table['CLASS_STAR']>0.7)
        field_table.remove_rows(np.where(bad))
        
        # Get some derived data
        field_table['e'] = np.sqrt(np.square(field_table['e1'])+np.square(field_table['e2']))
        field_table['sigma_arcsec'] = 3600.*field_table['FWHM_WORLD']/0.674490
        
        # Rename the y column i if it's present
        try:
            field_table.rename_column('MAG_y','MAG_i')
        except:
            pass
        
        # Discard columns we don't need
        field_table.keep_columns(['sigma_arcsec', 'MAG_i', 'Z_B', 'e'])
        
        # Output this field as a fits table
    
        output_filename = join(data_dir,"filtered_tables",field_name[0:-2]+"_mag_size.fits")
    
        try:
            remove(output_filename)
        except:
            pass
        field_table.write(output_filename,format="fits")
        
        print("Output " + output_filename + ".")
    
    print("Done!")
    
    return

if __name__ == "__main__":
    main(sys.argv)
