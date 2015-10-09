""" @file find_data_dir.py

    Created 8 Oct 2015

    Function to determine where the CFHTLenS data is stored.

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

from os.path import join, isfile

from magic_values import primary_data_dir, secondary_data_dir, fields_list_filename

def determine_data_dir(data_dir):
    
    # Check that the data directory is good
    if (data_dir is not None):
        qualified_fields_list_filename = join(data_dir, fields_list_filename)
        if (not isfile(qualified_fields_list_filename)):
            print "WARNING: Cannot find fields list in supplied data path. Checking default paths."
            data_dir = None
            
    if (data_dir is None):
        # Try primary path
        qualified_fields_list_filename = join(primary_data_dir, fields_list_filename)
        if (isfile(qualified_fields_list_filename)):
            data_dir = primary_data_dir
        else:
            qualified_fields_list_filename = join(secondary_data_dir, fields_list_filename)
            if (isfile(qualified_fields_list_filename)):
                data_dir = secondary_data_dir
            else:
                print "ERROR: Cannot find a usable data path. Try passing the path to the CFHTLenS " + "data at the command-line with --data_dir=/path/to/data/" # Try secondary path
    return data_dir
