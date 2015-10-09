""" @file get_fields.py

    Created 9 Oct 2015

    Function to get a list of field filenames

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

from os.path import join

from magic_values import fields_list_filename

def get_fields(data_dir):
    
    qualified_fields_list_filename = join(data_dir, fields_list_filename)
    
    field_names = []
    
    with open(qualified_fields_list_filename) as fields_list:
        for line in fields_list:
            for field_name in line.split():
                field_names.append(field_name)
    
    return field_names