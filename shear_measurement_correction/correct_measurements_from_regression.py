#!/usr/bin/env python

""" @file correct_measurements_from_regression.py

    Created 13 Aug 2015

    A script for generating a corrected shear measurement based on the
    results of calibration. Run this script with the --help option to see
    available options.

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

import click

from common import magic_values as mv

from shear_calibration.calibrate_from_regression import calibrate_from_regression

@click.command()
@click.option("--path", default=".", help="Root path where shear measurement data is contained.")
@click.option("--reg_file",default="stored_data.fits",
              help="'stored_data' file which contains regression data.")

@click.option("--input_tag", default="calibrated", help="Required pattern for input files")
@click.option("--output_tag", default="calibrated", help="Extra label to add to calibrated results files.")

@click.option("--processes", default=mv.max_num_threads, help="Number of parallel processes to use.")
@click.option("--strict", default=False, help="Whether or not to require tag to be in the expected place.")
def main(**kwargs):
    """ Main function for generating corrected shear estimates. Run this script with the
        --help option to see available options.
    """
    
    if(kwargs['output_tag'] is None):
        kwargs['tag_'] = ""
    else:
        if(len(kwargs['output_tag'])==0):
            kwargs['tag_'] = ""
        else:
            kwargs['tag_'] = "_" + kwargs['output_tag']
    
    calibrate_from_regression(**kwargs)

if __name__ == "__main__":
    main() # click will handle arguments
