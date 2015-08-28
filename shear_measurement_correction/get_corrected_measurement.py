#!/usr/bin/env python

""" @file get_corrected_measurement.py

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

from shear_calibration.calibrate_results import calibrate_all_results

@click.command()
@click.option("--path", default=".", help="Root path where shear measurement data is contained.")

@click.option("--m", default=0., help="Best-guess multiplicative bias parameter.")
@click.option("--c", default=0., help="Best-guess additive bias parameter.")
@click.option("--delta-m", "delta_m", default=0.,
              help="Error on estimate of multiplicative bias parameter.")
@click.option("--delta-c", "delta_c", default=0.,
              help="Error on estimate of additive bias parameter.")

@click.option("--tag", default="_calibrated", help="Extra label to add to calibrated results files.")

@click.option("--processes", default=mv.max_num_threads, help="Number of parallel processes to use.")
def main(**kwargs):
    """ Main function for generating corrected shear estimates. Run this script with the
        --help option to see available options.
    """
    
    calibrate_all_results(**kwargs)

if __name__ == "__main__":
    main() # click will handle arguments
