""" @file calibrate_results.py

    Created 13 Aug 2015

    Function to calibrate the shear measurements in a single measurements
    file.

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

import astropy.io.fits as fits

from common import magic_values as mv
from common import SHE_fits_format as sff

from calibrate_shear import calibrate_shear, get_error_of_calibrated_shear

def calibrate_results(filename_tuple, **kwargs):
    
    results_filename = filename_tuple[mv.rf_tuple_index]
    
    # Check if there's actually a results file to process; if not, return
    if(results_filename==None):
        return
    
    # Load in the results table
    fits_result_HDUlist = fits.open(results_filename)
    fits_result_table = fits_result_HDUlist[1].data
    
    # Calibrate the shear and error columns
    fits_result_table[sff.gal_g1_colname] = calibrate_shear(fits_result_table[sff.gal_g1_colname],
                                                            m=kwargs['m'], c=kwargs['c'],
                                                            delta_m=kwargs['delta_m'],
                                                            delta_c=kwargs['delta_c'])
    fits_result_table[sff.gal_g2_colname] = calibrate_shear(fits_result_table[sff.gal_g2_colname],
                                                            m=kwargs['m'], c=kwargs['c'],
                                                            delta_m=kwargs['delta_m'],
                                                            delta_c=kwargs['delta_c'])
    fits_result_table[sff.gal_g1_err_colname] = get_error_of_calibrated_shear(
                                                            fits_result_table[sff.gal_g1_err_colname],
                                                            m=kwargs['m'], c=kwargs['c'],
                                                            delta_m=kwargs['delta_m'],
                                                            delta_c=kwargs['delta_c'])
    fits_result_table[sff.gal_g1_err_colname] = get_error_of_calibrated_shear(
                                                            fits_result_table[sff.gal_g1_err_colname],
                                                            m=kwargs['m'], c=kwargs['c'],
                                                            delta_m=kwargs['delta_m'],
                                                            delta_c=kwargs['delta_c'])
    
    # Output the calibrated shears to a new file
    calibrated_results_filename = results_filename.replace(mv.fits_table_extension,
                                                           kwargs['tag'] + mv.fits_table_extension)
    fits_result_HDUlist.writeto(calibrated_results_filename, clobber=True)
    
    return