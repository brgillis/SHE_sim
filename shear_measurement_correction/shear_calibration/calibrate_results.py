""" @file calibrate_results.py

    Created 13 Aug 2015

    Functions to calibrate the shear measurements in a single measurements
    file or in all measurements files.

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

from SHE_SBE_common.get_filenames import get_filenames
from SHE_SBE_common.parmap import parmap
from SHE_SBE_common import magic_values as mv
from SHE_SBE_common import SHE_fits_format as sff

from calibrate_shear import calibrate_shear, get_error_of_calibrated_shear

def calibrate_results(filename_tuple, **kwargs):
    """ Calibrate the results in a single measurements file.
    
        Required kwargs: m1 <float>
                         m2 <float>
                         c1 <float>
                         c2 <float>
                         delta_m1 <float>
                         delta_m2 <float>
                         delta_c1 <float>
                         delta_c2 <float>
                         tag_ <string>
    """
    
    results_filename = filename_tuple[mv.rf_tuple_index]
    
    # Check if there's actually a results file to process; if not, return
    if(results_filename==None):
        return
    
    # Load in the results table
    fits_result_HDUlist = fits.open(results_filename)
    fits_result_table = fits_result_HDUlist[1].data
    
    # Calibrate the shear and error columns
    fits_result_table[sff.gal_g1_colname] = calibrate_shear(fits_result_table[sff.gal_g1_colname],
                                                            m=kwargs['m1'], c=kwargs['c1'],
                                                            delta_m=kwargs['delta_m1'],
                                                            delta_c=kwargs['delta_c1'])
    fits_result_table[sff.gal_g2_colname] = calibrate_shear(fits_result_table[sff.gal_g2_colname],
                                                            m=kwargs['m2'], c=kwargs['c2'],
                                                            delta_m=kwargs['delta_m2'],
                                                            delta_c=kwargs['delta_c2'])
    fits_result_table[sff.gal_g1_err_colname] = get_error_of_calibrated_shear(
                                                            fits_result_table[sff.gal_g1_err_colname],
                                                            m=kwargs['m1'], c=kwargs['c1'],
                                                            delta_m=kwargs['delta_m1'],
                                                            delta_c=kwargs['delta_c1'])
    fits_result_table[sff.gal_g1_err_colname] = get_error_of_calibrated_shear(
                                                            fits_result_table[sff.gal_g1_err_colname],
                                                            m=kwargs['m2'], c=kwargs['c2'],
                                                            delta_m=kwargs['delta_m2'],
                                                            delta_c=kwargs['delta_c2'])
    
    # Output the calibrated shears to a new file
    calibrated_results_filename = results_filename.replace(mv.result_tail,
                                                           kwargs['tag_'] + mv.result_tail)
    fits_result_HDUlist.writeto(calibrated_results_filename, clobber=True)
    
    return

def calibrate_all_results(**kwargs):
    """ Calibrate the results in a single measurements file.
    
        Required kwargs: path <string>
                         m1 <float>
                         m2 <float>
                         c1 <float>
                         c2 <float>
                         delta_m1 <float>
                         delta_m2 <float>
                         delta_c1 <float>
                         delta_c2 <float>
                         tag_ <string>
                         processes <int>
                         strict <bool>
    """
    
    # Get the filenames in the path
    filename_tuples = get_filenames(kwargs['path'],output_pattern=kwargs['input_tag'],
                                    strict=kwargs['strict'])
    
    # Define a wrapper function for calibrating results
    def calibrate_results_wrapper(filename_tuple):
        return calibrate_results(filename_tuple, **kwargs)
    
    # Do the shape estimation in parallel if more than one file
    if((len(filename_tuples)>1) and (kwargs['processes']>1)):
        parmap(calibrate_results_wrapper, filename_tuples, nprocs=kwargs['processes'])
    else:
        for filename_tuple in filename_tuples:
            calibrate_results_wrapper(filename_tuple)

    print("Finished calibrating shear measurements.")