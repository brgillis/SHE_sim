""" @file /disk2/brg/git/SHE_sim/shear_measurement_correction/shear_calibration/calibrate_from_regression.py

    Created 28 Aug 2015

    Contains a function to call calibrate_all_results using the measured
    m and c data from a regression.

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

from astropy.io import fits

from common import magic_values as mv
from shear_calibration.calibrate_results import calibrate_all_results

def calibrate_from_regression(**kwargs):
    """ Calls calibrate_all_results using the measured m and c data from a regression.
    
        Required kwargs: path <string> (path where measurements are located)
                         reg_file <string> (filename of regression data file)
                         tag <string> (Extra label to add to calibrated results files.)
                         processes <int>
    """
    
    regression_results = fits.open(kwargs['reg_file'])[1].data
    
    kwargs['m1'] = regression_results[mv.result_m_colname][mv.result_comp_1_pix_row_index]
    kwargs['m2'] = regression_results[mv.result_m_colname][mv.result_comp_2_pix_row_index]
    kwargs['c1'] = regression_results[mv.result_c_colname][mv.result_comp_1_pix_row_index]
    kwargs['c2'] = regression_results[mv.result_c_colname][mv.result_comp_2_pix_row_index]
    
    kwargs['delta_m1'] = regression_results[mv.result_m_stderr_colname][mv.result_comp_1_pix_row_index]
    kwargs['delta_m2'] = regression_results[mv.result_m_stderr_colname][mv.result_comp_2_pix_row_index]
    kwargs['delta_c1'] = regression_results[mv.result_c_stderr_colname][mv.result_comp_1_pix_row_index]
    kwargs['delta_c2'] = regression_results[mv.result_c_stderr_colname][mv.result_comp_2_pix_row_index]
    
    calibrate_all_results(**kwargs)