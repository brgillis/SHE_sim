""" @file make_mag_size_hist.py

    Created 8 Oct 2015

    Generate a 2D histogram of magnitude and size

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
import numpy as np

from astropy.io import fits

from find_data_dir import determine_data_dir

# Magic values
table_filename = "filtered_tables/W_all_mag_size.fits"
num_mag_bins = 20
num_size_bins = 40

mag_range = [15.,25.]
size_range = [-2.,1.]

hist_image_filename = "size_mag_hist.fits"

def main(argv):
    parser = argparse.ArgumentParser()
    
    # Data directory
    parser.add_argument("--data_dir",type=str, default=None,
                        help="The directory where CFHTLenS data can be found.")
    
    # Execute command-line parsing
    args = parser.parse_args()

    # Get the data directory to use
    desired_data_dir = args.data_dir
    data_dir = determine_data_dir(desired_data_dir)
    
    # Open the data table and get its data
    qualified_table_filename = join(data_dir,table_filename)
    data = fits.open(qualified_table_filename).data
    
    mags = data['MAG_i']
    sizes = data['sigma_arcsec']
    log_sizes = np.log10(sizes)
    
    hist2D, _xedges, _yedges = np.histogram2d(mags, log_sizes,
                                            bins=[num_mag_bins,num_size_bins],
                                            range=[mag_range,size_range],
                                            normed=False)
    
    # Normalize the histogram
    hist2D /= hist2D.sum()
    
    # Output it as a fits image
    hist_hdu = fits.PrimaryHDU(hist2D)
    
    hist_hdu.writeto(hist_image_filename)

if __name__ == "__main__":
    main(sys.argv)
