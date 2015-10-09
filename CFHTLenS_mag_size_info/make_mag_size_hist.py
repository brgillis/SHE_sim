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
import numpy as np
from os.path import join

from astropy.io import fits

from find_data_dir import determine_data_dir
from get_fields import get_fields

# Magic values
num_mag_bins = 20
num_size_bins = 40

mag_range = [18.,24.7]
log_size_range = [-0.3,1.]

from magic_values import hist_image_filename

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
    print("Using " + data_dir + " as data directory.")
    
    # Init empty combined hist
    combined_hist2D = None
    
    # Loop through all tables in this directory
    field_names = get_fields(data_dir)
    
    for field_name in field_names:
        # Get the table name
        table_filename = join(data_dir,"filtered_tables",field_name[0:-2] + "_mag_size.fits")
    
        # Load in the data
        data = fits.open(table_filename)[1].data
        
        mags = data['MAG_i']
        sizes = data['sigma_arcsec']
        log_sizes = np.log10(sizes)
        
        # Get the histogram
        hist2D, _xedges, _yedges = np.histogram2d(mags, log_sizes,
                                                bins=[num_mag_bins,num_size_bins],
                                                range=[mag_range,log_size_range],
                                                normed=False)
        
        # Add to combined data (or overwrite if the first)
        if combined_hist2D is None:
            combined_hist2D = hist2D
        else:
            combined_hist2D += hist2D
            
        print("Added data from file " + table_filename + ".")
    
    # Normalize the histogram
    combined_hist2D /= combined_hist2D.sum()
    
    # Make it a fits image
    hist_hdu = fits.PrimaryHDU(combined_hist2D)
    
    # Add min/max data to the header
    hist_hdu.header['MAG_MIN'] = mag_range[0]
    hist_hdu.header['MAG_MAX'] = mag_range[1]
    hist_hdu.header['SIZE_MIN'] = log_size_range[0]
    hist_hdu.header['SIZE_MAX'] = log_size_range[1]
    
    # Save it
    hist_hdu.writeto(join(data_dir,hist_image_filename),clobber=True)
    
    print("Done!")

if __name__ == "__main__":
    main(sys.argv)
