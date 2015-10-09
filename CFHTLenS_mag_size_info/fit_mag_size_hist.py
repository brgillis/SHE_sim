#!/usr/bin/env python

""" @file /disk2/brg/git/CFHTLenS_cat/CFHTLenS_spec_matching/fit_mag_size_hist.py

    Created 8 Oct 2015

    Fits gaussians to the mag-size histogram. 

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
from os import remove
from astropy.io import fits
from astropy.table import Table
from scipy.stats import norm
from scipy.optimize import fmin

import matplotlib
from matplotlib import pyplot
matplotlib.rcParams['ps.useafm'] = True
matplotlib.rcParams['pdf.use14corefonts'] = True
matplotlib.rcParams['text.usetex'] = True

from find_data_dir import determine_data_dir

from magic_values import hist_image_filename, mag_size_subdir

# Magic values
x_label = r"$\log_{\rm 10}\left(\sigma/{\rm arcsec}\right)$"
y_label = "p"
fontsize = 12

fit_filename = "mag_size_fits.fits"

def skewed_norm(x, a):
    return 2 * norm.pdf(x) * norm.cdf(a*x)

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
    
    mag_size_hist_filename = join(data_dir,mag_size_subdir,hist_image_filename)
    
    mag_size_hist = fits.open(mag_size_hist_filename)[0]
    
    mag_min = mag_size_hist.header['MAG_MIN']
    mag_max = mag_size_hist.header['MAG_MAX']
    size_min = mag_size_hist.header['SIZE_MIN']
    size_max = mag_size_hist.header['SIZE_MAX']
    
    mag_nbins = mag_size_hist.header['NAXIS2'] # 2 and 1 swapped since fits uses opposite ordering
    size_nbins = mag_size_hist.header['NAXIS1']
    
    size_binsize = (size_max-size_min)/size_nbins
    
    # Get an array of sizes for the middles of each bin
    sizes = np.linspace(start=size_min,stop=size_max,num=size_nbins,endpoint=False) + \
                0.5*(size_max-size_min)/size_nbins
                
    # Initialize an output table
    fit_table = Table(names=["mag_mid",
                             "size_p_mean",
                             "size_p_stddev",
                             "size_p_skew",
                             "size_s_mean",
                             "size_s_stddev",
                             "size_s_skew",
                             "size_ps_ratio"])
    
    # Initial parameter guess
    params_guess = np.array([0.17, # Primary mean
                             0.11, # Primary stddev
                             0.67, # Primary Skew
                             0.33, # Secondary mean
                             0.18, # Secondary stddev
                             0.94, # Secondary Skew
                             1.3]) # Primary-secondary scale
    
    # Loop through each magnitude bin, and fit a model to the size data for it
    for mag_index in xrange(mag_nbins):
        
        mag_mid = mag_min + (mag_index+0.5)*(mag_max-mag_min)/mag_nbins
        
        p_vals = mag_size_hist.data[mag_index,:]/mag_size_hist.data[mag_index,:].sum()/size_binsize
        
        def get_uChi_squared(params):
            
            ratio = params[6]
            
            p_mean = params[0]
            p_stddev = params[1]
            p_skew = params[2]
            
            p_xs = (sizes-p_mean)/p_stddev
            
            p_ps = skewed_norm(p_xs, p_skew)/p_stddev * (ratio/(1.+ratio))
            
            s_mean = params[3]
            s_stddev = params[4]
            s_skew = params[5]
            
            s_xs = (sizes-s_mean)/s_stddev
            
            s_ps = skewed_norm(s_xs, s_skew)/s_stddev * (1./(1.+ratio))
            
            p_predictions = p_ps + s_ps
            
            square_diffs = np.square((p_predictions - p_vals))
            
            uChi_squared =  (square_diffs/p_predictions).sum()
            
            return uChi_squared
        
        
        params = fmin(get_uChi_squared,params_guess,maxiter=5000,maxfun=10000)
    
        ratio = params[6]
        
        p_mean = params[0]
        p_stddev = params[1]
        p_skew = params[2]
        
        p_xs = (sizes-p_mean)/p_stddev
        
        p_ps = skewed_norm(p_xs, p_skew)/p_stddev * (ratio/(1.+ratio))
        
        s_mean = params[3]
        s_stddev = params[4]
        s_skew = params[5]
        
        s_xs = (sizes-s_mean)/s_stddev
        
        s_ps = skewed_norm(s_xs, s_skew)/s_stddev * (1./(1.+ratio))
        
        p_predictions = p_ps + s_ps
        
        # Add these params to the output table
        fit_table.add_row({"mag_mid": mag_mid,
                           "size_p_mean": params[0],
                           "size_p_stddev": params[1],
                           "size_p_skew": params[2],
                           "size_s_mean": params[3],
                           "size_s_stddev": params[4],
                           "size_s_skew": params[5],
                           "size_ps_ratio": params[6]})
        # Setup the figure
        fig = pyplot.figure()
        fig.subplots_adjust(wspace=0, hspace=0, bottom=0.1, right=0.95, top=0.95, left=0.12)
        
        ax = fig.add_subplot(1,1,1)
        ax.set_xlabel(x_label,fontsize=fontsize)
        ax.set_ylabel(y_label,fontsize=fontsize)
        
        # Plot the p values and predictions
        ax.plot(sizes,p_vals,color='r',marker='None',label="Actual")
        ax.plot(sizes,p_predictions,color='k',marker='None',label="Model")
        
        mag_label = "mag_" + str(mag_mid+0.05)[0:4]
        
        ax.set_title("MAG\_i = " + mag_label[4:])
        
        plot_filename = join(data_dir,mag_size_subdir,mag_label.replace(".","_") + "_size_pdf.eps")
        
        pyplot.savefig(plot_filename, format="eps", bbox_inches="tight", pad_inches=0.05)
        fig.show()
        
    # Output the fits table
    output_filename = join(data_dir,mag_size_subdir,fit_filename)
    try:
        remove(output_filename)
    except:
        pass
    fit_table.write(output_filename,format="fits")
            

if __name__ == "__main__":
    main(sys.argv)
