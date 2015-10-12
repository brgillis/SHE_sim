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
free_fit_filename = "mag_size_free_fits.fits"

def skewed_norm(x, a):
    return 2 * norm.pdf(x) * norm.cdf(a*x)
    
params_guess = [-8.46077932e-02, # 0: p_mean_log_slope
                8.05520159e-01, # 1: p_mean_log_intercept
                9.74303502e-02, # 2: p_stddev
                1.00897203e+00, # 3: p_skew
                -1.65903253e-01, # 4: ps_ratio_log_slope
                3.20201220e+00, # 5: ps_ratio_log_intercept
                5.88556617e+01, # 6: s_mean_pow
                2.25745045e+01, # 7: s_mean_knee
                -7.03205126e-02, # 8: s_mean_log_slope_low
                8.78148097e-01, # 9: s_mean_log_intercept_low
                -1.88560094e-01, # 10: s_mean_log_slope_high
                3.69534607e+00, # 11: s_mean_log_intercept_high
                4.41422174e+01, # 12: s_stddev_pow
                2.30583547e+01, # 13: s_stddev_knee
                -6.95018834e-01, # 14: s_stddev_log_low
                -3.64279452e-02, # 15: s_stddev_log_slope_high
                -6.86595310e-02, # 16: s_stddev_log_intercept_high
                1.54838862e+00, # 17: s_skew
                -1.03737338e+00, # 18: st_ratio_log_slope
                2.46766515e+01, # 19: st_ratio_log_intercept
                -1.70190582e-01, # 20: t_mean_log_slope
                4.54688588e+00, # 21: t_mean_log_intercept
                1.72665174e-01, # 22: t_stddev
                -3.53317861e+00] # 23: t_skew
    
params_lower = [-0.5, # 0: p_mean_log_slope
                0., # 1: p_mean_log_intercept
                0.05, # 2: p_stddev
                -10., # 3: p_skew
                -1., # 4: ps_ratio_log_slope
                0., # 5: ps_ratio_log_intercept
                1., # 6: s_mean_pow
                21., # 7: s_mean_knee
                -1., # 8: s_mean_log_slope_low
                0., # 9: s_mean_log_intercept_low
                -1., # 10: s_mean_log_slope_high
                0., # 11: s_mean_log_intercept_high
                1., # 12: s_stddev_pow
                21., # 13: s_stddev_knee
                -1., # 14: s_stddev_log_low
                -1., # 15: s_stddev_log_slope_high
                0., # 16: s_stddev_log_intercept_high
                -10., # 17: s_skew
                -2., # 18: st_ratio_log_slope
                0., # 19: st_ratio_log_intercept
                -1., # 20: t_mean_log_slope
                0., # 21: t_mean_log_intercept
                0.05, # 22: t_stddev
                -10] # 23: t_skew
    
params_upper = [0., # 0: p_mean_log_slope
                5., # 1: p_mean_log_intercept
                0.15, # 2: p_stddev
                10., # 3: p_skew
                0., # 4: ps_ratio_log_slope
                10., # 5: ps_ratio_log_intercept
                200., # 6: s_mean_pow
                24., # 7: s_mean_knee
                0., # 8: s_mean_log_slope_low
                10., # 9: s_mean_log_intercept_low
                0., # 10: s_mean_log_slope_high
                10., # 11: s_mean_log_intercept_high
                200., # 12: s_stddev_pow
                24., # 13: s_stddev_knee
                0., # 14: s_stddev_log_low
                0., # 15: s_stddev_log_slope_high
                10., # 16: s_stddev_log_intercept_high
                10., # 17: s_skew
                0., # 18: st_ratio_log_slope
                24., # 19: st_ratio_log_intercept
                0., # 20: t_mean_log_slope
                10., # 21: t_mean_log_intercept
                0.3, # 22: t_stddev
                10.] # 23: t_skew

def estimate_pdf_parameters(params, mag_mid):
    
    p_mean = 10.**(params[0]*mag_mid + params[1])
    p_stddev = params[2]*np.ones_like(mag_mid)
    p_skew = params[3]*np.ones_like(mag_mid)
    ps_ratio = 10.**(params[4]*mag_mid + params[5])
    
    s_mean_pow = params[6]*np.ones_like(mag_mid)
    xp = (mag_mid/params[7])**s_mean_pow
    s_mean = 10.**((params[8]*mag_mid + params[9])/(1+xp) + 
                     (params[10]*mag_mid + params[11])*xp/(1+xp))
    
    s_stddev_pow = params[12]*np.ones_like(mag_mid)
    xp = (mag_mid/params[13])**s_stddev_pow
    s_stddev = 10.**(params[14]/(1+xp) + 
                     (params[15]*mag_mid + params[16])*xp/(1+xp))
    s_skew = params[17]*np.ones_like(mag_mid)
    
    st_ratio = 10.**(params[18]*mag_mid + params[19])
    
    t_mean = (params[20]*mag_mid + params[21])
    t_stddev = params[22]*np.ones_like(mag_mid)
    t_skew = params[23]*np.ones_like(mag_mid)
        
    return p_mean, p_stddev, p_skew, ps_ratio, \
        s_mean, s_stddev, s_skew, st_ratio, \
        t_mean, t_stddev, t_skew
    

def estimate_size_pdf(sizes,
                      p_mean,
                      p_stddev,
                      p_skew,
                      ps_ratio,
                      s_mean,
                      s_stddev,
                      s_skew,
                      st_ratio,
                      t_mean,
                      t_stddev,
                      t_skew,
                      mag_mid=None,
                      fit_table=None,
                      fit_q=0.):
    
    p_xs = (sizes - p_mean) / p_stddev
    p_ps = skewed_norm(p_xs, p_skew) / p_stddev * (ps_ratio / (1. + ps_ratio))
    s_xs = (sizes - s_mean) / s_stddev
    s_ps = skewed_norm(s_xs, s_skew) / s_stddev * (1. / (1. + ps_ratio)) * (st_ratio / (1. + st_ratio))
    t_xs = (sizes - t_mean) / t_stddev
    t_ps = skewed_norm(t_xs, t_skew) / t_stddev * (1. / (1. + ps_ratio)) * (1. / (1. + st_ratio))
    p_predictions = p_ps + s_ps + t_ps
    
    # Add these params to the output table if desired
    if fit_table is not None:
        fit_table.add_row({"mag_mid":mag_mid,
                           "size_p_mean":p_mean,
                           "size_p_stddev":p_stddev,
                           "size_p_skew":p_skew,
                           "size_s_mean":s_mean,
                           "size_s_stddev":s_stddev,
                           "size_s_skew":s_skew,
                           "size_ps_ratio":ps_ratio,
                           "size_t_mean":t_mean,
                           "size_t_stddev":t_stddev,
                           "size_t_skew":t_skew,
                           "size_st_ratio":st_ratio,
                           "fit_q":fit_q})
    return p_predictions

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
    
    # Normalize each magnitude row of the histogram, and keep the data on the size of each row
    mag_row_sums = np.array([np.sum(mag_size_hist.data,axis=1)]).transpose()
    total_mag_row_sum = mag_row_sums.sum()
    normed_mag_size_hist = np.divide(mag_size_hist.data,mag_row_sums*size_binsize)
    
    
    # Get an array of sizes for the middles of each bin (keep size untransposed so it'll be orthogonal to mag axis)
    sizes = np.linspace(start=size_min,stop=size_max,num=size_nbins,endpoint=False) + \
                0.5*(size_max-size_min)/size_nbins
                
    # Get an array of mags for the middles of each bin
    mag_mids = np.array([np.linspace(start=mag_min,stop=mag_max,num=mag_nbins,endpoint=False) + \
                0.5*(mag_max-mag_min)/mag_nbins]).transpose()
                
    # Initialize output tables
    fit_table = Table(names=["mag_mid",
                             "size_p_mean",
                             "size_p_stddev",
                             "size_p_skew",
                             "size_s_mean",
                             "size_s_stddev",
                             "size_s_skew",
                             "size_ps_ratio",
                             "size_t_mean",
                             "size_t_stddev",
                             "size_t_skew",
                             "size_st_ratio",
                             "fit_q"])
    free_fit_table = Table(names=["mag_mid",
                             "size_p_mean",
                             "size_p_stddev",
                             "size_p_skew",
                             "size_s_mean",
                             "size_s_stddev",
                             "size_s_skew",
                             "size_ps_ratio",
                             "size_t_mean",
                             "size_t_stddev",
                             "size_t_skew",
                             "size_st_ratio",
                             "fit_q"])
    
    def get_uChi_squared_for_row(params,mag_index=None):
        
        p_predictions = estimate_size_pdf(sizes,
                  *params,
                  fit_table=None)
        
        # Check for zero p values
        p_predictions += 1e-9 # Impose floor
        
        if(mag_index is None):
            square_diffs = np.square((p_predictions - normed_mag_size_hist))
            uChi_squared =  np.sum(square_diffs/np.abs(p_predictions),axis=1)
        else:
            square_diffs = np.square((p_predictions - normed_mag_size_hist[mag_index]))
            uChi_squared =  np.sum(square_diffs/np.abs(p_predictions))
            

        return uChi_squared
    
    def get_total_uChi_squared(fit_params):
        
        row_params = estimate_pdf_parameters(fit_params, mag_mids)
        
        row_uChi_squareds = np.array([get_uChi_squared_for_row(row_params)]).transpose()
        
        total_uChi_squared = (row_uChi_squareds).sum()/total_mag_row_sum
        
        return total_uChi_squared
    
    # Fit the best set of parameters
    # best_fit_params = fmin(get_total_uChi_squared,params_guess,maxiter=50000,maxfun=100000)
    best_fit_params = params_guess
                             
    print("Guess cost: " + str(get_total_uChi_squared(params_guess)))
    print("fmin fit cost: " + str(get_total_uChi_squared(best_fit_params)))
    
    print("Best fit params: " + str(best_fit_params))
    
    # Loop through each magnitude bin, and plot fitted model for it
    for mag_index in xrange(mag_nbins):
        
        mag_mid = mag_mids[mag_index]
        
        row_params = np.array(estimate_pdf_parameters(best_fit_params, mag_mid)).transpose()[0]
        free_row_params = fmin(get_uChi_squared_for_row,row_params,
                               args=(mag_index,),maxfun=100000,maxiter=50000)
        
        fit_q = get_uChi_squared_for_row(row_params,mag_index)
        free_fit_q = get_uChi_squared_for_row(free_row_params,mag_index)
        
        p_vals = mag_size_hist.data[mag_index,:]/mag_size_hist.data[mag_index,:].sum()/size_binsize
        
        p_predictions = estimate_size_pdf(sizes, *row_params, fit_table=fit_table, fit_q=fit_q, mag_mid=mag_mid)
        free_p_predictions = estimate_size_pdf(sizes, *free_row_params, fit_table=free_fit_table,
                                               fit_q=free_fit_q, mag_mid=mag_mid)
        
        # Set up the figure
        fig = pyplot.figure()
        fig.subplots_adjust(wspace=0, hspace=0, bottom=0.1, right=0.95, top=0.95, left=0.12)
        
        ax = fig.add_subplot(1,1,1)
        ax.set_xlabel(x_label,fontsize=fontsize)
        ax.set_ylabel(y_label,fontsize=fontsize)
        
        # Plot the p values and predictions
        ax.plot(sizes,p_vals,color='k',marker='None',label="Actual")
        ax.plot(sizes,p_predictions,color='r',marker='None',label="Full model")
        ax.plot(sizes,free_p_predictions,color='b',marker='None',label="Free model")
        
        mag_label = "mag_" + str(mag_mid[0]+0.05)[0:4]
        
        ax.set_title("MAG\_i = " + mag_label[4:])
        
        plot_filename = join(data_dir,mag_size_subdir,mag_label.replace(".","_") + "_size_pdf.eps")
        
        #fig.show()
        
        pyplot.savefig(plot_filename, format="eps", bbox_inches="tight", pad_inches=0.05)
        
    # Output the fits tables
    output_filename = join(data_dir,mag_size_subdir,fit_filename)
    try:
        remove(output_filename)
    except:
        pass
    fit_table.write(output_filename,format="fits")
    output_filename = join(data_dir,mag_size_subdir,free_fit_filename)
    try:
        remove(output_filename)
    except:
        pass
    free_fit_table.write(output_filename,format="fits")
    
    print("Done!")
    
    return       

if __name__ == "__main__":
    main(sys.argv)
