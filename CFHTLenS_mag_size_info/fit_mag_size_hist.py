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
from copy import deepcopy

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
fit_with_penalty_filename = "mag_size_fits_with_penalty.fits"
free_fit_filename = "mag_size_free_fits.fits"

def skewed_norm(x, a):
    return 2 * norm.pdf(x) * norm.cdf(a*x)

def skewed_norm_mean(loc, scale, skew):
    return loc + scale*np.sqrt(2./np.pi)*skew/np.sqrt(1.+np.square(skew))
    
params_guess = [-8.15189167e-02, # 0: p_mean_log_slope
                7.30317554e-01, # 1: p_mean_log_intercept
                1.03156870e-01, # 2: p_stddev
                1.24490509e+00, # 3: p_skew
                -1.90804758e-01, # 4: ps_ratio_log_slope
                3.67507714e+00, # 5: ps_ratio_log_intercept
                4.58263408e+01, # 6: s_mean_pow
                2.28830774e+01, # 7: s_mean_knee
                -7.50196950e-02, # 8: s_mean_log_slope_low
                9.70583289e-01, # 9: s_mean_log_intercept_low
                -2.37784554e-01, # 10: s_mean_log_slope_high
                4.88726191e+00, # 11: s_mean_log_intercept_high
                5.13142587e+01, # 12: s_stddev_pow
                2.29953184e+01, # 13: s_stddev_knee
                -6.96087917e-01, # 14: s_stddev_log_low
                -3.82344399e-02, # 15: s_stddev_log_slope_high
                -2.57460848e-02, # 16: s_stddev_log_intercept_high
                1.53110775e+00, # 17: s_skew
                -1.01327579e+00, # 18: st_ratio_log_slope
                2.41047904e+01, # 19: st_ratio_log_intercept
                -1.70988798e-01, # 20: t_mean_log_slope
                4.56528741e+00, # 21: t_mean_log_intercept
                1.68014695e-01, # 22: t_stddev
                -3.37760457e+00] # 23: t_skew
    
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
                      p_loc,
                      p_scale,
                      p_skew,
                      ps_ratio,
                      s_loc,
                      s_scale,
                      s_skew,
                      st_ratio,
                      t_loc,
                      t_scale,
                      t_skew,
                      mag_mid=None,
                      fit_table=None,
                      fit_q=0.):
    
    p_xs = (sizes - p_loc) / p_scale
    p_ps = skewed_norm(p_xs, p_skew) / p_scale * (ps_ratio / (1. + ps_ratio))
    s_xs = (sizes - s_loc) / s_scale
    s_ps = skewed_norm(s_xs, s_skew) / s_scale * (1. / (1. + ps_ratio)) * (st_ratio / (1. + st_ratio))
    t_xs = (sizes - t_loc) / t_scale
    t_ps = skewed_norm(t_xs, t_skew) / t_scale * (1. / (1. + ps_ratio)) * (1. / (1. + st_ratio))
    p_predictions = p_ps + s_ps + t_ps
    
    p_mean = skewed_norm_mean(p_loc,p_scale,p_skew)
    s_mean = skewed_norm_mean(s_loc,s_scale,s_skew)
    t_mean = skewed_norm_mean(t_loc,t_scale,t_skew)
    
    p_fraction = ps_ratio/(1.+ps_ratio)
    st_fraction = 1.-p_fraction
    s_fraction = st_fraction*st_ratio/(1.+st_ratio)
    t_fraction = 1.-p_fraction-s_fraction
    
    st_mean = (s_fraction*s_mean + t_fraction*t_mean)/st_fraction
    
    # Add these params to the output table if desired
    if fit_table is not None:
        fit_table.add_row({"mag_mid":mag_mid,
                            "size_p_loc":p_loc,
                            "size_p_scale":p_scale,
                            "size_p_skew":p_skew,
                            "size_p_mean":p_mean,
                            "size_ps_ratio":ps_ratio,
                            "size_p_fraction":p_fraction,
                            "size_s_loc":s_loc,
                            "size_s_scale":s_scale,
                            "size_s_skew":s_skew,
                            "size_s_mean":s_mean,
                            "size_t_loc":t_loc,
                            "size_t_scale":t_scale,
                            "size_t_skew":t_skew,
                            "size_t_mean":t_mean,
                            "size_st_ratio":st_ratio,
                            "size_st_fraction":st_fraction,
                            "size_st_mean":st_mean,
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
                
    # Initialize free fit parameter vectors
    free_p_means = np.zeros((mag_nbins,1))
    free_st_means = np.zeros((mag_nbins,1))
    free_p_fractions = np.zeros((mag_nbins,1))
    free_means = np.zeros((mag_nbins,1))
                
    # Initialize output tables
    fit_table = Table(names=["mag_mid",
                             "size_p_loc",
                             "size_p_scale",
                             "size_p_skew",
                             "size_p_mean",
                             "size_ps_ratio",
                             "size_p_fraction",
                             "size_s_loc",
                             "size_s_scale",
                             "size_s_skew",
                             "size_s_mean",
                             "size_t_loc",
                             "size_t_scale",
                             "size_t_skew",
                             "size_t_mean",
                             "size_st_ratio",
                             "size_st_fraction",
                             "size_st_mean",
                             "fit_q"])
    fit_table_with_penalty = deepcopy(fit_table)
    free_fit_table = deepcopy(fit_table)
    
    def get_uChi_squared_for_row(params,
                                 mag_index=None,
                                 p_mean_penalty=0.,
                                 st_mean_penalty=0.,
                                 p_fraction_penalty=0.,
                                 mean_penalty=0.):
        
        p_predictions = estimate_size_pdf(sizes,
                  *params,
                  fit_table=None)
        
        # Check for zero p values
        p_predictions += 1e-9 # Impose floor
            
        # Determine the means of this estimate        
        p_mean = skewed_norm_mean(params[0],params[1],params[2])
        s_mean = skewed_norm_mean(params[4],params[5],params[6])
        t_mean = skewed_norm_mean(params[8],params[9],params[10])
        ps_ratio = params[3]
        st_ratio = params[7]
        
        p_fraction = ps_ratio/(1.+ps_ratio)
        s_fraction = (1.-p_fraction)*st_ratio/(1.+st_ratio)
        t_fraction = (1.-p_fraction-s_fraction)
        
        st_fraction = s_fraction + t_fraction
        st_mean = (s_fraction*s_mean + t_fraction*t_mean)/st_fraction
        mean = p_fraction*p_mean + st_fraction*st_mean
        
        if(mag_index is None):
            square_diffs = np.square((p_predictions - normed_mag_size_hist))
            p_mean_offset_penalty = (p_mean_penalty*np.abs(p_mean-free_p_means)).transpose()[0]
            st_mean_offset_penalty = (st_mean_penalty*np.abs(st_mean-free_st_means)).transpose()[0]
            p_fraction_offset_penalty = (p_fraction_penalty*np.abs(p_fraction-free_p_fractions)).transpose()[0]
            mean_offset_penalty = (mean_penalty*np.abs(mean-free_means)).transpose()[0]
            
            axis = 1
        else:
            square_diffs = np.square((p_predictions - normed_mag_size_hist[mag_index]))
            p_mean_offset_penalty = p_mean_penalty*np.abs(p_mean-free_p_means[mag_index])
            st_mean_offset_penalty = st_mean_penalty*np.abs(st_mean-free_st_means[mag_index])
            p_fraction_offset_penalty = p_fraction_penalty*np.abs(p_fraction-free_p_fractions[mag_index])
            mean_offset_penalty = mean_penalty*np.abs(mean-free_means[mag_index])
            
            axis = None
            
        uChi_squared =  np.sum(square_diffs/np.abs(p_predictions),axis=axis)

        return uChi_squared + p_fraction*p_mean_offset_penalty + st_fraction*st_mean_offset_penalty + \
            p_fraction_offset_penalty + mean_offset_penalty
    
    def get_total_uChi_squared(fit_params, *params):
        
        row_params = estimate_pdf_parameters(fit_params, mag_mids)
        
        row_uChi_squareds = np.array([get_uChi_squared_for_row(row_params,*params)]).transpose()
        
        total_uChi_squared = (row_uChi_squareds*mag_row_sums).sum()/total_mag_row_sum
        
        return total_uChi_squared
    
    # Fit the best set of parameters without testing against mean/fraction offsets
    # best_fit_params = params_guess
    best_fit_params = fmin(get_total_uChi_squared,params_guess,maxiter=50000,maxfun=100000)
    
    free_row_params_list = []
    
    # Loop through each magnitude bin, and plot fitted model for it
    for mag_index in xrange(mag_nbins):
        
        mag_mid = mag_mids[mag_index]
        
        row_params = np.array(estimate_pdf_parameters(best_fit_params, mag_mid)).transpose()[0]
        free_row_params = fmin(get_uChi_squared_for_row,row_params,
                               args=(mag_index,),maxfun=100000,maxiter=50000)
        
        free_row_params_list.append(free_row_params)
        
        # Append these results to the free parameter arrays
        free_p_mean = free_row_params[0]
        free_s_mean = free_row_params[4]
        free_t_mean = free_row_params[8]
        free_ps_ratio = free_row_params[3]
        free_st_ratio = free_row_params[7]
        
        free_p_fraction = free_ps_ratio/(1.+free_ps_ratio)
        free_s_fraction = (1.-free_p_fraction)*free_st_ratio/(1.+free_st_ratio)
        free_t_fraction = (1.-free_p_fraction-free_s_fraction)
        
        free_st_fraction = free_s_fraction + free_t_fraction
        free_st_mean = (free_s_fraction*free_s_mean + free_t_fraction*free_t_mean)/free_st_fraction
        free_mean = free_p_fraction*free_p_mean + free_st_fraction*free_st_mean
        
        free_p_means[mag_index] = free_p_mean
        free_st_means[mag_index] = free_st_mean
        free_p_fractions[mag_index] = free_p_fraction
        free_means[mag_index] = free_mean
        
    # Now that we have the means and fractions from the free fits, fit the overall model again,
    # now with a penalty for offsets in these
    
    penalty_params = (None,0.5,0.5,1.,3.)
    
    best_fit_params_with_penalty = fmin(get_total_uChi_squared,params_guess,maxiter=50000,maxfun=100000,
                                        args=penalty_params)
    
    print("Guess cost w/o penalty: " + str(get_total_uChi_squared(params_guess)))
    print("Initial fit cost w/o penalty: " + str(get_total_uChi_squared(best_fit_params)))
    print("Second fit cost w/o penalty: " + str(get_total_uChi_squared(best_fit_params_with_penalty)))
                             
    print("Guess cost with penalty: " + str(get_total_uChi_squared(params_guess,*penalty_params)))
    print("Initial fit cost with penalty: " + str(get_total_uChi_squared(best_fit_params,*penalty_params)))
    print("Second fit cost with penalty: " + str(get_total_uChi_squared(best_fit_params_with_penalty,*penalty_params)))
    
    print("Initial fit params: " + str(best_fit_params))
    print("Second fit params: " + str(best_fit_params_with_penalty))
    
        
    # Now go through and make plots
    
    for mag_index in xrange(mag_nbins):
        
        mag_mid = mag_mids[mag_index]
        
        row_params = np.array(estimate_pdf_parameters(best_fit_params, mag_mid)).transpose()[0]
        row_params_with_penalty = np.array(estimate_pdf_parameters(best_fit_params_with_penalty,
                                                                   mag_mid)).transpose()[0]
        free_row_params = free_row_params_list[mag_index]
        
        fit_q = get_uChi_squared_for_row(row_params,mag_index)
        fit_q_with_penalty = get_uChi_squared_for_row(row_params,*penalty_params)
        free_fit_q = get_uChi_squared_for_row(free_row_params,mag_index)
        
        p_vals = mag_size_hist.data[mag_index,:]/mag_size_hist.data[mag_index,:].sum()/size_binsize
        
        p_predictions = estimate_size_pdf(sizes, *row_params, fit_table=fit_table, fit_q=fit_q,
                                          mag_mid=mag_mid)
        p_predictions_with_penalty = estimate_size_pdf(sizes, *row_params_with_penalty,
                                                       fit_table=fit_table_with_penalty, fit_q=fit_q_with_penalty,
                                                       mag_mid=mag_mid)
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
        ax.plot(sizes,p_predictions_with_penalty,color='m',marker='None',label="Full model with penalty")
        ax.plot(sizes,free_p_predictions,color='b',marker='None',label="Free model")
        
        mag_label = "mag_" + str(mag_mid[0]+0.05)[0:4]
        
        ax.set_title("MAG\_i = " + mag_label[4:])
        
        ax.legend(numpoints=1,loc='upper right')
        
        fig.show()
        
        plot_filename = join(data_dir,mag_size_subdir,mag_label.replace(".","_") + "_size_pdf.eps")
        pyplot.savefig(plot_filename, format="eps", bbox_inches="tight", pad_inches=0.05)
        
    # Output the fits tables
    output_filename = join(data_dir,mag_size_subdir,fit_filename)
    try:
        remove(output_filename)
    except:
        pass
    fit_table.write(output_filename,format="fits")
    
    output_filename = join(data_dir,mag_size_subdir,fit_with_penalty_filename)
    try:
        remove(output_filename)
    except:
        pass
    fit_table_with_penalty.write(output_filename,format="fits")
    
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
