#!/usr/bin/env python

""" @file /disk2/brg/git/SHE_sim/SHE_SIM_gal_params/python/determine_n_cluster_fit.py

    Created 9 Feb 2016

    @TODO: File docstring

    ---------------------------------------------------------------------

    Copyright (C) 2016 Bryan R. Gillis

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

import matplotlib
from matplotlib import pyplot
import numpy as np
from scipy.optimize import fmin
import sys

matplotlib.rcParams['ps.useafm'] = True
matplotlib.rcParams['pdf.use14corefonts'] = True
matplotlib.rcParams['text.usetex'] = True

def main(argv):
    """ @TODO main docstring
    """

    # Data
    N_clusters = np.array((1161,1521,2248,2935,2456,2331,2364))
    zs = np.array((0.2,0.3,0.4,0.5,0.6,0.7,0.8))
    field_size_sq_deg = 154
    field_size_sq_arcsec = field_size_sq_deg*np.power(60.,4)
    z_bin_width = 0.1
    
    # Poisson error in number of clusters
    sigma_N_clusters = np.sqrt(N_clusters)
    
    # Normalize by field size
    n_clusters = N_clusters/field_size_sq_arcsec/z_bin_width
    sigma_n_clusters = sigma_N_clusters/field_size_sq_arcsec/z_bin_width
    
    # Determine guess values for the fit
    guess_n_scale = n_clusters[0]/np.square(zs[0])
    guess_z_median = (zs[0]+zs[-1])/2
    guess = np.array((guess_n_scale,guess_z_median))
    
    def get_n_of_z(z,n_scale,z_median):
        return n_scale * np.square(z) * np.exp( -np.power(1.412*z/z_median,1.5) )
    
    def get_chi2_of_fit(scale_and_median):
        chi2 = np.sum( np.square( ( get_n_of_z(zs,scale_and_median[0],scale_and_median[1])
                                     - n_clusters ) / sigma_n_clusters ) )
        return chi2
    
    fit_n_scale, fit_z_median = fmin(get_chi2_of_fit,guess)
    
    chi2 = np.sum( np.square( ( fit_n_scale * np.square(zs) * np.exp( -np.power(1.412*zs/fit_z_median,1.5) )
                                     - n_clusters ) / sigma_n_clusters ) )
    
    print("n_scale = " + str(fit_n_scale) + " clusters per arcsec^2 per unit redshift")
    print("z_median = " + str(fit_z_median))
    print("chi^2 = " + str(chi2))
    
    sample_zs = np.linspace(0.,1.,101)
    fit_ns = get_n_of_z(sample_zs,fit_n_scale,fit_z_median)
    
    fig = pyplot.figure()
    fig.subplots_adjust(wspace=0.5, hspace=0, bottom=0.1, right=0.95, top=0.95, left=0.12)
        
    ax = fig.add_subplot(1,1,1)
    ax.set_xlabel("z",labelpad=10,fontsize=16)
    ax.set_ylabel("N",labelpad=10,fontsize=16)
    
    ax.plot(sample_zs,fit_ns,color='r')
    
    ax.errorbar(zs,n_clusters,sigma_n_clusters,color='k',marker='.')
    
    pyplot.show()

if __name__ == "__main__":
    main(sys.argv)
