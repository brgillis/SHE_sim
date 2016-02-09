#!/usr/bin/env python

""" @file /disk2/brg/git/SHE_sim/SHE_SIM_gal_params/python/determine_n_galaxy_fit.py

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
from astropy.io import ascii

matplotlib.rcParams['ps.useafm'] = True
matplotlib.rcParams['pdf.use14corefonts'] = True
matplotlib.rcParams['text.usetex'] = True

def main(argv):
    """ @TODO main docstring
    """

    # Data
    table = ascii.read("/disk2/brg/git/CFHTLenS_cat/Data/full_tables/W1m0m0_i.dat")
    gal_zs = table["Z_B"]
    field_size_sq_deg = 0.893539028901
    
    z_edges = np.linspace(0.3,2.0,18)
    zs = (z_edges[:-1] + z_edges[1:]) / 2.
    
    N_gals, _ = np.histogram(gal_zs,z_edges)
    
    field_size_sq_arcsec = field_size_sq_deg*np.power(60.,2)
    z_bin_width = z_edges[1] - z_edges[0]
    
    # Poisson error in number of clusters
    sigma_N_gals = np.sqrt(N_gals)
    
    # Normalize by field size
    n_gals = N_gals/field_size_sq_arcsec/z_bin_width
    sigma_n_gals = sigma_N_gals/field_size_sq_arcsec/z_bin_width
    
    # Determine guess values for the fit
    guess_n_scale = n_gals[0]/np.square(zs[0])
    guess_z_median = (zs[0]+zs[-1])/2
    guess = np.array((guess_n_scale,guess_z_median))
    
    def get_n_of_z(z,n_scale,z_median):
        return n_scale * np.square(z) * np.exp( -np.power(1.412*z/z_median,1.5) )
    
    def get_chi2_of_fit(scale_and_median):
        chi2 = np.sum( np.square( ( get_n_of_z(zs,scale_and_median[0],scale_and_median[1])
                                     - n_gals ) / sigma_n_gals ) )
        return chi2
    
    fit_n_scale, fit_z_median = fmin(get_chi2_of_fit,guess)
    
    chi2 = np.sum( np.square( ( fit_n_scale * np.square(zs) * np.exp( -np.power(1.412*zs/fit_z_median,1.5) )
                                     - n_gals ) / sigma_n_gals ) )
    
    print("n_scale = " + str(fit_n_scale) + " galaxies per arcmin^2 per unit redshift")
    print("z_median = " + str(fit_z_median))
    print("chi^2 = " + str(chi2))
    
    sample_zs = np.linspace(0.,2.,101)
    fit_ns = get_n_of_z(sample_zs,fit_n_scale,fit_z_median)
    
    fig = pyplot.figure()
    fig.subplots_adjust(wspace=0.5, hspace=0, bottom=0.1, right=0.95, top=0.95, left=0.12)
        
    ax = fig.add_subplot(1,1,1)
    ax.set_xlabel(r"$z$",labelpad=10,fontsize=16)
    ax.set_ylabel(r"$N/\mathrm{arcmin}^2/\mathrm{unit}\;\mathrm{redshift}$",labelpad=10,fontsize=16)
    
    ax.plot(sample_zs,fit_ns,color='r',label="Fit")
    
    ax.errorbar(zs,n_gals,sigma_n_gals,color='k',marker='.',label='CFHTLenS W1m0m0')
    
    ax.legend(loc='upper left')
    
    pyplot.savefig("n_galaxy_fit.png", format="png", bbox_inches="tight", pad_inches=0.05)
    
    pyplot.show()
    
    return

if __name__ == "__main__":
    main(sys.argv)
