#!/usr/bin/env python

""" @file plot_mock_bias.py

    Created 27 Oct 2015

    Generates mock bias plots for presentations

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

import numpy as np
import sys

import matplotlib
import matplotlib.pyplot as pyplot
matplotlib.rcParams['ps.useafm'] = True
matplotlib.rcParams['pdf.use14corefonts'] = True
matplotlib.rcParams['text.usetex'] = True

def main(argv):
    
    figsize = (4,3)
    sigma_label = r"$\log_{10}\left( \sigma/arcsec \right)$"
    m_label = r"$m$"
    count_label = "Fraction of galaxies"
    fontsize = 12
    
    output_format = "eps"
    file_ext = "." + output_format
    
    sigma_bounds = [-0.85,-0.15]
    sigma_bounds_bar = [-0.9,-0.1]
    m_bounds = [-0.02, 0.062]
    
    sigma = np.linspace(start=-0.8,stop=-0.2,num=7)
    count_1 = np.array([100.,200.,100.,50.,25.,10.,5.])
    count_2 = np.array([50.,100.,200.,100.,50.,25.,10.])
    m = np.array([-0.01,-0.015,-0.011,-0.001,0.01,0.025,0.05])
    
    count_1 /= np.sum(count_1)
    count_2 /= np.sum(count_2)
    
    m_mean_1 = np.sum(m*count_1)
    m_mean_2 = np.sum(m*count_2)
    
    m_corr_1 = (1+m)*(1-m_mean_1+np.square(m_mean_1))-1
    
    m_corr_1_mean_1 = np.sum(m_corr_1*count_1)
    m_corr_1_mean_2 = np.sum(m_corr_1*count_2)
        
    # count plots
    for counts, bar_color, figname in \
        ((count_1,"r","mock_counts_1"),
         (count_2,"m","mock_counts_2"),
         ):
        
        # Setup the figure
        fig = pyplot.figure(figsize=figsize)
        fig.subplots_adjust(wspace=0, hspace=0, bottom=0.1, right=0.95, top=0.95, left=0.12)
        
        ax = fig.add_subplot(1,1,1)
        ax.set_xlabel(sigma_label,fontsize=fontsize)
        ax.set_ylabel(count_label,fontsize=fontsize)
        
        ax.set_xlim(sigma_bounds)
        #ax.set_ylim([0,1.1*np.max(counts)])
        
        # Plot stuff
        ax.bar(sigma,counts,color=bar_color)
        
        filename = figname + file_ext
        
        # Save the figure
        pyplot.savefig(filename, format=output_format, bbox_inches="tight", pad_inches=0.05)
        
        print("Finished plotting to " + filename + ".")
        
    # count plots
    for counts, bar_color, figname in \
        ((count_1,"r","mock_counts_1"),
         (count_2,"m","mock_counts_2"),
         ):
        
        # Setup the figure
        fig = pyplot.figure(figsize=figsize)
        fig.subplots_adjust(wspace=0, hspace=0, bottom=0.1, right=0.95, top=0.95, left=0.12)
        
        ax = fig.add_subplot(1,1,1)
        ax.set_xlabel(sigma_label,fontsize=fontsize)
        ax.set_ylabel(count_label,fontsize=fontsize)
        
        ax.set_xlim(sigma_bounds_bar)
        #ax.set_ylim([0,1.1*np.max(counts)])
        
        # Plot stuff
        ax.bar(sigma-0.05,counts,color=bar_color,width=0.1)
        
        filename = figname + file_ext
        
        # Save the figure
        pyplot.savefig(filename, format=output_format, bbox_inches="tight", pad_inches=0.05)
        
        print("Finished plotting to " + filename + ".")
    
    # m plots
    for m_vals, m_mean, m_mean_label, m_mean_color, m_vals_label, figname in \
        ((m,m_mean_1,r"$\left\langle m \right\rangle $", "r", r"$m$", "mock_bias_c1"),
         (m_corr_1,m_corr_1_mean_1,r"$\left\langle m \right\rangle $", "r", r"$m_{\rm corrected}$", "mock_bias_c1_corr1"),
         (m,m_mean_2,r"$\left\langle m \right\rangle ' $", "m", r"$m$", "mock_bias_c2"),
         (m_corr_1,m_corr_1_mean_2,r"$\left\langle m \right\rangle ' $", "m", r"$m_{\rm corrected}$", "mock_bias_c2_corr1"),
         ):
        
        # Setup the figure
        fig = pyplot.figure(figsize=figsize)
        fig.subplots_adjust(wspace=0, hspace=0, bottom=0.1, right=0.95, top=0.95, left=0.12)
        
        ax = fig.add_subplot(1,1,1)
        ax.set_xlabel(sigma_label,fontsize=fontsize)
        ax.set_ylabel(m_label,fontsize=fontsize)
        sigma_bounds_bar
        ax.set_xlim(sigma_bounds)
        ax.set_ylim(m_bounds)
        
        # Plot stuff
        ax.plot(sigma_bounds,[0,0],color='k',linestyle='solid',linewidth=0.5,label=None)

        ax.plot(sigma,m_vals,color='r',marker='o',linestyle='None',label=m_vals_label)
        
        ax.plot(sigma_bounds,[m_mean,m_mean],color=m_mean_color,linestyle='dashed',label=m_mean_label)
        
        ax.legend(loc="upper left",numpoints=1)
        
        filename = figname + file_ext
        
        # Save the figure
        pyplot.savefig(filename, format=output_format, bbox_inches="tight", pad_inches=0.05)
        
        print("Finished plotting to " + filename + ".")
    
    return

if __name__ == "__main__":
    main(sys.argv)
