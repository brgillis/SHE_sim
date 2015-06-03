'''
Created on 2 Jun 2015

@author: brg
'''

import numpy as np
import sys

from astropy.io import fits

import matplotlib
import matplotlib.pyplot as pyplot
matplotlib.rcParams['ps.useafm'] = True
matplotlib.rcParams['pdf.use14corefonts'] = True
matplotlib.rcParams['text.usetex'] = True


def main(argv):
    
    team_names = ("Control", "1", "2", "3")
    
    low_SN_data_files = ("/home/brg/git/SHE_sim/Generate_GalSim_Images/trunk/benchmark_low_SN/stored_data.fits",
                          "/home/brg/Data/SBE/MegaLUT/low_SN/stored_data.fits",
                          "/home/brg/Data/SBE/Oxford/low_SN/stored_data.fits",
                          "/home/brg/Data/SBE/gfit/high_SN/stored_data.fits")
    high_SN_data_files = ("/home/brg/git/SHE_sim/Generate_GalSim_Images/trunk/benchmark_high_SN/stored_data.fits",
                          "/home/brg/Data/SBE/MegaLUT/high_SN/stored_data.fits",
                          "/home/brg/Data/SBE/Oxford/high_SN/stored_data.fits",
                          "/home/brg/Data/SBE/gfit/high_SN/stored_data.fits")
    
    num_teams = len(team_names)
    m_limits = (0.001,0.3)
    c_limits = (0.5e-5,1.5e-3)
    
    m_target = 2e-3
    c_target = 1e-5
    
    fontsize = 24
    SN_shift = 0.05
    
    # Results
    
    low_SN_m1_pix = []
    low_SN_m2_pix = []
    low_SN_c1_pix = []
    low_SN_c2_pix = []
    low_SN_m1_stderr_pix = []
    low_SN_m2_stderr_pix = []
    low_SN_c1_stderr_pix = []
    low_SN_c2_stderr_pix = []
    low_SN_m1_psf = []
    low_SN_m2_psf = []
    low_SN_c1_psf = []
    low_SN_c2_psf = []
    low_SN_m1_stderr_psf = []
    low_SN_m2_stderr_psf = []
    low_SN_c1_stderr_psf = []
    low_SN_c2_stderr_psf = []
    
    high_SN_m1_pix = []
    high_SN_m2_pix = []
    high_SN_c1_pix = []
    high_SN_c2_pix = []
    high_SN_m1_stderr_pix = []
    high_SN_m2_stderr_pix = []
    high_SN_c1_stderr_pix = []
    high_SN_c2_stderr_pix = []
    high_SN_m1_psf = []
    high_SN_m2_psf = []
    high_SN_c1_psf = []
    high_SN_c2_psf = []
    high_SN_m1_stderr_psf = []
    high_SN_m2_stderr_psf = []
    high_SN_c1_stderr_psf = []
    high_SN_c2_stderr_psf = []
    
    for team_data_file in low_SN_data_files:
        data = fits.open(team_data_file)[2].data
        
        low_SN_m1_pix.append(data.field("m")[0])
        low_SN_m2_pix.append(data.field("m")[1])
        low_SN_m1_psf.append(data.field("m")[2])
        low_SN_m2_psf.append(data.field("m")[3])
        low_SN_c1_pix.append(data.field("c")[0])
        low_SN_c2_pix.append(data.field("c")[1])
        low_SN_c1_psf.append(data.field("c")[2])
        low_SN_c2_psf.append(data.field("c")[3])
        low_SN_m1_stderr_pix.append(data.field("m_stderr")[0])
        low_SN_m2_stderr_pix.append(data.field("m_stderr")[1])
        low_SN_m1_stderr_psf.append(data.field("m_stderr")[2])
        low_SN_m2_stderr_psf.append(data.field("m_stderr")[3])
        low_SN_c1_stderr_pix.append(data.field("c_stderr")[0])
        low_SN_c2_stderr_pix.append(data.field("c_stderr")[1])
        low_SN_c1_stderr_psf.append(data.field("c_stderr")[2])
        low_SN_c2_stderr_psf.append(data.field("c_stderr")[3])
    
    for team_data_file in high_SN_data_files:
        data = fits.open(team_data_file)[2].data
        
        high_SN_m1_pix.append(data.field("m")[0])
        high_SN_m2_pix.append(data.field("m")[1])
        high_SN_m1_psf.append(data.field("m")[2])
        high_SN_m2_psf.append(data.field("m")[3])
        high_SN_c1_pix.append(data.field("c")[0])
        high_SN_c2_pix.append(data.field("c")[1])
        high_SN_c1_psf.append(data.field("c")[2])
        high_SN_c2_psf.append(data.field("c")[3])
        high_SN_m1_stderr_pix.append(data.field("m_stderr")[0])
        high_SN_m2_stderr_pix.append(data.field("m_stderr")[1])
        high_SN_m1_stderr_psf.append(data.field("m_stderr")[2])
        high_SN_m2_stderr_psf.append(data.field("m_stderr")[3])
        high_SN_c1_stderr_pix.append(data.field("c_stderr")[0])
        high_SN_c2_stderr_pix.append(data.field("c_stderr")[1])
        high_SN_c1_stderr_psf.append(data.field("c_stderr")[2])
        high_SN_c2_stderr_psf.append(data.field("c_stderr")[3])
    
    # Get the mean abs for each value
    low_SN_m_pix = 0.5*(np.abs(low_SN_m1_pix)+np.abs(low_SN_m2_pix))
    low_SN_m_psf = 0.5*(np.abs(low_SN_m1_psf)+np.abs(low_SN_m2_psf))
    low_SN_c_pix = 0.5*(np.abs(low_SN_c1_pix)+np.abs(low_SN_c2_pix))
    low_SN_c_psf = 0.5*(np.abs(low_SN_c1_psf)+np.abs(low_SN_c2_psf))
    
    high_SN_m_pix = 0.5*(np.abs(high_SN_m1_pix)+np.abs(high_SN_m2_pix))
    high_SN_m_psf = 0.5*(np.abs(high_SN_m1_psf)+np.abs(high_SN_m2_psf))
    high_SN_c_pix = 0.5*(np.abs(high_SN_c1_pix)+np.abs(high_SN_c2_pix))
    high_SN_c_psf = 0.5*(np.abs(high_SN_c1_psf)+np.abs(high_SN_c2_psf))
    
    # Get the errors of the mean abses
    low_SN_m_stderr_pix = 0.5*np.sqrt(np.square(low_SN_m1_stderr_pix)+np.square(low_SN_m2_stderr_pix))
    low_SN_m_stderr_psf = 0.5*np.sqrt(np.square(low_SN_m1_stderr_psf)+np.square(low_SN_m2_stderr_psf))
    low_SN_c_stderr_pix = 0.5*np.sqrt(np.square(low_SN_c1_stderr_pix)+np.square(low_SN_c2_stderr_pix))
    low_SN_c_stderr_psf = 0.5*np.sqrt(np.square(low_SN_c1_stderr_psf)+np.square(low_SN_c2_stderr_psf))
    
    high_SN_m_stderr_pix = 0.5*np.sqrt(np.square(high_SN_m1_stderr_pix)+np.square(high_SN_m2_stderr_pix))
    high_SN_m_stderr_psf = 0.5*np.sqrt(np.square(high_SN_m1_stderr_psf)+np.square(high_SN_m2_stderr_psf))
    high_SN_c_stderr_pix = 0.5*np.sqrt(np.square(high_SN_c1_stderr_pix)+np.square(high_SN_c2_stderr_pix))
    high_SN_c_stderr_psf = 0.5*np.sqrt(np.square(high_SN_c1_stderr_psf)+np.square(high_SN_c2_stderr_psf))
    
    # Set up indices for x-axis
    team_indices = np.array(range(num_teams))
    
    # Plot values
    for (low_SN_y_array, low_SN_y_stderr_array, high_SN_y_array, high_SN_y_stderr_array, y_limits, filename, y_label) in \
                    ((low_SN_m_pix,low_SN_m_stderr_pix,high_SN_m_pix,high_SN_m_stderr_pix,m_limits,"SBE_m_pix.eps",r"$\overline{|m|}$ in pixel frame"),
                     (low_SN_m_psf,low_SN_m_stderr_psf,high_SN_m_psf,high_SN_m_stderr_psf,m_limits,"SBE_m_psf.eps",r"$\overline{|m|}$ in psf frame"),
                     (low_SN_c_pix,low_SN_c_stderr_pix,high_SN_c_pix,high_SN_c_stderr_pix,c_limits,"SBE_c_pix.eps",r"$\overline{|c|}$ in pixel frame"),
                     (low_SN_c_psf,low_SN_c_stderr_psf,high_SN_c_psf,high_SN_c_stderr_psf,c_limits,"SBE_c_psf.eps",r"$\overline{|c|}$ in psf frame")):
        
        # Setup the figure
        fig = pyplot.figure()
        fig.subplots_adjust(wspace=0, hspace=0, bottom=0.1, right=0.95, top=0.95, left=0.12)
        
        ax = fig.add_subplot(1,1,1)
        ax.set_xlabel("Team",fontsize=fontsize)
        ax.set_ylabel(y_label,fontsize=fontsize)
        
        # Plot the points
        ax.errorbar(team_indices-SN_shift,low_SN_y_array,color='r',marker='d',linestyle='None',
                    markersize=10,label="Low S/N",yerr=low_SN_y_stderr_array)
        ax.errorbar(team_indices+SN_shift,high_SN_y_array,color='b',marker='o',linestyle='None',
                    markersize=10,label="High S/N",yerr=high_SN_y_stderr_array)
        
        # Draw the target line
        if(y_limits==m_limits):
            target = m_target
        else:
            target = c_target
        ax.plot([-1,num_teams],[target,target],label=None,color="k",linestyle="dotted")
            
        # Draw the legend
        ax.legend(numpoints=1,loc='lower left')
        
        # Set up the axes
        ax.set_xlim(-1,num_teams)
        ax.set_xticks(range(num_teams))
        ax.set_xticklabels(team_names,fontsize=fontsize)
        
        ax.set_yscale("log", nonposy='clip')
        ax.set_ylim(y_limits)
        ax.tick_params(axis='y', which='major', labelsize=fontsize)
        
        
        # Save the figure
        pyplot.savefig(filename, format="eps", bbox_inches="tight", pad_inches=0.05)
        
        print("Finished plotting to " + filename + ".")
        

    
    # Plot errors
    for (low_SN_y_stderr_array, high_SN_y_stderr_array, y_limits, filename, y_label) in \
                    ((low_SN_m_stderr_pix,high_SN_m_stderr_pix,m_limits,"SBE_m_err_pix.eps",r"$\sigma(|m|)$ in pixel frame"),
                     (low_SN_m_stderr_psf,high_SN_m_stderr_psf,m_limits,"SBE_m_err_psf.eps",r"$\sigma(|m|)$ in psf frame"),
                     (low_SN_c_stderr_pix,high_SN_c_stderr_pix,c_limits,"SBE_c_err_pix.eps",r"$\sigma(|c|)$ in pixel frame"),
                     (low_SN_c_stderr_psf,high_SN_c_stderr_psf,c_limits,"SBE_c_err_psf.eps",r"$\sigma(|c|)$ in psf frame")):
        
        # Setup the figure
        fig = pyplot.figure()
        fig.subplots_adjust(wspace=0, hspace=0, bottom=0.1, right=0.95, top=0.95, left=0.12)
        
        ax = fig.add_subplot(1,1,1)
        ax.set_xlabel("Team",fontsize=fontsize)
        ax.set_ylabel(y_label,fontsize=fontsize)
        
        # Plot the points
        ax.errorbar(team_indices-SN_shift,low_SN_y_stderr_array,color='r',marker='d',linestyle='None',
                    markersize=10,label="Low S/N")
        ax.errorbar(team_indices+SN_shift,high_SN_y_stderr_array,color='b',marker='o',linestyle='None',
                    markersize=10,label="High S/N")
        
        # Draw the target line
        if(y_limits==m_limits):
            target = m_target
        else:
            target = c_target
        ax.plot([-1,num_teams],[target,target],label=None,color="k",linestyle="dotted")
            
        # Draw the legend
        ax.legend(numpoints=1,loc='lower left')
        
        # Set up the axes
        ax.set_xlim(-1,num_teams)
        ax.set_xticks(range(num_teams))
        ax.set_xticklabels(team_names,fontsize=fontsize)
        
        ax.set_yscale("log", nonposy='clip')
        ax.set_ylim(y_limits)
        ax.tick_params(axis='y', which='major', labelsize=fontsize)
        
        
        # Save the figure
        pyplot.savefig(filename, format="eps", bbox_inches="tight", pad_inches=0.05)
        
        print("Finished plotting to " + filename + ".")
    

if __name__ == "__main__":
    main(sys.argv)