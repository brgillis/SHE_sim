'''
Created on 2 Jun 2015

@author: brg
'''

import numpy as np
import click

from astropy.io import fits

import matplotlib
import matplotlib.pyplot as pyplot
matplotlib.rcParams['ps.useafm'] = True
matplotlib.rcParams['pdf.use14corefonts'] = True
matplotlib.rcParams['text.usetex'] = True

@click.command()
@click.option("--mode",default="SBE",help="'SBE' for SBE plotting, 'calibration' for calibration")
def main(**kwargs):
    
    if(kwargs['mode']=="SBE"):
    
        team_names = ("Control (KSB)", "1", "2", "3", "4")
        
        low_SN_data_files = ("/disk2/brg/Program_Files/workspace/Generate_GalSim_Images/benchmark_low_SN_v3/bias_measurements_GalSim_KSB.fits",
                              "/home/brg/Data/SBE/MegaLUT/low_SN/stored_data.fits",
                              "/home/brg/Data/SBE/Oxford/low_SN/stored_data.fits",
                              "/home/brg/Data/SBE/gfit/low_SN/stored_data.fits",
                              "/home/brg/Data/SBE/FDNT/low_SN/stored_data.fits")
        high_SN_data_files = ("/disk2/brg/Program_Files/workspace/Generate_GalSim_Images/benchmark_high_SN_v3/bias_measurements_GalSim_KSB.fits",
                              "/home/brg/Data/SBE/MegaLUT/high_SN/stored_data.fits",
                              "/home/brg/Data/SBE/Oxford/high_SN/stored_data.fits",
                              "/home/brg/Data/SBE/gfit/high_SN/stored_data.fits",
                              "/home/brg/Data/SBE/FDNT/high_SN/stored_data.fits")
        m_limits = (0.001,0.3)
        c_limits = (0.5e-5,1.5e-3)
        m_err_limits = (0,0.01)
        c_err_limits = (0,1e-4)
        
        x_label = "Team"
        
    elif(kwargs['mode']=="calibration"):
    
        team_names = ("Uncalibrated", "Calibrated")
        
        low_SN_data_files = ("/disk2/brg/Program_Files/workspace/Generate_GalSim_Images/benchmark_low_SN_v3/bias_measurements_GalSim_KSB.fits",
                             "/disk2/brg/Program_Files/workspace/Generate_GalSim_Images/benchmark_low_SN_v3/bias_measurements_GalSim_KSB_calibrated.fits")
        high_SN_data_files = ("/disk2/brg/Program_Files/workspace/Generate_GalSim_Images/benchmark_high_SN_v3/bias_measurements_GalSim_KSB.fits",
                              "/disk2/brg/Program_Files/workspace/Generate_GalSim_Images/benchmark_high_SN_v3/bias_measurements_GalSim_KSB_calibrated.fits")
        m_limits = (0.00001,0.3)
        c_limits = (0.5e-7,1.5e-3)
        m_err_limits = (0,0.012)
        c_err_limits = (0,1.2e-4)
        
        x_label = "Type"
    
    num_teams = len(team_names)
    
    m_target = 1e-4
    c_target = 5e-4
    
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
        data = fits.open(team_data_file)[-1].data #FIXME - temporary kludge till formats are normalized
        
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
        
        label = "SBE"
        
        y_m_label_head = r"$\sigma(|m|)/(1+m)$"
        y_c_label_head = r"$\sigma(|c|)/(1+m)$"
    
    for team_data_file in high_SN_data_files:
        data = fits.open(team_data_file)[-1].data #FIXME - temporary kludge till formats are normalized
        
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
        
        label = "calibration"
        
        y_m_label_head = r"$\sigma(|m|)$"
        y_c_label_head = r"$\sigma(|c|)$"
    
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
    
    # Get calibration errors
    if(kwargs['mode'] != "calibration"):
        low_SN_m_cal_err_pix = 0.5*np.sqrt(np.square(low_SN_m1_stderr_pix/np.add(1,low_SN_m1_pix))+
                                           np.square(low_SN_m2_stderr_pix/np.add(1,low_SN_m2_pix)))
        low_SN_m_cal_err_psf = 0.5*np.sqrt(np.square(low_SN_m1_stderr_psf/np.add(1,low_SN_m1_psf))+
                                           np.square(low_SN_m2_stderr_psf/np.add(1,low_SN_m2_psf)))
        low_SN_c_cal_err_pix = 0.5*np.sqrt(np.square(low_SN_c1_stderr_pix/np.add(1,low_SN_m1_pix))+
                                           np.square(low_SN_c2_stderr_pix/np.add(1,low_SN_m2_pix)))
        low_SN_c_cal_err_psf = 0.5*np.sqrt(np.square(low_SN_c1_stderr_psf/np.add(1,low_SN_m1_psf))+
                                           np.square(low_SN_c2_stderr_psf/np.add(1,low_SN_m2_psf)))
        
        high_SN_m_cal_err_pix = 0.5*np.sqrt(np.square(high_SN_m1_stderr_pix/np.add(1,high_SN_m1_pix))+
                                            np.square(high_SN_m2_stderr_pix/np.add(1,high_SN_m2_pix)))
        high_SN_m_cal_err_psf = 0.5*np.sqrt(np.square(high_SN_m1_stderr_psf/np.add(1,high_SN_m1_psf))+
                                            np.square(high_SN_m2_stderr_psf/np.add(1,high_SN_m2_psf)))
        high_SN_c_cal_err_pix = 0.5*np.sqrt(np.square(high_SN_c1_stderr_pix/np.add(1,high_SN_m1_pix))+
                                            np.square(high_SN_c2_stderr_pix/np.add(1,high_SN_m2_pix)))
        high_SN_c_cal_err_psf = 0.5*np.sqrt(np.square(high_SN_c1_stderr_psf/np.add(1,high_SN_m1_psf))+
                                            np.square(high_SN_c2_stderr_psf/np.add(1,high_SN_m2_psf)))
    else:
        low_SN_m_cal_err_pix = low_SN_m_stderr_pix
        low_SN_m_cal_err_psf = low_SN_m_stderr_psf
        low_SN_c_cal_err_pix = low_SN_c_stderr_pix
        low_SN_c_cal_err_psf = low_SN_c_stderr_psf
        
        high_SN_m_cal_err_pix = high_SN_m_stderr_pix
        high_SN_m_cal_err_psf = high_SN_m_stderr_psf
        high_SN_c_cal_err_pix = high_SN_c_stderr_pix
        high_SN_c_cal_err_psf = high_SN_c_stderr_psf
    
    # Set up indices for x-axis
    team_indices = np.array(range(num_teams))
    
    # Plot values
    for (low_SN_y_array, low_SN_y_stderr_array, high_SN_y_array, high_SN_y_stderr_array, y_limits, filename, y_label) in \
                    ((low_SN_m_pix,low_SN_m_stderr_pix,high_SN_m_pix,high_SN_m_stderr_pix,m_limits,label + "_m_pix.eps",r"$\overline{|m|}$ in pixel frame"),
                     (low_SN_m_psf,low_SN_m_stderr_psf,high_SN_m_psf,high_SN_m_stderr_psf,m_limits,label + "_m_psf.eps",r"$\overline{|m|}$ in psf frame"),
                     (low_SN_c_pix,low_SN_c_stderr_pix,high_SN_c_pix,high_SN_c_stderr_pix,c_limits,label + "_c_pix.eps",r"$\overline{|c|}$ in pixel frame"),
                     (low_SN_c_psf,low_SN_c_stderr_psf,high_SN_c_psf,high_SN_c_stderr_psf,c_limits,label + "_c_psf.eps",r"$\overline{|c|}$ in psf frame")):
        
        # Setup the figure
        fig = pyplot.figure()
        fig.subplots_adjust(wspace=0, hspace=0, bottom=0.1, right=0.95, top=0.95, left=0.12)
        
        ax = fig.add_subplot(1,1,1)
        ax.set_xlabel(x_label,fontsize=fontsize)
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
                    ((low_SN_m_cal_err_pix,high_SN_m_cal_err_pix,m_err_limits,label + "_m_err_pix.eps",y_m_label_head + " in pixel frame"),
                     (low_SN_m_cal_err_psf,high_SN_m_cal_err_psf,m_err_limits,label + "_m_err_psf.eps",y_m_label_head + " in psf frame"),
                     (low_SN_c_cal_err_pix,high_SN_c_cal_err_pix,c_err_limits,label + "_c_err_pix.eps",y_c_label_head + " in pixel frame"),
                     (low_SN_c_cal_err_psf,high_SN_c_cal_err_psf,c_err_limits,label + "_c_err_psf.eps",y_c_label_head + " in psf frame")):
        
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
        if kwargs['mode']!="calibration":
            if(y_limits==m_err_limits):
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
        
        #ax.set_yscale("log", nonposy='clip')
        ax.set_ylim(y_limits)
        ax.tick_params(axis='y', which='major', labelsize=fontsize)
        
        
        # Save the figure
        pyplot.savefig(filename, format="eps", bbox_inches="tight", pad_inches=0.05)
        
        print("Finished plotting to " + filename + ".")
        
    return
    

if __name__ == "__main__":
    main()