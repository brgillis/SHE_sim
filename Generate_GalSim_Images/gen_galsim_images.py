"""gen_galsim_images.py
   Created by Bryan Gillis, March 2014
   Last edited by brg, 21 April 2014
   
   Requirements: GalSim toolkit (and its requirements).
   GalSim can be downloaded from https://github.com/GalSim-developers/GalSim, and see
   https://github.com/GalSim-developers/GalSim/wiki/Installation-Instructions for its
   installation requirements.

   This is a script to generate a set of seeded random images using the GalSim toolkit,
   suitable for Euclid benchmarking tests of shear-measurement algorithms.
   This script should be invoked with a command that points it to an appropriate
   configuration script, for instance through:
   
   $ python gen_galsim_images.py my_config.cfg

   See the sample configuration files included with this script for the proper format.
   It is recommended that you start with an appropriate configuration script and modify
   it to suit your purposes.
   
   If a configuration file is not passed to this script, it will use the set of default
   configuration values assigned in the load_default_configurations(...) function below.
   
   NOTE: If you got this script or the configuration file from someone else, check the
   output directory and change it if necessary. It can use absolute paths, and the
   folders will be created if necessary, which may result in creating a file structure
   you don't want.
"""

import sys
import math
import numpy           as np
import galsim
import multiprocessing as mtp
import subprocess
    
# Default values
default_random_seed = 8241573
default_shear_type =   'shear'
default_image_type =   '32f'
default_gain = 1.7

def load_default_configurations(config_dict):
    """This function loads a default set of configuration parameters. If you wish to run
       this script without a configuration file, you can edit the parameters here. If you
       do so, ensure that all lines are entered in lower-case, which is what the program
       will be expecting.
    """
    
    print "No configuration script loaded. The script will proceed using the set of"
    print "configuration parameters coded into it. These can be viewed and edited in"
    print "the 'load_default_configurations' function within the script."
    
    # Bookkeeping parameters
    config_dict['config_mode']             =   "distribution"
    config_dict['output_folder_name']      =   "output"
    config_dict['thread_folder_name_base'] =   "thread"
    config_dict['output_name_base']        =   "galaxy_image"
    config_dict['num_threads']             =   1
    config_dict['num_files']               =   1
    config_dict['grid_or_cube']            =   "grid"
    if(config_dict['grid_or_cube'] == 'cube'):
        config_dict['num_per_file']        =   10000
    elif(config_dict['grid_or_cube'] == 'grid'):
        config_dict['num_per_row']         =   10
        config_dict['num_per_col']         =   10
        config_dict['num_per_file']        =   config_dict['num_per_row'] * config_dict['num_per_col']
    else:
        return bad_config_format()
    config_dict['pixel_scale']             =   0.1
    config_dict['sample_scale']            =   0.05
    config_dict['image_size_x_pix']        =   200
    config_dict['image_size_y_pix']        =   200
    config_dict['stamp_padding_x_pix']     =   1
    config_dict['stamp_padding_y_pix']     =   1
    config_dict['image_type']              =   "32f"
    config_dict['shear_type']              =   "shear-angle"
    config_dict['use_flux_or_s/n']         =   "flux"
    config_dict['init_random_seed_factor'] =   0
    config_dict['suppress_noise']          =   0
    config_dict['shape_noise_cancellation']=   1
    
    # Check if other parameters are consistent with shape noise cancellation if it's enabled
    if((config_dict['shape_noise_cancellation']!=0) and
       (config_dict['shear_type']!='shear-angle') and 
       (config_dict['shear_type']!='ellipticity-angle') ):
        if(config_dict['shear_type']!='ellipticity'):
            config_dict['shear_type'] = 'ellipticity-angle'
        else:
            config_dict['shear_type'] = 'shear-angle'
        print("ERROR: Shear type is incompatible with shape noise cancellation.")
        print("It has been changed to " +  config_dict['shear_type'] + ". Be warned that this")
        print("will likely result in improper values for shapes and shears. Either disable")
        print("shape noise cancellation or correct it to one of the '-angle' types.")
        
    if((config_dict['shape_noise_cancellation']!=0) and
       (config_dict['num_per_file'] % 2 == 1)):
        print("WARNING: Shape noise cancellation requires an even number of images per")
        print("file. The number of images per file has been doubled to compensate.")
        if(config_dict['grid_or_cube'] == 'cube'):
            config_dict['num_per_file'] *= 2
        elif(config_dict['grid_or_cube'] == 'grid'):
            config_dict['num_per_row'] *= 2
            config_dict['num_per_file'] *= 2
        
        
    # Check whether we're using the simple or complex config_dict['config_mode']
    if(config_dict['config_mode'] == 'simple'):
        
        # Using simple mode, so we'll just use fixed values for all parameters.
        # Only the noise realisation will differ between images.
        
        config_dict['psf_stddev_arcsec']      = 0.2
        config_dict['psf_shape_1']            = 0
        config_dict['psf_shape_2']            = 0
        config_dict['galaxy_stddev_arcsec']   = 0.3
        config_dict['galaxy_flux_or_SN']      = 40000
        config_dict['galaxy_shape_1']         = 0.2
        config_dict['galaxy_shape_2']         = 120
        config_dict['galaxy_shear_1']         = 0.03
        config_dict['galaxy_shear_2']         = 45
        config_dict['sky_level_subtracted']   = 100
        config_dict['sky_level_unsubtracted'] = 0
        config_dict['read_noise']             = 0.3
        config_dict['gain']                   = 1.7
        
    elif(config_dict['config_mode'] == 'distribution'):
        
        # Using distribution mode, so we'll read in using the get_dist function,
        # which determines the distribution type from the second line and the
        # parameters for it from the third. The first line tells whether it'll
        # be randomly generated once per file or once per postage stamp.
        
        config_dict['psf_stddev_arcsec_dist'] = \
            get_dist("LogNormal-Peak", "-0.7 0.005", "File")
        config_dict['psf_shape_1_dist'] = \
            get_dist("Contracted_Rayleigh", "0.01 0.9 4", "File")
        config_dict['psf_shape_2_dist'] = \
            get_dist("Uniform", "0 180", "File")
        config_dict['galaxy_stddev_arcsec_dist'] = \
            get_dist("LogNormal-Peak", "-0.5 0.15", "Stamp")
        config_dict['galaxy_SN_dist'] = \
            get_dist("Fixed", "40000", "File")
        config_dict['galaxy_shape_1_dist'] = \
            get_dist("Contracted_Rayleigh", "0.25 0.9 4", "Stamp")
        config_dict['galaxy_shape_2_dist'] = \
            get_dist("Uniform", "0 180", "Stamp")
        config_dict['galaxy_shear_1_dist'] = \
            get_dist("Contracted_Rayleigh", "0.01 0.9 4", "File")
        config_dict['galaxy_shear_2_dist'] = \
            get_dist("Uniform", "0 180", "File")
        config_dict['sky_level_subtracted_dist'] = \
            get_dist("Fixed", "100", "File")
        config_dict['sky_level_unsubtracted_dist'] = \
            get_dist("Fixed", "0", "File")
        config_dict['read_noise_dist'] = \
            get_dist("Fixed", "0.3", "File")
        config_dict['gain_dist'] = \
            get_dist("Fixed", "1.7", "File")
            
        if(config_dict['shape_noise_cancellation']!=0):
            # Enforce that galaxy and psf shape are generated once per file
            config_dict['psf_shape_1_dist'] = (config_dict['psf_shape_1_dist'][0],
                                                  config_dict['psf_shape_1_dist'][1],
                                                  'file')
            config_dict['psf_shape_2_dist'] = (config_dict['psf_shape_2_dist'][0],
                                                  config_dict['psf_shape_2_dist'][1],
                                                  'file')
            config_dict['galaxy_shape_1_dist'] = (config_dict['galaxy_shape_1_dist'][0],
                                                  config_dict['galaxy_shape_1_dist'][1],
                                                  'file')
            config_dict['galaxy_shape_2_dist'] = (config_dict['galaxy_shape_2_dist'][0],
                                                  config_dict['galaxy_shape_2_dist'][1],
                                                  'file')

def main(argv):
    
    # Magic numbers
    
    # The smallest possible configuration file length
    config_min_length              = 15
    
    # Lists of allowed values for certain configuration parameters
    allowed_shear_types = ['ellipticity', 'shear', 'ellipticity-angle', 'shear-angle']
    allowed_image_types = ['16i', '32i', '32f', '64f']
    
    # Hard bounds for input
    min_num_files    = 1
    min_num_per_file = 1
    min_image_size   = 2
    min_pixel_scale  = 0.0001
        
    # Initialise the configuration dictionary
    config_dict = {};
    
    # Check that a configuration file name was passed at command line
    if(len(argv)) <= 1:
        load_default_configurations(config_dict)
        
    else:
    
        # Open and read in the config file
        config_file_name = argv[1]
        with open(config_file_name, 'r') as config_file:
        
            # Read in the file, except for comment lines
            config_lines = []
            for config_line in config_file:
                if((config_line[0] != '#') and (len(config_line.strip()) > 0)):
                    config_lines.append(config_line.strip())
                    
            # The reading here is pretty primative: It assumes a specific line order, and takes values
            # from the end of each line. There's no check for in-line comments or that each line
            # actually has the proper parameter in it.
                    
            # Check it has at least the minimum length. Otherwise we know something's wrong
            if(len(config_lines) < config_min_length):
                return bad_config_format()
            
            # Get configuration file version from first line
            version = str(config_lines[0].split()[-1])
            
            # Check which configuration file version it is. This allows backwards-compatibility with
            # older config files.
            if(version == '1.0'):
                
                if(load_config_1_0(config_dict, config_lines)): return
                    
            elif(version == '1.1'):
                
                if(load_config_1_1(config_dict, config_lines)): return
                
            elif(version == '1.2'):
                
                if(load_config_1_2(config_dict, config_lines)): return
                
            elif(version == '1.3'): 
                
                if(load_config_1_3(config_dict, config_lines)): return
                
            else:
                print "ERROR: Unsupported config file version or improperly formatted"
                print "configuration file. Please check sample files and ensure it has"
                print "one of the proper formats and versions."
                return
            
            # Check the values we read in against hard bounds and adjust if necessary,
            # printing a warning.
            
            if(config_dict['num_files'] < min_num_files):
                config_dict['num_files'] = min_num_files
                print "WARNING: Adjusted num_files to minimum of " + str(min_num_files)
                
            if(config_dict['num_per_file'] < min_num_per_file):
                config_dict['num_per_file'] = min_num_per_file
                if(config_dict['grid_or_cube'] != 'cube'):
                    config_dict['num_per_row'] = 1
                    config_dict['num_per_col'] = 1
                print "WARNING: Adjusted num_per_file to minimum of " + str(min_num_per_file)
                
            if(config_dict['pixel_scale'] < min_pixel_scale):
                config_dict['pixel_scale'] = min_pixel_scale
                print "WARNING: Adjusted pixel_scale to minimum of " + str(min_pixel_scale)
                print "Check you're using units of arcsec/pixel!"
                
            if(config_dict['sample_scale'] < min_pixel_scale):
                config_dict['sample_scale'] = min_pixel_scale
                print "WARNING: Adjusted sample_scale to minimum of " + str(min_pixel_scale)
                print "Check you're using units of arcsec/pixel!"
                
            if(config_dict['image_size_x_pix'] < min_image_size):
                config_dict['image_size_x_pix'] = min_image_size
                print "WARNING: Adjusted config_dict['image_size_x_pix'] to minimum of " + str(min_image_size)
                
            if(config_dict['image_size_y_pix'] < min_image_size):
                config_dict['image_size_y_pix'] = min_image_size
                print "WARNING: Adjusted config_dict['image_size_y_pix'] to minimum of " + str(min_image_size)
                
            shear_type_found = False
            for im_type in allowed_shear_types:
                if(config_dict['shear_type'] == im_type):
                    shear_type_found = True
                    break
            if not shear_type_found:
                config_dict['shear_type'] = default_shear_type
                print "WARNING: Shear im_type not known. It's been adjusted to"
                print "the default im_type (" + default_shear_type + ")."
                
            image_type_found = False
            for im_type in allowed_image_types:
                if(config_dict['image_type'] == im_type):
                    image_type_found = True
                    break
            if not image_type_found:
                config_dict['image_type'] = default_image_type
                print "WARNING: Image im_type not known. It's been adjusted to"
                print "the default im_type (" + default_image_type + ")."
            
            # End reading in the configuration file
        
    # Ensure the base output folder exists
    # At present there's no checking before trying to create it, so an error will
    # be printed if it already exists. This can be disregarded.
    
    # Check if the folder path was given with a slash at the end. If so, trim it
    if(config_dict['output_folder_name'][-1] == '/'):
        config_dict['output_folder_name'] = config_dict['output_folder_name'][0:-1]
        
    # Ensure that the output folder exists
    cmd = 'mkdir -p ' + config_dict['output_folder_name']
    subprocess.call(cmd, shell=True)
            
    # And now we allow the script to split into separate threads for image creation
    processes = []
    if(config_dict['num_threads'] <= 1):
        # For debugging purposes, we separate this to a simple function call
        generator(config_dict, 0)
    else:
        # Call each thread
        for thread_num in xrange(config_dict['num_threads']):
            p = mtp.Process(target=generator, args=( config_dict, thread_num ) )
            p.start()
            processes.append(p)
            
        # Join each thread back up so they can finish.
        for p in processes:
            p.join()
                
    print 'Execution complete.'
            
            
def generator(config_dict, thread_num=0):
    """This function does the heavy lifting of actually generating images with the GalSim toolkit.
       It's called by main and passed a dictionary with configuration parameters to use, and told
       which thread number it's being called from. It uses the thread number to set its rng seed
       and to determine where to save its files.
    """
    
    # Magic numbers for hard parameter bounds
    
    shear_max      = 0.99
    stddev_min     = 1e-9
    SN_min         = 1e-9
    sky_level_min  = 0
    read_noise_min = 0
    gain_min       = 0.00001
    
    # Set up this thread
    
    # Set up the initial random seed for this thread
    seed_factor = config_dict['num_files']*config_dict['num_per_file'] # So it's large enough for no overlapping seeds
    init_random_seed = config_dict['init_random_seed_factor'] * seed_factor * config_dict['num_threads'] + \
                                                                seed_factor * thread_num
        
    # Ensure that the folder for this thread exists
    thread_folder_name = config_dict['output_folder_name']      + \
                         '/'                                    + \
                         config_dict['thread_folder_name_base'] + \
                         '_'                                    + \
                         str(thread_num) 
    cmd = 'mkdir -p ' + thread_folder_name
    subprocess.call(cmd, shell=True)
    
    file_name_base = thread_folder_name + '/' + config_dict['output_name_base'] + '_'
    
    # Loop through each of the files we want to create, generating images for them
    for i in xrange(0, config_dict['num_files']):

        # Setup for the file
        if(config_dict['grid_or_cube'] == 'cube'):
            # We'll keep a stack of individual stamps
            images = []
            
        else:
            # We'll generate one big image, with stamps on a grid
            full_x_size = config_dict['image_size_x_pix']    *  config_dict['num_per_row'] + \
                          config_dict['stamp_padding_x_pix'] * (config_dict['num_per_row'] - 1)
            full_y_size = config_dict['image_size_y_pix']    *  config_dict['num_per_col'] + \
                          config_dict['stamp_padding_y_pix'] * (config_dict['num_per_col'] - 1)
            
            # Create the image object, using the appropriate method for the image type
            if(config_dict['image_type'] == '16i'):
                full_image = galsim.ImageC(full_x_size , full_y_size, scale=config_dict['sample_scale'])
            elif(config_dict['image_type'] == '32i'):
                full_image = galsim.ImageI(full_x_size , full_y_size, scale=config_dict['sample_scale'])
            elif(config_dict['image_type'] == '32f'):
                full_image = galsim.ImageF(full_x_size , full_y_size, scale=config_dict['sample_scale'])
            elif(config_dict['image_type'] == '64f'):
                full_image = galsim.ImageD(full_x_size , full_y_size, scale=config_dict['sample_scale'])
            else:
                print "ERROR: Bad image type slipped through somehow."
                return
            
            # Set up ix and iy as running counters for where in the grid we are
            ix = 0
            iy = 0
            
        # Initialise our array of output data strings, which we'll later print out
        data_strings = []
        
        # Seed numpy.random and the galsim deviate for the parameters that are generated per-file
        # Note that this is the same seed as is used for the first in each file, but since
        # no parameter will be generated in both places, this won't cause issues
        np.random.seed(init_random_seed + i * config_dict['num_per_file'])
        rng = galsim.BaseDeviate(init_random_seed + i * config_dict['num_per_file'])
        
        # For the parameters that are randomly generated once per file, generate them now
        if(config_dict['config_mode'] == 'distribution'):
            if(config_dict['psf_stddev_arcsec_dist'][2]      == 'file'):
                config_dict['psf_stddev_arcsec']                 = get_rand_from_dist(config_dict['psf_stddev_arcsec_dist'], rng)
            if(config_dict['psf_shape_1_dist'][2]            == 'file'):
                config_dict['psf_shape_1']                       = get_rand_from_dist(config_dict['psf_shape_1_dist'], rng)
            if(config_dict['psf_shape_2_dist'][2]            == 'file'):
                config_dict['psf_shape_2']                       = get_rand_from_dist(config_dict['psf_shape_2_dist'], rng)
            if(config_dict['galaxy_stddev_arcsec_dist'][2]   == 'file'):
                config_dict['galaxy_stddev_arcsec']              = get_rand_from_dist(config_dict['galaxy_stddev_arcsec_dist'], rng)
            if(config_dict['galaxy_SN_dist'][2]              == 'file'):
                config_dict['galaxy_flux_or_SN']                         = get_rand_from_dist(config_dict['galaxy_SN_dist'], rng)
            if(config_dict['galaxy_shape_1_dist'][2]         == 'file'):
                config_dict['galaxy_shape_1']                    = get_rand_from_dist(config_dict['galaxy_shape_1_dist'], rng)
            if(config_dict['galaxy_shape_2_dist'][2]         == 'file'):
                config_dict['galaxy_shape_2']                    = get_rand_from_dist(config_dict['galaxy_shape_2_dist'], rng)
            if(config_dict['galaxy_shear_1_dist'][2]         == 'file'):
                config_dict['galaxy_shear_1']                    = get_rand_from_dist(config_dict['galaxy_shear_1_dist'], rng)
            if(config_dict['galaxy_shear_2_dist'][2]         == 'file'):
                config_dict['galaxy_shear_2']                    = get_rand_from_dist(config_dict['galaxy_shear_2_dist'], rng)
            if(config_dict['sky_level_subtracted_dist'][2]   == 'file'):
                config_dict['sky_level_subtracted']              = get_rand_from_dist(config_dict['sky_level_subtracted_dist'], rng)
            if(config_dict['sky_level_unsubtracted_dist'][2] == 'file'):
                config_dict['sky_level_unsubtracted']            = get_rand_from_dist(config_dict['sky_level_unsubtracted_dist'], rng)
            if(config_dict['read_noise_dist'][2]             == 'file'):
                config_dict['read_noise']                        = get_rand_from_dist(config_dict['read_noise_dist'], rng)
            if(config_dict['gain_dist'][2]                   == 'file'):
                config_dict['gain']                              = get_rand_from_dist(config_dict['gain_dist'], rng)
                
        # If shape noise cancellation is being applied, we'll use the angle for galaxy shape we generated
        # as the initial angle. Also, set up parameters for the cancellation
        if(config_dict['shape_noise_cancellation']!=0):
            angle_step = 90. / (config_dict['num_per_file']/2.)
            config_dict['galaxy_shape_2'] += 90 - angle_step
            first_in_pair = True

        # Loop through each of the stamps we want to create within this file
        for j in xrange(0, config_dict['num_per_file']):
            
            # Seed numpy.random and the galsim deviate for the parameters that are generated per-stamp
            np.random.seed(init_random_seed + i * config_dict['num_per_file'] + j)
            rng = galsim.BaseDeviate(init_random_seed + i * config_dict['num_per_file'] + j)
            
            # Handling for shape noise cancellation
            if(config_dict['shape_noise_cancellation']!=0):
                if(first_in_pair):
                    first_in_pair = False
                    
                    # Increment angle
                    config_dict['galaxy_shape_2'] -= 90
                    config_dict['galaxy_shape_2'] += angle_step
                    
                    # Generate values for this pair as appropriate
                    if(config_dict['config_mode'] == 'distribution'):
                        if(config_dict['psf_stddev_arcsec_dist'][2]      != 'file'):
                            config_dict['psf_stddev_arcsec']                 = get_rand_from_dist(config_dict['psf_stddev_arcsec_dist'], rng)
                        if(config_dict['psf_shape_1_dist'][2]            != 'file'):
                            config_dict['psf_shape_1']                         = get_rand_from_dist(config_dict['psf_shape_1_dist'], rng)
                        if(config_dict['psf_shape_2_dist'][2]            != 'file'):
                            config_dict['psf_shape_2']                         = get_rand_from_dist(config_dict['psf_shape_2_dist'], rng)
                        if(config_dict['galaxy_stddev_arcsec_dist'][2]   != 'file'):
                            config_dict['galaxy_stddev_arcsec']                = get_rand_from_dist(config_dict['galaxy_stddev_arcsec_dist'], rng)
                        if(config_dict['galaxy_SN_dist'][2]              != 'file'):
                            config_dict['galaxy_flux_or_SN']                   = get_rand_from_dist(config_dict['galaxy_SN_dist'], rng)
                        if(config_dict['galaxy_shape_1_dist'][2]         != 'file'):
                            config_dict['galaxy_shape_1']                      = get_rand_from_dist(config_dict['galaxy_shape_1_dist'], rng)
                        if(config_dict['galaxy_shear_1_dist'][2]         != 'file'):
                            config_dict['galaxy_shear_1']                      = get_rand_from_dist(config_dict['galaxy_shear_1_dist'], rng)
                        if(config_dict['galaxy_shear_2_dist'][2]         != 'file'):
                            config_dict['galaxy_shear_2']                     = get_rand_from_dist(config_dict['galaxy_shear_2_dist'], rng)
                        if(config_dict['sky_level_subtracted_dist'][2]   != 'file'):
                            config_dict['sky_level_subtracted']                = get_rand_from_dist(config_dict['sky_level_subtracted_dist'], rng)
                        if(config_dict['sky_level_unsubtracted_dist'][2] != 'file'):
                            config_dict['sky_level_unsubtracted']              = get_rand_from_dist(config_dict['sky_level_unsubtracted_dist'], rng)
                        if(config_dict['read_noise_dist'][2]             != 'file'):
                            config_dict['read_noise']                         = get_rand_from_dist(config_dict['read_noise_dist'], rng)
                        if(config_dict['gain_dist'][2]                   != 'file'):
                            config_dict['gain']                               = get_rand_from_dist(config_dict['gain_dist'], rng)
                            
                else:
                    # For the second in the pair, simply adjust the shape angle
                    # The other values will remain the same as for the last image
                    first_in_pair = True
                    config_dict['galaxy_shape_2'] += 90
                    
            else:
                # Regular generation
            
                # For the parameters that are randomly generated once per stamp, generate them now
                if(config_dict['config_mode'] == 'distribution'):
                    if(config_dict['psf_stddev_arcsec_dist'][2]      != 'file'):
                        config_dict['psf_stddev_arcsec']                 = get_rand_from_dist(config_dict['psf_stddev_arcsec_dist'], rng)
                    if(config_dict['psf_shape_1_dist'][2]            != 'file'):
                        config_dict['psf_shape_1']                         = get_rand_from_dist(config_dict['psf_shape_1_dist'], rng)
                    if(config_dict['psf_shape_2_dist'][2]            != 'file'):
                        config_dict['psf_shape_2']                         = get_rand_from_dist(config_dict['psf_shape_2_dist'], rng)
                    if(config_dict['galaxy_stddev_arcsec_dist'][2]   != 'file'):
                        config_dict['galaxy_stddev_arcsec']                = get_rand_from_dist(config_dict['galaxy_stddev_arcsec_dist'], rng)
                    if(config_dict['galaxy_SN_dist'][2]              != 'file'):
                        config_dict['galaxy_flux_or_SN']                           = get_rand_from_dist(config_dict['galaxy_SN_dist'], rng)
                    if(config_dict['galaxy_shape_1_dist'][2]         != 'file'):
                        config_dict['galaxy_shape_1']                      = get_rand_from_dist(config_dict['galaxy_shape_1_dist'], rng)
                    if(config_dict['galaxy_shape_2_dist'][2]         != 'file'):
                        config_dict['galaxy_shape_2']                      = get_rand_from_dist(config_dict['galaxy_shape_2_dist'], rng)
                    if(config_dict['galaxy_shear_1_dist'][2]         != 'file'):
                        config_dict['galaxy_shear_1']                      = get_rand_from_dist(config_dict['galaxy_shear_1_dist'], rng)
                    if(config_dict['galaxy_shear_2_dist'][2]         != 'file'):
                        config_dict['galaxy_shear_2']                     = get_rand_from_dist(config_dict['galaxy_shear_2_dist'], rng)
                    if(config_dict['sky_level_subtracted_dist'][2]   != 'file'):
                        config_dict['sky_level_subtracted']                = get_rand_from_dist(config_dict['sky_level_subtracted_dist'], rng)
                    if(config_dict['sky_level_unsubtracted_dist'][2] != 'file'):
                        config_dict['sky_level_unsubtracted']              = get_rand_from_dist(config_dict['sky_level_unsubtracted_dist'], rng)
                    if(config_dict['read_noise_dist'][2]             != 'file'):
                        config_dict['read_noise']                         = get_rand_from_dist(config_dict['read_noise_dist'], rng)
                    if(config_dict['gain_dist'][2]                   != 'file'):
                        config_dict['gain']                               = get_rand_from_dist(config_dict['gain_dist'], rng)
                
                # Check against hard bounds for the generated values, and enforce them if necessary
                
                # Shears/ellipticities ( quad sum must be < 1 )
                if((config_dict['shear_type'] == 'shear') or (config_dict['shear_type'] == 'ellipticity')):
                    
                    e_tot = np.sqrt(config_dict['psf_shape_1'] ** 2    + config_dict['psf_shape_2'] ** 2   )
                    if(e_tot > shear_max):
                        config_dict['psf_shape_1']    *= shear_max / e_tot
                        config_dict['psf_shape_2']    *= shear_max / e_tot
                        print "WARNING: Hard boundaries had to be enforced for psf shape. If this is happening"
                        print "frequently, consider using a different distribution which will not cross the"
                        print "hard maximum of e or g = " + str(shear_max)
                        
                    e_tot = np.sqrt(config_dict['galaxy_shape_1'] ** 2 + config_dict['galaxy_shape_2'] ** 2)
                    if(e_tot > shear_max):
                        config_dict['galaxy_shape_1'] *= shear_max / e_tot
                        config_dict['galaxy_shape_2'] *= shear_max / e_tot
                        print "WARNING: Hard boundaries had to be enforced for galaxy shape. If this is happening"
                        print "frequently, consider using a different distribution which will not cross the"
                        print "hard maximum of e or g = " + str(shear_max)
                        
                    e_tot = np.sqrt(config_dict['galaxy_shear_1'] ** 2 + config_dict['galaxy_shear_2'] ** 2)
                    if(e_tot > shear_max):
                        config_dict['galaxy_shear_1'] *= shear_max / e_tot
                        config_dict['galaxy_shear_2'] *= shear_max / e_tot
                        print "WARNING: Hard boundaries had to be enforced for galaxy shear. If this is happening"
                        print "frequently, consider using a different distribution which will not cross the"
                        print "hard maximum of e or g = " + str(shear_max)
                else:
                    
                    if(config_dict['psf_shape_1']    > shear_max):
                        config_dict['psf_shape_1']    = shear_max
                        print "WARNING: Hard boundary had to be enforced for psf shape. If this is happening"
                        print "frequently, consider using a different distribution which will not cross the"
                        print "hard maximum of e or g = " + str(shear_max)
                        
                    if(config_dict['galaxy_shape_1'] > shear_max):
                        config_dict['galaxy_shape_1'] = shear_max
                        print "WARNING: Hard boundary had to be enforced for galaxy shape. If this is happening"
                        print "frequently, consider using a different distribution which will not cross the"
                        print "hard maximum of e or g = " + str(shear_max)
                        
                    if(config_dict['galaxy_shear_1'] > shear_max):
                        config_dict['galaxy_shear_1'] = shear_max
                        print "WARNING: Hard boundary had to be enforced for galaxy shear. If this is happening"
                        print "frequently, consider using a different distribution which will not cross the"
                        print "hard maximum of e or g = " + str(shear_max)
                    
                # stddev, SN, and sky level ( must be > 0 )
                if(config_dict['galaxy_stddev_arcsec']   < stddev_min):
                    config_dict['galaxy_stddev_arcsec']      = stddev_min
                    print "WARNING: Hard boundary had to be enforced for galaxy sigma. If this is happening"
                    print "frequently, consider using a different distribution which will not cross the"
                    print "hard minimum of sigma = " + str(stddev_min)
                if(config_dict['psf_stddev_arcsec']      < stddev_min):
                    config_dict['psf_stddev_arcsec']         = stddev_min
                    print "WARNING: Hard boundary had to be enforced for psf sigma. If this is happening"
                    print "frequently, consider using a different distribution which will not cross the"
                    print "hard minimum of sigma = " + str(stddev_min)
                if(config_dict['galaxy_flux_or_SN']              < SN_min):
                    config_dict['galaxy_flux_or_SN']                 = SN_min
                    print "WARNING: Hard boundary had to be enforced for galaxy S/N. If this is happening"
                    print "frequently, consider using a different distribution which will not cross the"
                    print "hard minimum of S/N = " + str(SN_min)
                if(config_dict['sky_level_unsubtracted'] < sky_level_min):
                    config_dict['sky_level_unsubtracted']    = sky_level_min
                    print "WARNING: Hard boundary had to be enforced for unsubtracted sky level. If this is happening"
                    print "frequently, consider using a different distribution which will not cross the"
                    print "hard minimum of unsubtracted sky level = " + str(sky_level_min)
                if(config_dict['sky_level_subtracted']   < sky_level_min):
                    config_dict['sky_level_subtracted']      = sky_level_min
                    print "WARNING: Hard boundary had to be enforced for subtracted sky level. If this is happening"
                    print "frequently, consider using a different distribution which will not cross the"
                    print "hard minimum of subtracted sky level = " + str(sky_level_min)
                if(config_dict['read_noise']             < read_noise_min):
                    config_dict['read_noise']                = read_noise_min
                    print "WARNING: Hard boundary had to be enforced for read noise. If this is happening"
                    print "frequently, consider using a different distribution which will not cross the"
                    print "hard minimum of read noise = " + str(read_noise_min)
                if(config_dict['gain']                   < gain_min):
                    config_dict['gain']                      = gain_min
                    print "WARNING: Hard boundary had to be enforced for gain. If this is happening"
                    print "frequently, consider using a different distribution which will not cross the"
                    print "hard minimum of gain = " + str(gain_min)

            # Calculate the per-pixel sky level
            sky_level_subtracted_pixel   = config_dict['sky_level_subtracted']   * config_dict['sample_scale'] ** 2
            sky_level_unsubtracted_pixel = config_dict['sky_level_unsubtracted'] * config_dict['sample_scale'] ** 2
            
            # Determine the flux (total counts on image) of the galaxy from the S/N if necessary
            if(config_dict['use_flux_or_s/n']=='s/n'):
                gal_flux = get_flux(config_dict['galaxy_flux_or_SN']     ,
                                    config_dict['galaxy_stddev_arcsec']  ,
                                    config_dict['psf_stddev_arcsec']     ,
                                    config_dict['sky_level_subtracted']  ,
                                    config_dict['sky_level_unsubtracted'],
                                    config_dict['read_noise']            ,
                                    config_dict['sample_scale']          ,
                                    config_dict['gain'])
            else:
                gal_flux = config_dict['galaxy_flux_or_SN']

            # Set up the unsheared profiles for the psf and galaxy, both initially as Gaussian
            psf = galsim.Gaussian(flux=1.      , sigma=config_dict['psf_stddev_arcsec']   )  # PSF flux should always = 1
            gal = galsim.Gaussian(flux=gal_flux, sigma=config_dict['galaxy_stddev_arcsec'])
            
            # Apply shearing to the psf and galaxy, with the method determined by the "shear_type" parameter
            if(config_dict['shear_type'] == 'shear'):
                psf.applyShear(g1=config_dict['psf_shape_1']   , g2=config_dict['psf_shape_2']   )
                gal.applyShear(g1=config_dict['galaxy_shape_1'], g2=config_dict['galaxy_shape_2'])
                gal.applyShear(g1=config_dict['galaxy_shear_1'], g2=config_dict['galaxy_shear_2'])
                
            elif(config_dict['shear_type'] == 'ellipticity'):
                psf.applyShear(e1=config_dict['psf_shape_1']   , e2=config_dict['psf_shape_2']   )
                gal.applyShear(e1=config_dict['galaxy_shape_1'], e2=config_dict['galaxy_shape_2'])
                gal.applyShear(e1=config_dict['galaxy_shear_1'], e2=config_dict['galaxy_shear_2'])
                
            elif(config_dict['shear_type'] == 'shear-angle'):
                psf.applyShear(g=config_dict['psf_shape_1']   , beta=config_dict['psf_shape_2']    * galsim.degrees)
                gal.applyShear(g=config_dict['galaxy_shape_1'], beta=config_dict['galaxy_shape_2'] * galsim.degrees)
                gal.applyShear(g=config_dict['galaxy_shear_1'], beta=config_dict['galaxy_shear_2'] * galsim.degrees)
                
            elif(config_dict['shear_type'] == 'ellipticity-angle'):
                psf.applyShear(e=config_dict['psf_shape_1']   , beta=config_dict['psf_shape_2']    * galsim.degrees)
                gal.applyShear(e=config_dict['galaxy_shape_1'], beta=config_dict['galaxy_shape_2'] * galsim.degrees)
                gal.applyShear(e=config_dict['galaxy_shear_1'], beta=config_dict['galaxy_shear_2'] * galsim.degrees)
                
            else:
                print "ERROR: Bad shear type slipped through somehow."
                return
            
            # Set up the pixel profile (a top-hat)
            pix = galsim.Pixel(config_dict['pixel_scale'])

            # Convolve the galaxy, psf, and pixel profile to determine the final (well, before noise) pixelized image
            final = galsim.Convolve([gal, psf, pix])
            
            if(config_dict['grid_or_cube'] == 'cube'):
                # For generating a datacube, set up the image to draw now
                
                if(config_dict['image_type'] == '16i'):
                    gal_image = galsim.ImageC(config_dict['image_size_x_pix'],
                                              config_dict['image_size_y_pix'],
                                              scale=config_dict['sample_scale'])
                    
                elif(config_dict['image_type'] == '32i'):
                    gal_image = galsim.ImageI(config_dict['image_size_x_pix'],
                                              config_dict['image_size_y_pix'],
                                              scale=config_dict['sample_scale'])
                    
                elif(config_dict['image_type'] == '32f'):
                    gal_image = galsim.ImageF(config_dict['image_size_x_pix'],
                                              config_dict['image_size_y_pix'],
                                              scale=config_dict['sample_scale'])
                    
                elif(config_dict['image_type'] == '64f'):
                    gal_image = galsim.ImageD(config_dict['image_size_x_pix'],
                                              config_dict['image_size_y_pix'],
                                              scale=config_dict['sample_scale'])
                    
                else:
                    print "ERROR: Bad image type slipped through somehow."
                    return
            else:
                # For generating a grid, set up the boundaries for this stamp and increment the counter of where in the grid we are
                bounds = galsim.BoundsI(ix      * config_dict['image_size_x_pix'] + ix * config_dict['stamp_padding_x_pix'] + 1,
                                       (ix + 1) * config_dict['image_size_x_pix'] + ix * config_dict['stamp_padding_x_pix']    ,
                                        iy      * config_dict['image_size_y_pix'] + iy * config_dict['stamp_padding_x_pix'] + 1,
                                       (iy + 1) * config_dict['image_size_y_pix'] + iy * config_dict['stamp_padding_x_pix']    ) 
                gal_image = full_image[bounds]
                ix += 1
                if(ix >= config_dict['num_per_row']):
                    # End of the row, so go to the beginning of the next one
                    ix = 0
                    iy += 1
                    
            # Draw the image. Note the commented out version - for some reason GalSim installed differently on a remote machine
            # I used. If this part spits out an error, try switching which of the following two lines is commented out.
            final.draw(gal_image, dx=config_dict['sample_scale'])
            # final.draw(image, scale=config_dict['sample_scale'])
            
            if(config_dict['grid_or_cube'] == 'cube'):
                # If we're doing a datacube, we'll apply the unsubtracted sky level and noise to each image here
                # and add in this image to the list of images
                
                gal_image += sky_level_unsubtracted_pixel
                if(config_dict['suppress_noise'] == 0):
                    noise = galsim.CCDNoise(rng,
                                            gain=config_dict['gain'],
                                            read_noise=config_dict['read_noise'],
                                            sky_level=sky_level_subtracted_pixel)
                    gal_image.addNoise(noise)
                images.append(gal_image)
            
            # Record all data used for this galaxy
            data_string = str(config_dict['psf_stddev_arcsec'])      + '\t' + \
                          str(config_dict['psf_shape_1'])            + '\t' + \
                          str(config_dict['psf_shape_2'])            + '\t' + \
                          str(config_dict['galaxy_stddev_arcsec'])   + '\t' + \
                          str(config_dict['galaxy_flux_or_SN'])              + '\t' + \
                          str(config_dict['galaxy_shape_1'])         + '\t' + \
                          str(config_dict['galaxy_shape_2'])         + '\t' + \
                          str(config_dict['galaxy_shear_1'])         + '\t' + \
                          str(config_dict['galaxy_shear_2'])         + '\t' + \
                          str(config_dict['sky_level_subtracted'])   + '\t' + \
                          str(config_dict['sky_level_unsubtracted']) + '\t' + \
                          str(config_dict['read_noise'])             + '\t' + \
                          str(config_dict['gain'])                   + '\n'
                
            data_strings.append(data_string)
            
            # End generating each individual stamp
            
        # Output the image
        image_file_name = file_name_base + str(i) + '.fits'
        if(config_dict['grid_or_cube'] == 'cube'):
            galsim.fits.writeMulti(images, image_file_name)
            
        else:
            # If we're doing a grid, we'll apply the unsubtracted sky level and noise to the full image
            # here. This will actually add noise to the borders, note, so you won't be able to rely on them
            # being zero.
            full_image += sky_level_unsubtracted_pixel
            if(config_dict['suppress_noise'] == 0):
                noise = galsim.CCDNoise(rng                                 ,
                                        gain=config_dict['gain']            ,
                                        read_noise=config_dict['read_noise'],
                                        sky_level=sky_level_subtracted_pixel)
                full_image.addNoise(noise)
            galsim.fits.write(full_image, image_file_name)
        
        # Output the datafile
        text_file_name = file_name_base + str(i) + '.dat'
        with open(text_file_name, 'w') as data_file:
            # Write the header line first
            data_file.write('# PSF_sigma_arcsec\t'          +
                              'PSF_shape_1\t'               +
                              'PSF_shape_2\t'               + 
                              'Galaxy_sigma_arcsec\t'       + 
                              'Galaxy_SN\tGalaxy_shape_1\t' +
                              'Galaxy_shape_2\t'            +
                              'Galaxy_shear_1\t'            + 
                              'Galaxy_shear_2\t'            +
                              'Sky_level_subtracted\t'      +
                              'Sky_level_unsubtracted\t'    + 
                              'Read_noise\t'                +
                              'Gain\n'                      )
            for s in data_strings:
                data_file.write(s)

        # Allow group access to the files
        cmd = 'chmod g+rw ' + text_file_name
        subprocess.call(cmd, shell=True)
        cmd = 'chmod g+rw ' + image_file_name
        subprocess.call(cmd, shell=True)
        


def bad_config_format():
    """A simple function to output a boilerplate error for a general problem with the configuration file.
    """
    
    print "ERROR: Improperly formatted configuration file. Please check"
    print "sample files and ensure it has one of the proper formats."
    
    # Returns true to alert that there's been an error
    return True


def old_config_version(ver):
    """A boilerplate alert that the config file version being used is old, and some default values
       will have to be assumed.
    """
    print "Using configuration file for version " + ver + " of the script. This is an old version."
    print "Default values will be used for new parameters added in more recent versions."
    return


def new_config_version(ver):
    """A boilerplate note that the config file version being used is up to date.
    """
    print "Using configuration file for version " + ver + " of the script."
    print "This is the most recent version - no default values need to be assumed."
    return


def get_flux(galaxy_SN             ,
             galaxy_stddev_arcsec  ,
             psf_stddev_arcsec     ,
             sky_level_subtracted  ,
             sky_level_unsubtracted,
             read_noise            ,
             pixel_scale           ,
             gain                  ):
    """This function calculates flux from galaxy SN. This depends on the particular definition of SN,
       and will most likely need to be rewritten for the correct definition.
 
       The definition used here is a simplified tophat model, with the tophat containing
       ~50% of the light for a circular galaxy and PSF. This definition is used simply
       because it's easy to analytically invert.
    """

    # # Old algorithm below, which got counts and ADUs a bit mixed up, and used a model for 95% of light,
    # # neglecting PSF. This is retained for reference as a large number of images were created on Hector
    # # using this formula, and this can be used to help rescale the stated S/Ns for those images.
    #
    # size_of_gal = math.pi * 1.96 * (galaxy_stddev_arcsec) ** 2
    # sky_noise_behind_galaxy = math.sqrt( (sky_level_subtracted + sky_level_unsubtracted) * size_of_gal)
    # read_noise_behind_galaxy = (read_noise/gain) *  math.sqrt( size_of_gal ) / pixel_scale
    #
    # background_noise = math.sqrt( sky_noise_behind_galaxy ** 2 + config_dict['read_noise']_behind_galaxy ** 2 )
    #
    # flux = (1 / 0.95) * 0.5 * config_dict['galaxy_flux_or_SN'] * (config_dict['galaxy_flux_or_SN'] + math.sqrt(4*background_noise**2+config_dict['galaxy_flux_or_SN']**2))

    # Estimate the half-light radius of the galaxy (using the magic number 0.674490, which represents the
    # sigma for a Gaussian which contains half the distribution), and use it to calculate the area within
    # the half-light aperture in arcsec
    size_of_gal = math.pi * 0.674490 * ((galaxy_stddev_arcsec) ** 2 + (psf_stddev_arcsec) ** 2)
    
    # Calculate the sky noise and read noise in the half-light aperture, remembering that noise scales with
    # the sqrt of area.
    
    # Sky level is given initially in ADU/arcsec^2. So we convert to counts by multiplying by gain, then
    # take the sqrt of counts to get noise (assuming it's Poisson). Then we scale by sqrt(size_of_gal).  
    sky_noise_behind_galaxy = math.sqrt((sky_level_subtracted + sky_level_unsubtracted) * gain * size_of_gal)
    
    # Read noise is given initially in counts/pixel. So we convert to counts per arcsec (using just pixel_scale,
    # not pixel_scale^2 since it scales with sqrt(area), and then scale by sqrt(size_of_gal)
    read_noise_behind_galaxy = read_noise * math.sqrt(size_of_gal) / pixel_scale
    
    # Total noise is sum in quadrature of the two components
    background_noise = math.sqrt(sky_noise_behind_galaxy ** 2 + read_noise_behind_galaxy ** 2)
    
    # The galaxy's S/N is calculated from both its own Poisson noise and the background noise within its
    # half-light aperture. This can be analytically inverted to give the expression below.
    half_flux = 0.5 * galaxy_SN * (galaxy_SN + math.sqrt(4 * background_noise ** 2 + galaxy_SN ** 2))
    
    flux = half_flux * 2
    
    return flux


def get_dist(type_line   ,
             params_line ,
             mode_line=''):
    
    """This function takes in two or three lines from the configuration file in order to set
       up a tuple containing the necessary information for the distribution from which to
       draw a random variable.
       
       The function adjusts to read in the expected number of parameters for each different
       distribution dist_type. It also performs a few minor corrections: All strings are fixed to
       lower case (so they can be read in in any case), and values that can only ever be
       positive are read in as their absolute values.
    """
    
    # Read in the dist_type of distribution to use
    dist_type = str(type_line.split()[-1]).lower()
    
    if(dist_type == 'fixed'):
        # Fixed to a single value.
        # Parameters: Value
        params = (float(params_line.split()[-1]))
        
    elif(dist_type == 'gaussian'):
        # Gaussian distribution
        # Parameters: Mean, Sigma
        params = (float(params_line.split()[-2])        ,
                  np.abs(float(params_line.split()[-1])))
        
    elif(dist_type == 'truncated_gaussian'):
        # Truncated Gaussian distribution
        # Parameters: Mean, Sigma, Min, Max
        params = (float(params_line.split()[-4])        ,
                  np.abs(float(params_line.split()[-3])),
                  float(params_line.split()[-2])        ,
                  float(params_line.split()[-1])        )
        
    elif(dist_type == 'lognormal-peak'):
        # Gaussian distribution in log10-space
        # Parameters: Peak, Sigma (both in log10-space)
        params = (float(params_line.split()[-2])        ,
                  np.abs(float(params_line.split()[-1])))
        
    elif(dist_type == 'lognormal-mean'):
        # Gaussian distribution in log10-space
        # Parameters: Mean, Sigma (both in log10-space)
        params = (float(params_line.split()[-2])        ,
                  np.abs(float(params_line.split()[-1])))
        
    elif(dist_type == 'rayleigh'):
        # Rayleigh distribution
        # Parameters: Sigma
        params = (np.abs(float(params_line.split()[-1])))
        
    elif(dist_type == 'truncated_rayleigh'):
        # Truncated Rayleigh distribution
        # Parameters: Sigma, Min, Max
        params = (np.abs(float(params_line.split()[-3])),
                  np.abs(float(params_line.split()[-2])),
                  np.abs(float(params_line.split()[-1])))
        
    elif(dist_type == 'contracted_rayleigh'):
        # Contracted Rayleigh distribution (using Bryan's formula to shrink in the tail)
        # Parameters: Sigma, Max, Contraction Power
        params = (np.abs(float(params_line.split()[-3])),
                  np.abs(float(params_line.split()[-2])),
                  np.abs(float(params_line.split()[-1])))
        
    elif(dist_type == 'softened_rayleigh'):
        # Softened Rayleigh distribution (using Rachel and Barney's formula to softly truncate)
        # Parameters: Sigma, Beginning of softening zone, End of softening zone (max)
        # WARNING: The current implementation of this uses GalSim functions, which encounter
        # various glitches for extremely small probablities (ie. sigma << 1). Until this is
        # fixed, be careful about using this method with sigma smaller than 0.2ish.
        params = (np.abs(float(params_line.split()[-3])),
                  np.abs(float(params_line.split()[-2])),
                  np.abs(float(params_line.split()[-1])))
        
    elif(dist_type == 'uniform'):
        # Uniform distribution
        # Parameters: Min, Max
        params = (float(params_line.split()[-2]),
                  float(params_line.split()[-1]))
        
    elif(dist_type == 'loguniform'):
        # Uniform distribution in log10-space
        # Parameters: Min, Max (both in log10-space)
        params = (float(params_line.split()[-2]),
                  float(params_line.split()[-1]))
        
    else:  # Assume uniform otherwise
        print "WARNING: Unrecognised distribution dist_type. Assuming uniform (min, max)."
        print "If this is not correct, halt the program and correct."
        dist_type = 'uniform'
        params = (float(params_line.split()[-2]),
                  float(params_line.split()[-1]))
        
    # Read in the 'mode'
    # This determines whether a variable is generated randomly once per file, or
    # once per postage stamp. Allowed options are 'stamp' and 'file'
    if(mode_line == ''):
        mode = 'stamp' # For backwards-compatibility, when this was the only option
        
    else:
        mode = str(mode_line.split()[-1]).lower()
        
    if((mode != 'file') and (mode != 'stamp')):
        print "WARNING: Unrecognised mode. Assuming 'stamp.' If this is not correct,"
        print "halt the program and correct."
        mode = 'stamp'
        
    return (dist_type, params, mode)


def get_rand_from_dist(dist, last_deviate=None):
    
    """This function generates a random variable from a tuple representing the distribution
       from which to draw the random variable. The distribution tuple should have been set
       up from the get_dist function (see above). Documentation for how to set up distributions
       in the config file can be found in the documentation block for that function.
    """
    
    # Magic numbers
    attempt_max = 10000 # Maximum number of attempts allowed to generate a random number
                        # within truncation boundaries before the algorithm gives up and
                        # returns (min_val+max_val)/2
    
    mode = dist[0]
    if(mode.lower() == 'fixed'):
        return dist[1]
    
    elif(mode.lower() == 'gaussian'):
        mean   = dist[1][0]
        stddev = dist[1][1]
        return np.random.randn() * stddev + mean
    
    elif(mode.lower() == 'truncated_gaussian'):
        
        # Setup
        mean   = dist[1][0]
        stddev = dist[1][1]
        min_val    = dist[1][2]
        max_val    = dist[1][3]
        
        # Sanity checks
        if(min_val > max_val):
            min_val, max_val = max_val, min_val
            
        if(min_val == max_val):
            return min_val
        
        # Try values until one works
        good_value = False
        attempt_counter = 0
        
        while((not good_value)and(attempt_counter < attempt_max)):
            test_result = np.random.randn() * stddev + mean
            if((test_result >= min_val)and(test_result <= max_val)):
                good_value = True
            else:
                attempt_counter += 1
                
        if(good_value):
            # We found a value between min_val and max_val
            return test_result
        else:
            # Failsafe so the program won't crash. This will generally only happen for either extreme
            # min_val and max_val values, for which the user should have chosen a better distribution to model
            # it, or for very small (max_val-min_val), for which the below is a decent return value. Even in
            # the latter case, the user probably should have chosen a better distribution.
            print "WARNING: Failsafe for truncated Gaussian triggered. If this is happening"
            print "frequently, consider using a different distribution."
            return (min_val + max_val) / 2  
        
    elif(mode.lower() == 'lognormal-peak'):
        mean   = dist[1][0]
        stddev = dist[1][1]
        return 10 ** (np.random.randn() * stddev + mean)
    
    elif(mode.lower() == 'lognormal-mean'):
        mean   = dist[1][0]
        stddev = dist[1][1]
        
        # This uses the correction factor I integrated out to adjust to give the correct mean
        return 10 ** (np.random.randn() * stddev + mean) * np.exp(-((stddev * np.log(10)) ** 2) / 2)
    
    elif(mode.lower() == 'rayleigh'):
        sigma = dist[1]
        return np.random.rayleigh(sigma)
    elif(mode.lower() == 'truncated_rayleigh'):
        sigma = dist[1][0]
        min_val   = dist[1][1]
        max_val   = dist[1][2]
        if(min_val > max_val):
            min_val, max_val = max_val, min_val
        if(min_val == max_val):
            return min_val
        good_value = False
        attempt_counter = 0
        while((not good_value) and (attempt_counter < attempt_max)):
            test_result = np.random.rayleigh(sigma)
            if((test_result >= min_val)and(test_result <= max_val)):
                good_value = True
            else:
                attempt_counter += 1
        if(good_value):
            # We found a value between min_val and max_val
            return test_result
        else:
            # Failsafe so the program won't crash. This will generally only happen for either extreme
            # min_val and max_val values, for which the user should have chosen a better distribution to model
            # it, or for very small (max_val-min_val), for which the below is a decent return value. Even in
            # the latter case, the use probably should have chosen a better distribution.
            print "WARNING: Failsafe for truncated Rayleigh triggered. If this is happening"
            print "frequently, consider using a different distribution."
            return (min_val + max_val) / 2  
        
    elif(mode.lower() == 'contracted_rayleigh'):
        sigma = dist[1][0]
        max_val   = dist[1][1]
        p     = dist[1][2]
        
        # Generate an initial random Rayleigh variable
        first_result = np.random.rayleigh(sigma)
        
        # Transform it via Bryan's formula to rein in large values to be less than the max_val
        return (first_result / np.power(1 + np.power(first_result / max_val, p), 1.0 / p))
    
    elif(mode.lower() == 'softened_rayleigh'):
        sigma        = dist[1][0]
        shear_soften = dist[1][1]
        shear_max    = dist[1][2]
        soften_param = 2.*(shear_max - shear_soften)
        
        # To draw from a Rayleigh distribution, we could simply use the GalSim Weibull deviate class with a=2.
        # So we would instantiate the RNG directly and draw from it without specifying a distribution.
        # Once we want to soften the distribution, we do have to specify it directly as a galsim.LookupTable

        # Set up a grid of shear values.
        gvals = np.linspace(0., shear_max, int(1000 * shear_max) + 1)
        
        # Define the Rayleigh distribution.
        pvals = (gvals / sigma ** 2) * np.exp(-gvals ** 2 / (2 * sigma ** 2))
        
        # Apply the softening.
        pvals[gvals > shear_soften] *= np.cos(np.pi * (gvals[gvals > shear_soften] - shear_soften) / soften_param)
        
        # Instantiate the LookupTable.
        gdist_tab = galsim.LookupTable(gvals, pvals)
        
        # Instantiate the random deviate.
        dd = galsim.DistDeviate(last_deviate, gdist_tab)
        
        return dd()
    
    elif(mode.lower() == 'uniform'):
        dist_min = dist[1][0]
        dist_max = dist[1][1]
        return np.random.random() * (dist_max - dist_min) + dist_min
    
    elif(mode.lower() == 'loguniform'):
        dist_min = dist[1][0]
        dist_max = dist[1][1]
        return 10 ** (np.random.random() * (dist_max - dist_min) + dist_min)
    
    else:  # Assume uniform
        print "WARNING: Unrecognized distribution mode. Assuming uniform (min_val, max_val)."
        dist_min = dist[1][0]
        dist_max = dist[1][1]
        return np.random.random() * (dist_max - dist_min) + dist_min
    
    
def load_config_1_0(config_dict, config_lines):
    """Function to load configuration settings for a version 1.0 config file.
    """
    
    # Magic numbers for minimum lengths
    config_simple_length_1_0       = 15
    config_distribution_length_1_0 = 24
    config_min_length_1_0          = 15
    
    old_config_version('1.0')
            
    # Check that it meets the minimum length requirement for this version
    if(len(config_lines) < config_min_length_1_0):
        return bad_config_format()
    
    # Read in bookkeeping/global parameters
    config_dict['config_mode']      =   str(config_lines[1].split()[-1]).lower()
    config_dict['output_name_base'] =   str(config_lines[2].split()[-1])
    config_dict['num_files']        =   int(config_lines[3].split()[-1])
    config_dict['num_per_file']     =   int(config_lines[4].split()[-1])
    config_dict['pixel_scale']      = float(config_lines[5].split()[-1])
    config_dict['image_size_x_pix'] =   int(config_lines[6].split()[-1])
    config_dict['image_size_y_pix'] =   int(config_lines[7].split()[-1])
    
    # Check whether we're using the simple or complex config_mode
    if(config_dict['config_mode'] == 'simple'):
        if(len(config_lines) < config_simple_length_1_0):
            return bad_config_format()
        
        # Using simple mode, so we'll just use fixed values for all parameters.
        # Only the noise realisation will differ between images.
        
        config_dict['galaxy_stddev_arcsec'] = float(config_lines[ 8].split()[-1])
        config_dict['psf_stddev_arcsec']    = float(config_lines[ 9].split()[-1])
        config_dict['galaxy_flux_or_SN']            = float(config_lines[10].split()[-1])
        config_dict['galaxy_shape_1']       = float(config_lines[11].split()[-1])
        config_dict['galaxy_shape_2']         = float(config_lines[12].split()[-1])
        config_dict['galaxy_shear_1']         = float(config_lines[13].split()[-1])
        config_dict['galaxy_shear_2']         = float(config_lines[14].split()[-1])
        config_dict['sky_level_subtracted'] = float(config_lines[15].split()[-1])
        
    elif(config_dict['config_mode'] == 'distribution'):
        if(len(config_lines) < config_distribution_length_1_0):
            return bad_config_format()
        
        # Using distribution mode, so we'll read in using the get_dist function,
        # which determines the distribution type from the first line and the
        # parameters for it from the second.
        
        config_dict['galaxy_stddev_arcsec_dist'] = \
            get_dist(config_lines[8], config_lines[9])
        config_dict['psf_stddev_arcsec_dist'] = \
            get_dist(config_lines[10], config_lines[11])
        config_dict['galaxy_SN_dist'] = \
            get_dist(config_lines[12], config_lines[13])
        config_dict['galaxy_shape_1_dist'] = \
            get_dist(config_lines[14], config_lines[15])
        config_dict['galaxy_shape_2_dist'] = \
            get_dist(config_lines[16], config_lines[17])
        config_dict['galaxy_shear_1_dist'] = \
            get_dist(config_lines[18], config_lines[19])
        config_dict['galaxy_shear_2_dist'] = \
            get_dist(config_lines[20], config_lines[21])
        config_dict['sky_level_subtracted_dist'] = \
            get_dist(config_lines[22], config_lines[23])
    else:
        return bad_config_format()
    
    # Fill in values unsupported by this version
    load_back_config_1_0(config_dict)
    
    return False
    
    
def load_config_1_1(config_dict, config_lines):
    """Function to load configuration settings for a version 1.1 config file.
    """
    
    # Magic numbers for minimum lengths
    config_simple_length_1_1       = 22
    config_distribution_length_1_1 = 32
    config_min_length_1_1          = 22
    
    old_config_version('1.1')
    
    # Check that it meets the minimum length requirement for this version
    if(len(config_lines) < config_min_length_1_1):
        return bad_config_format()
    
    config_dict['config_mode']      =   str(config_lines[1].split()[-1]).lower()
    config_dict['output_name_base'] =   str(config_lines[2].split()[-1])
    config_dict['num_files']        =   int(config_lines[3].split()[-1])
    config_dict['num_per_file']     =   int(config_lines[4].split()[-1])
    config_dict['pixel_scale']      = float(config_lines[5].split()[-1])
    config_dict['image_size_x_pix'] =   int(config_lines[6].split()[-1])
    config_dict['image_size_y_pix'] =   int(config_lines[7].split()[-1])
    config_dict['image_type']       =   str(config_lines[8].split()[-1]).lower()
    config_dict['shear_type']       =   str(config_lines[9].split()[-1]).lower()
    
    # Check whether we're using the simple or complex config_mode                
    if(config_dict['config_mode'] == 'simple'):
        if(len(config_lines) < config_simple_length_1_1):
            return bad_config_format()
        
        # Using simple mode, so we'll just use fixed values for all parameters.
        # Only the noise realisation will differ between images.
        
        config_dict['galaxy_stddev_arcsec']    = float(config_lines[10].split()[-1])
        config_dict['psf_stddev_arcsec']       = float(config_lines[11].split()[-1])
        config_dict['galaxy_flux_or_SN']               = float(config_lines[12].split()[-1])
        config_dict['galaxy_shape_1']          = float(config_lines[13].split()[-1])
        config_dict['galaxy_shape_2']          = float(config_lines[14].split()[-1])
        config_dict['galaxy_shear_1']          = float(config_lines[15].split()[-1])
        config_dict['galaxy_shear_2']          = float(config_lines[16].split()[-1])
        config_dict['sky_level_subtracted']    = float(config_lines[17].split()[-1])
        config_dict['sky_level_unsubtracted']  = float(config_lines[18].split()[-1])
        config_dict['read_noise']              = float(config_lines[19].split()[-1])
        config_dict['init_random_seed_factor'] =   int(config_lines[20].split()[-1])
        config_dict['suppress_noise']          =   int(config_lines[21].split()[-1])
        
    elif(config_dict['config_mode'] == 'distribution'):
        if(len(config_lines) < config_distribution_length_1_1):
            return bad_config_format()
        
        # Using distribution mode, so we'll read in using the get_dist function,
        # which determines the distribution type from the first line and the
        # parameters for it from the second.
        
        config_dict['galaxy_stddev_arcsec_dist'] = \
            get_dist(config_lines[10], config_lines[11])
        config_dict['psf_stddev_arcsec_dist'] = \
            get_dist(config_lines[12], config_lines[13])
        config_dict['galaxy_SN_dist'] = \
            get_dist(config_lines[14], config_lines[15])
        config_dict['galaxy_shape_1_dist'] = \
            get_dist(config_lines[16], config_lines[17])
        config_dict['galaxy_shape_2_dist'] = \
            get_dist(config_lines[18], config_lines[19])
        config_dict['galaxy_shear_1_dist'] = \
            get_dist(config_lines[20], config_lines[21])
        config_dict['galaxy_shear_2_dist'] = \
            get_dist(config_lines[22], config_lines[23])
        config_dict['sky_level_subtracted_dist'] = \
            get_dist(config_lines[24], config_lines[25])
        config_dict['sky_level_unsubtracted_dist'] = \
            get_dist(config_lines[26], config_lines[27])
        config_dict['read_noise_dist'] = \
            get_dist(config_lines[28], config_lines[29])
            
        config_dict['init_random_seed_factor'] = int(config_lines[30].split()[-1])
        config_dict['suppress_noise']          = int(config_lines[31].split()[-1])

    else:
        return bad_config_format()
    
    # Fill in values unsupported by this version
    load_back_config_1_1(config_dict)
    
    
def load_config_1_2(config_dict, config_lines):
    """Function to load configuration settings for a version 1.2 config file.
    """
    
    # Magic numbers for minimum lengths
    config_simple_length_1_2       = 32
    config_distribution_length_1_2 = 58
    config_min_length_1_2          = 32
    
    old_config_version('1.2')
            
    # Check that it meets the minimum length requirement for this version
    if(len(config_lines) < config_min_length_1_2):
        return bad_config_format()
    
    # Read in bookkeeping parameters
    config_dict['config_mode']             =   str(config_lines[ 1].split()[-1]).lower()
    config_dict['output_folder_name']      =   str(config_lines[ 2].split()[-1])
    config_dict['thread_folder_name_base'] =   str(config_lines[ 3].split()[-1])
    config_dict['output_name_base']        =   str(config_lines[ 4].split()[-1])
    config_dict['num_threads']             =   int(config_lines[ 5].split()[-1])
    config_dict['num_files']               =   int(config_lines[ 6].split()[-1])
    config_dict['grid_or_cube']            =   str(config_lines[ 7].split()[-1]).lower()
    if(config_dict['grid_or_cube'] == 'cube'):
        config_dict['num_per_file']        =   int(config_lines[ 8].split()[-1])
    elif(config_dict['grid_or_cube'] == 'grid'):
        config_dict['num_per_row']         =   int(config_lines[ 8].split()[-2])
        config_dict['num_per_col']         =   int(config_lines[ 8].split()[-1])
        config_dict['num_per_file']        =   config_dict['num_per_row'] * config_dict['num_per_col']
    else:
        return bad_config_format()
    config_dict['pixel_scale']             = float(config_lines[ 9].split()[-1])
    config_dict['sample_scale']            = float(config_lines[10].split()[-1])
    config_dict['image_size_x_pix']        =   int(config_lines[11].split()[-1])
    config_dict['image_size_y_pix']        =   int(config_lines[12].split()[-1])
    config_dict['stamp_padding_x_pix']     =   int(config_lines[13].split()[-1])
    config_dict['stamp_padding_y_pix']     =   int(config_lines[14].split()[-1])
    config_dict['image_type']                      =   str(config_lines[15].split()[-1]).lower()
    config_dict['shear_type']                      =   str(config_lines[16].split()[-1]).lower()
    config_dict['init_random_seed_factor'] =   int(config_lines[17].split()[-1])
    config_dict['suppress_noise']          =   int(config_lines[18].split()[-1])
    
    # Check whether we're using the simple or complex config_dict['config_mode']
    if(config_dict['config_mode'] == 'simple'):
        if(len(config_lines) < config_simple_length_1_2):
            return bad_config_format()
        
        # Using simple mode, so we'll just use fixed values for all parameters.
        # Only the noise realisation will differ between images.
        
        config_dict['psf_stddev_arcsec']      = float(config_lines[19].split()[-1])
        config_dict['psf_shape_1']            = float(config_lines[20].split()[-1])
        config_dict['psf_shape_2']            = float(config_lines[21].split()[-1])
        config_dict['galaxy_stddev_arcsec']   = float(config_lines[22].split()[-1])
        config_dict['galaxy_flux_or_SN']              = float(config_lines[23].split()[-1])
        config_dict['galaxy_shape_1']         = float(config_lines[24].split()[-1])
        config_dict['galaxy_shape_2']         = float(config_lines[25].split()[-1])
        config_dict['galaxy_shear_1']         = float(config_lines[26].split()[-1])
        config_dict['galaxy_shear_2']         = float(config_lines[27].split()[-1])
        config_dict['sky_level_subtracted']   = float(config_lines[28].split()[-1])
        config_dict['sky_level_unsubtracted'] = float(config_lines[29].split()[-1])
        config_dict['read_noise']             = float(config_lines[30].split()[-1])
        config_dict['gain']                   = float(config_lines[31].split()[-1])
        
    elif(config_dict['config_mode'] == 'distribution'):
        if(len(config_lines) < config_distribution_length_1_2):
            return bad_config_format()
        
        # Using distribution mode, so we'll read in using the get_dist function,
        # which determines the distribution type from the first line and the
        # parameters for it from the second.
        
        config_dict['psf_stddev_arcsec_dist'] = \
            get_dist(config_lines[20], config_lines[21], config_lines[19])
        config_dict['psf_shape_1_dist'] = \
            get_dist(config_lines[23], config_lines[24], config_lines[22])
        config_dict['psf_shape_2_dist'] = \
            get_dist(config_lines[26], config_lines[27], config_lines[25])
        config_dict['galaxy_stddev_arcsec_dist'] = \
            get_dist(config_lines[29], config_lines[30], config_lines[28])
        config_dict['galaxy_SN_dist'] = \
            get_dist(config_lines[32], config_lines[33], config_lines[31])
        config_dict['galaxy_shape_1_dist'] = \
            get_dist(config_lines[35], config_lines[36], config_lines[34])
        config_dict['galaxy_shape_2_dist'] = \
            get_dist(config_lines[38], config_lines[39], config_lines[37])
        config_dict['galaxy_shear_1_dist'] = \
            get_dist(config_lines[41], config_lines[42], config_lines[40])
        config_dict['galaxy_shear_2_dist'] = \
            get_dist(config_lines[44], config_lines[45], config_lines[43])
        config_dict['sky_level_subtracted_dist'] = \
            get_dist(config_lines[47], config_lines[48], config_lines[46])
        config_dict['sky_level_unsubtracted_dist'] = \
            get_dist(config_lines[50], config_lines[51], config_lines[49])
        config_dict['read_noise_dist'] = \
            get_dist(config_lines[53], config_lines[54], config_lines[52])
        config_dict['gain_dist'] = \
            get_dist(config_lines[56], config_lines[57], config_lines[55])
        
    else:
        return bad_config_format()
    
    # Fill in values unsupported by this version
    load_back_config_1_2(config_dict)
    
    
def load_config_1_3(config_dict, config_lines):
    """Function to load configuration settings for a version 1.3 config file.
    """
    
    # Magic numbers for minimum lengths
    config_simple_length_1_3       = 34
    config_distribution_length_1_3 = 60
    config_min_length_1_3          = 34
            
    new_config_version('1.3')
    
    # Check that it meets the minimum length requirement for this version
    if(len(config_lines) < config_min_length_1_3):
        return bad_config_format()
    
    # Initialize counter for which line we're on
    # This is used to make it easier to add lines in future versions
    lc = 1
    
    # Read in bookkeeping parameters
    config_dict['config_mode']             =   str(config_lines[lc].split()[-1]).lower()
    lc += 1
    config_dict['output_folder_name']      =   str(config_lines[lc].split()[-1])
    lc += 1
    config_dict['thread_folder_name_base'] =   str(config_lines[lc].split()[-1])
    lc += 1
    config_dict['output_name_base']        =   str(config_lines[lc].split()[-1])
    lc += 1
    config_dict['num_threads']             =   int(config_lines[lc].split()[-1])
    lc += 1
    config_dict['num_files']               =   int(config_lines[lc].split()[-1])
    lc += 1
    config_dict['grid_or_cube']            =   str(config_lines[lc].split()[-1]).lower()
    lc += 1
    if(config_dict['grid_or_cube'] == 'cube'):
        config_dict['num_per_file']        =   int(config_lines[lc].split()[-1])
    elif(config_dict['grid_or_cube'] == 'grid'):
        config_dict['num_per_row']         =   int(config_lines[lc].split()[-2])
        config_dict['num_per_col']         =   int(config_lines[lc].split()[-1])
        config_dict['num_per_file']        =   config_dict['num_per_row'] * config_dict['num_per_col']
    else:
        return bad_config_format()
    lc += 1
    config_dict['pixel_scale']             = float(config_lines[lc].split()[-1])
    lc += 1
    config_dict['sample_scale']            = float(config_lines[lc].split()[-1])
    lc += 1
    config_dict['image_size_x_pix']        =   int(config_lines[lc].split()[-1])
    lc += 1
    config_dict['image_size_y_pix']        =   int(config_lines[lc].split()[-1])
    lc += 1
    config_dict['stamp_padding_x_pix']     =   int(config_lines[lc].split()[-1])
    lc += 1
    config_dict['stamp_padding_y_pix']     =   int(config_lines[lc].split()[-1])
    lc += 1
    config_dict['image_type']              =   str(config_lines[lc].split()[-1]).lower()
    lc += 1
    config_dict['shear_type']              =   str(config_lines[lc].split()[-1]).lower()
    lc += 1
    config_dict['use_flux_or_s/n']         =   str(config_lines[lc].split()[-1]).lower()
    lc += 1
    config_dict['init_random_seed_factor'] =   int(config_lines[lc].split()[-1])
    lc += 1
    config_dict['suppress_noise']          =   int(config_lines[lc].split()[-1])
    lc += 1
    config_dict['shape_noise_cancellation']=   int(config_lines[lc].split()[-1])
    lc += 1
    
    # Check if other parameters are consistent with shape noise cancellation if it's enabled
    if((config_dict['shape_noise_cancellation']!=0) and
       (config_dict['shear_type']!='shear-angle') and 
       (config_dict['shear_type']!='ellipticity-angle') ):
        if(config_dict['shear_type']!='ellipticity'):
            config_dict['shear_type'] = 'ellipticity-angle'
        else:
            config_dict['shear_type'] = 'shear-angle'
        print("ERROR: Shear type is incompatible with shape noise cancellation.")
        print("It has been changed to " +  config_dict['shear_type'] + ". Be warned that this")
        print("will likely result in improper values for shapes and shears. Either disable")
        print("shape noise cancellation or correct it to one of the '-angle' types.")
        
    if((config_dict['shape_noise_cancellation']!=0) and
       (config_dict['num_per_file'] % 2 == 1)):
        print("WARNING: Shape noise cancellation requires an even number of images per")
        print("file. The number of images per file has been doubled to compensate.")
        if(config_dict['grid_or_cube'] == 'cube'):
            config_dict['num_per_file'] *= 2
        elif(config_dict['grid_or_cube'] == 'grid'):
            config_dict['num_per_row'] *= 2
            config_dict['num_per_file'] *= 2
        
        
    # Check whether we're using the simple or complex config_dict['config_mode']
    if(config_dict['config_mode'] == 'simple'):
        if(len(config_lines) < config_simple_length_1_3):
            return bad_config_format()
        
        # Using simple mode, so we'll just use fixed values for all parameters.
        # Only the noise realisation will differ between images.
        
        config_dict['psf_stddev_arcsec']      = float(config_lines[21].split()[-1])
        lc += 1
        config_dict['psf_shape_1']            = float(config_lines[22].split()[-1])
        lc += 1
        config_dict['psf_shape_2']            = float(config_lines[23].split()[-1])
        lc += 1
        config_dict['galaxy_stddev_arcsec']   = float(config_lines[24].split()[-1])
        lc += 1
        config_dict['galaxy_flux_or_SN']      = float(config_lines[25].split()[-1])
        lc += 1
        config_dict['galaxy_shape_1']         = float(config_lines[26].split()[-1])
        lc += 1
        config_dict['galaxy_shape_2']         = float(config_lines[27].split()[-1])
        lc += 1
        config_dict['galaxy_shear_1']         = float(config_lines[28].split()[-1])
        lc += 1
        config_dict['galaxy_shear_2']         = float(config_lines[29].split()[-1])
        lc += 1
        config_dict['sky_level_subtracted']   = float(config_lines[30].split()[-1])
        lc += 1
        config_dict['sky_level_unsubtracted'] = float(config_lines[31].split()[-1])
        lc += 1
        config_dict['read_noise']             = float(config_lines[32].split()[-1])
        lc += 1
        config_dict['gain']                   = float(config_lines[33].split()[-1])
        lc += 1
        
    elif(config_dict['config_mode'] == 'distribution'):
        if(len(config_lines) < config_distribution_length_1_3):
            return bad_config_format()
        
        # Using distribution mode, so we'll read in using the get_dist function,
        # which determines the distribution type from the first line and the
        # parameters for it from the second.
        
        config_dict['psf_stddev_arcsec_dist'] = \
            get_dist(config_lines[lc+1], config_lines[lc+2], config_lines[lc])
        lc += 3
        config_dict['psf_shape_1_dist'] = \
            get_dist(config_lines[lc+1], config_lines[lc+2], config_lines[lc])
        lc += 3
        config_dict['psf_shape_2_dist'] = \
            get_dist(config_lines[lc+1], config_lines[lc+2], config_lines[lc])
        lc += 3
        config_dict['galaxy_stddev_arcsec_dist'] = \
            get_dist(config_lines[lc+1], config_lines[lc+2], config_lines[lc])
        lc += 3
        config_dict['galaxy_SN_dist'] = \
            get_dist(config_lines[lc+1], config_lines[lc+2], config_lines[lc])
        lc += 3
        config_dict['galaxy_shape_1_dist'] = \
            get_dist(config_lines[lc+1], config_lines[lc+2], config_lines[lc])
        lc += 3
        config_dict['galaxy_shape_2_dist'] = \
            get_dist(config_lines[lc+1], config_lines[lc+2], config_lines[lc])
        lc += 3
        config_dict['galaxy_shear_1_dist'] = \
            get_dist(config_lines[lc+1], config_lines[lc+2], config_lines[lc])
        lc += 3
        config_dict['galaxy_shear_2_dist'] = \
            get_dist(config_lines[lc+1], config_lines[lc+2], config_lines[lc])
        lc += 3
        config_dict['sky_level_subtracted_dist'] = \
            get_dist(config_lines[lc+1], config_lines[lc+2], config_lines[lc])
        lc += 3
        config_dict['sky_level_unsubtracted_dist'] = \
            get_dist(config_lines[lc+1], config_lines[lc+2], config_lines[lc])
        lc += 3
        config_dict['read_noise_dist'] = \
            get_dist(config_lines[lc+1], config_lines[lc+2], config_lines[lc])
        lc += 3
        config_dict['gain_dist'] = \
            get_dist(config_lines[lc+1], config_lines[lc+2], config_lines[lc])
        lc += 3
            
        if(config_dict['shape_noise_cancellation']!=0):
            # Enforce that galaxy and psf shape are generated once per file
            config_dict['psf_shape_1_dist'] = (config_dict['psf_shape_1_dist'][0],
                                                  config_dict['psf_shape_1_dist'][1],
                                                  'file')
            config_dict['psf_shape_2_dist'] = (config_dict['psf_shape_2_dist'][0],
                                                  config_dict['psf_shape_2_dist'][1],
                                                  'file')
            config_dict['galaxy_shape_1_dist'] = (config_dict['galaxy_shape_1_dist'][0],
                                                  config_dict['galaxy_shape_1_dist'][1],
                                                  'file')
            config_dict['galaxy_shape_2_dist'] = (config_dict['galaxy_shape_2_dist'][0],
                                                  config_dict['galaxy_shape_2_dist'][1],
                                                  'file')
    else:
        return bad_config_format()   
    
def load_back_config_1_0(config_dict):
    """Function to load configuration settings that were introduced after version 1.0
    """
    
    config_dict['shear_type']                  = default_shear_type
    config_dict['image_type']                  = default_image_type
    config_dict['read_noise']             = 0
    config_dict['sky_level_unsubtracted'] = 0
    config_dict['read_noise_dist']             = ('fixed', 0, 'file')
    config_dict['sky_level_unsubtracted_dist'] = ('fixed', 0, 'file')
    config_dict['init_random_seed_factor']     = 1
    config_dict['suppress_noise']              = 0
    
    load_back_config_1_1(config_dict)
    
    
def load_back_config_1_1(config_dict):
    """Function to load configuration settings that were introduced after version 1.1
    """
    config_dict['output_folder_name']      = "./"
    config_dict['thread_folder_name_base'] = "thread"
    config_dict['num_threads']             = 1
    config_dict['sample_scale']            = config_dict['pixel_scale']
    config_dict['grid_or_cube']            = "grid"
    config_dict['stamp_padding_x_pix']     = 1
    config_dict['stamp_padding_y_pix']     = 1
    config_dict['psf_shape_1']             = 0
    config_dict['psf_shape_2']             = 0
    config_dict['gain']                    = default_gain
    config_dict['psf_shape_1_dist']        = ('fixed', 0, 'file')
    config_dict['psf_shape_2_dist']        = ('fixed', 0, 'file')
    config_dict['gain_dist']               = ('fixed', default_gain, 'file')
    
    load_back_config_1_2(config_dict)
    
    
def load_back_config_1_2(config_dict):
    """Function to load configuration settings that were introduced after version 1.2
    """
    config_dict['shape_noise_cancellation'] = 0
    config_dict['use_flux_or_s/n'] = 's/n'
    load_back_config_1_3(config_dict)
    return
    
    
def load_back_config_1_3(config_dict):
    """Function to load configuration settings that were introduced after version 1.3
    """
    return
    

if __name__ == "__main__":
    main(sys.argv)
