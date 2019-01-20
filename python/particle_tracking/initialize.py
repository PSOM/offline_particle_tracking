#
# contains multiple initialize files
# make your own for your project.

def initialize_default():
    ''' The function sets all the important parameters.
    '''
    
    
    import warnings
    import os,sys
    from time import time
    import pandas as pd
    import numpy as np
    import xarray as xr
    import datetime
    
    model = {} # create dictionary object
    
    ## PARAMETERS
    
    #####################
    ## PSOM parameters
    #####################
    
    # Enable more descriptive output
    model['verbose'] = True
    if ~model['verbose']:
        warnings.simplefilter('ignore')
    # Enable parallelized particle loop
    model['parallel'] = False
    # Day of year where simulation starts
    model['start_day'] = 1
    # Model timestep in second (needed to convert from timestamp to DOY
    model['timestep'] = 288
    model['delta_x'] = 1000 # meters
    model['delta_y'] = 1000 # meters
    # Path of folder containing model outputs
    model['path'] = '../../psom_output/output_2013asiri_06b/' 
    # constant to play with model output temporal resolution
    # "2" means every "2"nd model output file
    model['outputskip'] = 1 

    # Define periodicities in model domain
    model['periodic_ew'] = 1 
    model['periodic_ns'] = 0 

    #####################
    ## Particle parameters
    #####################

    #----------------
    # SEEDING
    #----------------
    particle = pd.DataFrame([])# create pandas Series object

    # Enable depth-keeping particles
    particle.isobaric = True
    
    # day of year for particle initialization (i.e., first seeding)
    particle.initime = 10
    # frequency of particles seeding (in days)
    particle.inifreq = 2**-2 
    # Number of seeding events
    particle.ininumber = 1
    # Type of particle seeding ('dynamic' or 'static')
    particle.initype = 'static' 

    # Number of particle classes to be released. Please refer to iniparticles
    # to customize the sinking velocity of each particle class
    particle.numofclasses = 1

    # Parameters to customize the initialization of the particles (in meters)
    # Only applied if seeding type is 'static'
    particle.istart = 188000
    particle.irange = 2000
    particle.irez = 10

    particle.jstart = 100000
    particle.jrange = 2000
    particle.jrez = 10

    particle.kstart = -100
    particle.krange = 100
    particle.krez = 1

    #----------------
    # ADVECTION
    #----------------
    # timestep of particle tracking in days
    particle.timestep = 2**-4
    # length of particle tracking experiment in days
    particle.length = 15
    # Define particle-tracking direction ('forward' or 'backward')
    particle.direction = 'forward' 

    #----------------
    # OUTPUT
    #----------------
    # frequency of particles output (in days)
    particle.outfreq = 2**-3
    # Format of output ('csv' or 'sqlite')
    particle.outputformat = 'csv' 

    particle.outputdir = './committee_meeting_3d/' 
    particle.outputfilename = 'vertical_parti%d_dk' %np.int(sys.argv[1])

    if os.path.exists(os.path.join(particle.outputdir,particle.outputfilename)):
        sys.exit('Warning file already exists!')


    #####################
    ## CORE CODE
    #####################
    time_total = datetime.datetime.now()

    # Change a few things if backward-tracking
    if particle.direction=='backward':
        particle.timestep = -particle.timestep
        particle.length = -particle.length
        
    return model, particle