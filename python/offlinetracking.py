"""
This routine conducts offline particle tracking based on imported 3D velocity field.



To Do
--------
4) translate advectparticles

5) Translate backward tracking.
6) Translate dynamic seeding
"""

import pandas as pd
import numpy as np
from saveconfiguration import saveconfiguration
from getvelocityfield import getvelocityfield
from iniparticles import iniparticles
from getpartivar import getpartivar
from saveparti import saveparti
from advectparticles import advectparticles
from advectparticles_parallel import advectparticles_parallel
import sys
sys.path.insert(0, 'testcase')
from myvelocityfield import myvelocityfield
import datetime
from joblib import Parallel, delayed
#import multiprocessing

## PARAMETERS
#
model = {} # create pandas Series object
#####################
## PSOM parameters
#####################
# Day of year where simulation starts
model['start_day'] = 0 
# Model timestep in second (needed to convert from timestamp to DOY
model['timestep'] = 216 
# Path of folder containing model outputs
#model['path'] = '/Volumes/garuda/Mathieu/output_Papa1km/particle_tracking/' 
model['path'] = '/Volumes/garuda-1/Mathieu/output_Papa_Dx=500_from_1km/K1/highrez/' 
# constant to play with model output temporal resolution
# "2" means every other model output file
# "4" means every 4 model output file
# "1" means every file
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

# day of year for particle initialization (i.e., first seeding)
particle.initime = 60.5
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
particle.istart = 50000
particle.irange = 1000
particle.irez = 2000

particle.jstart = 150000
particle.jrange = 1000
particle.jrez = 2000

particle.kstart = -10
particle.krange = 1
particle.krez = 3

#----------------
# ADVECTION
#----------------
# timestep of particle tracking in days
particle.timestep = 2**-4
# length of particle tracking experiment in days
particle.length = 5
# Define particle-tracking direction ('forward' or 'backward')
particle.direction = 'forward' 

#----------------
# OUTPUT
#----------------
# frequency of particles output (in days)
particle.outfreq = 2**-3
# Format of output ('csv' or 'sqlite')
particle.outputformat = 'sqlite' 

particle.outputdir = '/Users/mathieudever/Desktop/particle_test/PYTHON' 
particle.outputfilename = 'pyton_test_1parti' 


#####################
## CORE CODE
#####################
time_total = datetime.datetime.now()

# Change a few things if backward-tracking
if particle.direction=='backward':
    particle.timestep = -particle.timestep
    particle.length = -particle.length

# Print the model configuration to a text-file
saveconfiguration(model,particle)


#pool = multiprocessing.Pool()
#%%
# Record the total number of advective timesteps
Titer = np.arange(particle.initime,particle.initime+particle.length+particle.timestep,particle.timestep).size

# Main Loop
for tt in np.arange(particle.initime,particle.initime+particle.length+particle.timestep,particle.timestep):
    print('\n')
    print('%2.2f %% DONE (doy %d)' %( (tt-particle.initime)/(particle.length)*100,tt) ) 
    
    # Get current time with datetime
    time_loop = datetime.datetime.now()
    
    ## Get the velocity field relevant to the time step considered, and 
    # already interpolated onto the time step
    model = getvelocityfield(model,tt)
   

    ## Uncomment to run test case with homogeneous flow
    model  = myvelocityfield(model)
    
    ## Advect particles backward in time if backward particle tracking
    
    # Skip a time step in the case of backward particle tracking. this
    # is because the forward case uses the velocity field at time t to
    # compute trajectories from t to t+Dt, so backward tracking should
    # use the same velocity field to compute the trajectory from t+Dt
    # to t.
    """
    if particle.direction=='backward' and tt == particle.initime:
        None # Skip a time step
    elif particle.direction=='backward':
        None #run advectparticles
    """
    
    ## Seed particle when required
    
    # If day of year matches the seeding frequency (threshold is to account
    # for rouding error is inifreq is not a power of 2).
    if (tt-particle.initime)%particle.inifreq<=1e-10 and particle.ininumber != 0:

        # If seeding is set to 'dynamic', the iniparticle routine uses the
        # current particle position of the 1st particle to seed (at a 
        # constant depth)
        
        """
        if particle.initype=='dynamic') and tt != particle.initime:
            particle.istart = parti.x(parti.id==1) 
            particle.irange = 1 
            particle.irez = 1 
            
            particle.jstart = parti.y(parti.id==1) 
            particle.jrange = 1 
            particle.jrez = 1 
            
            particle.kstart = -50 
            particle.krange = 1 
            particle.krez = 1 
        """
        
        # Seed particles
        try:
            parti
        except:
            parti=pd.DataFrame([])
        parti = iniparticles(tt,particle,parti)
        particle.ininumber = particle.ininumber -1 
    
    # printlay total number of particle in the experiment
    print('Total number of particles is %d.' %len(parti.x))
    
    ## Get the model variables interpolated onto the particle's position
    # Interpolate model variables onto particle's position in 4D
    parti = getpartivar(tt,parti,model)
    
    ## Save particle position when required

    if (tt-particle.initime)%particle.outfreq<=1e-10:
        saveparti(tt,particle,parti)
    
    ## Advect particles forward in time if forward particle tracking
    if particle.direction=='forward':
        
        ## Unparallel advection
        #parti = advectparticles(tt,particle,model,parti)

        ## Parallel advection
        print('PROJECT PARTICLE DATA... ')
        Parallel(n_jobs=8, prefer="threads")(delayed(advectparticles_parallel)(ii,parti,particle.timestep,particle.direction,model) for ii in range(0,len(parti.x)))
        print('Removed %d particles' % np.count_nonzero(parti.isnull()))
        if np.count_nonzero(parti.isnull()) == parti.shape[0]:
            print('END-OF-RUN No particle left to track')
            parti = parti.dropna()
            exit(0)
        else:
            parti = parti.dropna()

    ## printlays some info to user
    
    # End timer
    # Displays time for each time loop and projects the remaining time
    telapsed = datetime.datetime.now()-time_loop 
    print('This step took %s minutes' %(str( telapsed ).split('.')[0]))
    print('At this rate, you have %3.2f minutes left' %(telapsed.seconds*(Titer-tt)/60))
    
totaltime = datetime.datetime.now()-time_total
print('Entire run took %s minutes' %(str( totaltime ).split('.')[0])) 
