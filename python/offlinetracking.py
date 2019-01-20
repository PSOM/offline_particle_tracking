"""
This routine conducts offline particle tracking based on imported 3D velocity field.

"""

# python packages
import os,sys
import pandas as pd
import numpy as np
import xarray as xr
import dask.dataframe as dd
from time import time
import warnings

# custom functions
from initialize_comm import __initialize__
from particle_tracking.model_module import *
from particle_tracking.parti_module import *
# select initialization file for your experiment:
from particle_tracking.initialize import initialize_experiment00 as __initialize__

# start timer for full script
time_total = time()

# set parameters in model and particle
model, particle = __initialize__()

# Print the model configuration to a text-file
saveconfiguration(model,particle)

# Record the total number of advective timesteps
Titer = np.arange(particle.initime,particle.initime+particle.length+particle.timestep,particle.timestep).size

# Main Loop
iteration = [particle.initime+t*particle.timestep for t in range(int( particle.length/particle.timestep ))]
for tcount,tt in enumerate(iteration):
    # start loop timer
    time_loop = time()
    #print('\n')
    print('%2.2f %% DONE (doy %3.2f)' %( (tt-particle.initime)/(particle.length)*100, tt)  ) 
    
    ## Get the velocity field relevant to the time step considered, and 
    # already interpolated onto the time step
    model = getvelocityfield(model,tt)
    # create 3D interpolants of the model fields
    interpolants = modelinterpolants(tt,model)
   
    ## Advect particles backward in time if backward particle tracking
    
    # Skip a time step in the case of backward particle tracking. this
    # is because the forward case uses the velocity field at time t to
    # compute trajectories from t to t+Dt, so backward tracking should
    # use the same velocity field to compute the trajectory from t+Dt
    # to t.
    
    #if particle.direction=='backward' and tt == particle.initime:
    #    None # Skip a time step
    #elif particle.direction=='backward':
    #    None #run advectparticles
    
    ## Seed particle when required
    
    # If day of year matches the seeding frequency (threshold is to account
    # for rouding error is inifreq is not a power of 2).
    if (tt-particle.initime) %particle.inifreq<=1e-10 and particle.ininumber != 0:

        # If seeding is set to 'dynamic', the iniparticle routine uses the
        # current particle position of the 1st particle to seed (at a 
        # constant depth)
        
        #if particle.initype=='dynamic') and tt != particle.initime:
        #    particle.istart = parti.x(parti.id==1) 
        #    particle.irange = 1 
        #    particle.irez = 1 
        #    
        #    particle.jstart = parti.y(parti.id==1) 
        #    particle.jrange = 1 
        #    particle.jrez = 1 
        #    
        #    particle.kstart = -50 
        #    particle.krange = 1 
        #    particle.krez = 1 
        #
        
        # Seed particles
        try:
            parti
        except:
            parti=pd.DataFrame([])
        parti = iniparticles(tt,particle,parti,model)
        particle.ininumber = particle.ininumber -1 
    
    # printlay total number of particle in the experiment
    if model['verbose']:
        print('Total number of particles is %d.' %len(parti.x))
    
    ## Get the model variables interpolated onto the particle's position
    # Interpolate model variables onto particle's position in 4D
    parti = getpartivar(parti,tt,model,interpolants) # this is a vectorized opertion (fast)
    
    ## Save particle position when required
    if (tt-particle.initime)%particle.outfreq<=1e-10:
        saveparti(parti,tt,particle,model)
    
    ## Advect particles forward in time if forward particle tracking (parallel)
    #parti_dask = dd.from_pandas(parti,npartitions=100)
    if model['verbose']:
        print('PROJECT PARTICLE DATA... ')
        
    def apply_advectparticles_to_df(df,timestep,direction,model): 
        return df.apply(advectparticles,axis=1,args=(timestep,direction,model))
    
    if model['parallel']:
        
        parti = dd.from_pandas(parti,npartitions=8).map_partitions(lambda df : df.apply(
                advectparticles,axis=1,args=(thetimestep,thedirection,model)),meta=pd.DataFrame(dtype=float,columns=parti.columns))
        parti = client.persist(parti)
        parti = parti.compute()
        
    else:
        parti = apply_advectparticles_to_df(parti,particle.timestep,particle.direction,model)
    
    if model['verbose']:
        print('DONE')
        print('#################')
        print(' ')
    
    if model['verbose']:
        print('Removed %d particles' % np.count_nonzero(parti.isnull()))
#     if np.count_nonzero(parti.isnull()) == parti.shape[0]:
#         print('END-OF-RUN No particle left to track')
#         parti = parti.dropna()
#         exit(0)
#     else:
#         parti = parti.dropna()
    

    if model['verbose']:
        # Displays time for each time loop and projects the remaining time
        telapsed = (time()-time_loop)/60 # in minutes
        print('This step took %04.2f minutes' %(telapsed) )
        print('At this rate, you have %4.2f minutes left' %(telapsed*(Titer-tcount+1)))
    
print('Entire run took %03.2f minutes' %((time()-time_total)/60)) 