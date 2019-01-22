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
from particle_tracking.parti_module import advectparticles
from particle_tracking.test_module import *

def test_advection():
    # set parameters in model and particle
    model, particle = initialize_test()

    # Record the total number of advective timesteps
    Titer = np.arange(particle['initime'],particle['initime']+particle['length']+particle['timestep'],particle['timestep']).size

    # Main Loop
    iteration = [particle['initime']+t*particle['timestep'] for t in range(int( particle['length']/particle['timestep'] ))]
    for tcount,tt in enumerate(iteration):
        
        ## Get the velocity field relevant to the time step considered, and 
        # already interpolated onto the time step
        model = getvelocityfield_test(model,tt)
        # create 3D interpolants of the model fields
        interpolants = modelinterpolants_test(tt,model)

        # If day of year matches the seeding frequency (threshold is to account
        # for rouding error is inifreq is not a power of 2).
        if (tt-particle['initime']) %particle['inifreq']<=1e-10 and particle['ininumber'] != 0:

            # Seed particles
            try:
                parti
            except:
                parti=pd.DataFrame([])
            parti = iniparticles_test(tt,particle,parti,model)
            particle['ininumber'] = particle['ininumber'] -1 

        ## Get the model variables interpolated onto the particle's position
        # Interpolate model variables onto particle's position in 4D
        parti = getpartivar_test(parti,tt,model,interpolants) # this is a vectorized opertion (fast)

        def apply_advectparticles_to_df(df,timestep,direction,model): 
            return df.apply(advectparticles,axis=1,args=(timestep,direction,model))

        parti = apply_advectparticles_to_df(parti,particle['timestep'],particle['direction'],model)

pwd
