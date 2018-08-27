#!/usr/bin/env python3
# -*- coding: utf-8 -*-
def iniparticles(tt,particle,parti):
    """
    This routine initializes the particles. It is only called when seeding
    particles. It creates the variables, and assigns all intrinsec particles
    characteristics: the DOY, ID, and the vertical sinking velocity assigned
    to the particle. The time-dependant variables at the particles' location 
    are assigned in another routine: "getpartivar.m".
    
    Questions:
        1) Can doy ever be not empy when this is called?
    
    """
    import sys
    import numpy as np
    import pandas as pd
    
    print('SEEDING PARTICLES...')
    
    # Number of particle to seed
    nprtoseed = np.int( np.ceil(particle.irange/particle.irez)*np.ceil(particle.jrange/particle.jrez)*np.ceil(particle.krange/particle.krez)*particle.numofclasses)
    print('Seeding %d particles!' %nprtoseed)
    
    #if tt == particle.initime:
        # If first seed, create the structure for particle data
    seeding=pd.DataFrame(index=np.arange(nprtoseed))
        
    # Particle day of year
    seeding['doy'] = tt*np.ones((nprtoseed,1))
    
    # Particle ID
    blop = (abs((tt-particle.initime)/particle.inifreq*nprtoseed)+\
                   np.arange(1,abs(tt-particle.initime)/particle.inifreq*nprtoseed+nprtoseed+1)).astype(int)
    #blop = length(parti.x)+1:length(parti.x)+1+nprtoseed
    seeding['id'] = blop
    
    xtemp=[]
    ytemp=[]
    ztemp=[]
    wsinktemp=[]
    for theclass in range( particle.numofclasses ): # careful with zero indexing
        # Particle seeding location
        [x,y,z] = np.meshgrid( np.arange(particle.istart,particle.istart+particle.irange,particle.irez),\
            np.arange(particle.jstart,particle.jstart+particle.jrange,particle.jrez),\
            np.arange(particle.kstart,particle.kstart+particle.krange,particle.krez))
        
        xtemp = np.append(xtemp,x)
        ytemp = np.append(ytemp,y)
        ztemp = np.append(ztemp,z)
    
        # Prescribe sinking velocity (in m/s)
        if theclass == 0:
            wsinktemp = np.append(wsinktemp,np.zeros((nprtoseed//particle.numofclasses,1)))
        elif theclass == 1:
            wsinktemp = np.append(wsinktemp,-1/86400*np.ones((nprtoseed//particle.numofclasses,1)))
        elif theclass == 2:
            wsinktemp= np.append(wsinktemp,-10/86400*np.ones((nprtoseed//particle.numofclasses,1)))
        elif theclass == 3:
            wsinktemp = np.append(wsinktemp,-50/86400*np.ones((nprtoseed//particle.numofclasses,1)))
        else:
            print('non-defined sinkinge velocity class') 
            sys.exit(0)
            
    seeding['x'] = xtemp        
    seeding['y'] = ytemp
    seeding['z'] = ztemp

    # Particle velocity at timestep
    # Must be zero when seeded
    seeding['u'] = np.zeros((nprtoseed,1))
    seeding['v'] = np.zeros((nprtoseed,1))
    seeding['w'] = np.zeros((nprtoseed,1))
    seeding['wsink'] = wsinktemp
    seeding['wtotal'] = seeding.w+seeding.wsink
    
    # Get particle characteristics
    # Must be zero when seeded
    seeding['s'] = np.zeros((nprtoseed,1))
    seeding['t'] = np.zeros((nprtoseed,1))
    seeding['rho'] = np.zeros((nprtoseed,1))
    seeding['pv'] = np.zeros((nprtoseed,1))
    seeding['vor'] = np.zeros((nprtoseed,1))
    
    parti = pd.concat([parti,seeding])
    
    
    print('DONE')
    print('#################')
    print(' ')
    return parti