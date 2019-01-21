def initialize_test():
    ''' The function sets all the important parameters for a test.
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
    particle = {} # create dict object

    # Enable depth-keeping particles
    particle['isobaric'] = True
    
    # day of year for particle initialization (i.e., first seeding)
    particle['initime'] = 10
    # frequency of particles seeding (in days)
    particle['inifreq'] = 2**-2 
    # Number of seeding events
    particle['ininumber'] = 1
    # Type of particle seeding ('dynamic' or 'static')
    particle['initype'] = 'static' 

    # Number of particle classes to be released. Please refer to iniparticles
    # to customize the sinking velocity of each particle class
    particle['numofclasses'] = 1

    # Parameters to customize the initialization of the particles (in meters)
    # Only applied if seeding type is 'static'
    particle['istart'] = 48000
    particle['irange'] = 4000
    particle['irez'] = 1000

    particle['jstart'] = 48000
    particle['jrange'] = 4000
    particle['jrez'] = 1000

    particle['kstart'] = -100
    particle['krange'] = 100
    particle['krez'] = 1

    #----------------
    # ADVECTION
    #----------------
    # timestep of particle tracking in days
    particle['timestep'] = 2**-4
    # length of particle tracking experiment in days
    particle['length'] = 1
    # Define particle-tracking direction ('forward' or 'backward')
    particle['direction'] = 'forward' 

    #----------------
    # OUTPUT
    #----------------
    # frequency of particles output (in days)
    particle['outfreq'] = 2**-3
    # Format of output ('csv' or 'sqlite')
    #particle.outputformat = 'csv' 

    #particle.outputdir = './committee_meeting_3d/' 
    #particle.outputfilename = 'vertical_parti%d_dk' %np.int(sys.argv[1])

    # if os.path.exists(os.path.join(particle.outputdir,particle.outputfilename)):
    #    sys.exit('Warning file already exists!')


    #####################
    ## CORE CODE
    #####################
    time_total = datetime.datetime.now()

    # Change a few things if backward-tracking
    if particle['direction']=='backward':
        particle['timestep'] = -particle['timestep']
        particle['length'] = -particle['length']
        
    return model, particle

def getpartivar_test(parti,tt,model,interpolants):
    """
    This routine interpolate the velocty field in time and space onto the
    particles location using a 4-D linear interpolant. The velocity field has
    been interpolated linearly in time already in "getvelocityfield.m", so
    only the hydrographic variables are interpolated in time as well as in
    space.
    
    note:
    the above is not true anymore. i interpolate every field in time in getvelocityfield
    
    Parameter
    ------------
    tt: int
        Day of year
    parti: pandas DataFrame
        parti is being updated with the interplated values of the variables at the particles' location.
    """
    
    if model['verbose']:
        print('INTERPOLATE VARIABLES ONTO PARTICLES... ')
    parti['doy'] = tt

    parti['u'] = interpolants['Fu']((parti.x,parti.y,parti.z))
    parti['v'] = interpolants['Fv']((parti.x,parti.y,parti.z)) 
    parti['w'] = interpolants['Fw']((parti.x,parti.y,parti.z))

    ## Relative vorticity
    parti['wtotal'] = parti.wsink+parti.w
    
    ##
    if model['verbose']:
        print('DONE')
        print('#################')
        print(' ')

    return parti

def iniparticles_test(tt,particle,parti,model):
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
    nprtoseed = np.int( np.ceil(particle['irange']/particle['irez'])*\
                       np.ceil(particle['jrange']/particle['jrez'])*\
                       np.ceil(particle['krange']/particle['krez'])*particle['numofclasses'])
    print('Seeding %d particles!' %nprtoseed)
    
    #if tt == particle.initime:
        # If first seed, create the structure for particle data
    seeding=pd.DataFrame(index=np.arange(nprtoseed))
        
    # Particle day of year
    if nprtoseed==1:
        seeding['doy'] = tt*np.ones((nprtoseed,1))[0]
    else:
        seeding['doy'] = tt*np.ones((nprtoseed,1))
    
    # Particle ID
    blop = (abs((tt-particle['initime'])/particle['inifreq']*nprtoseed)+\
                   np.arange(1,abs(tt-particle['initime'])/particle['inifreq']*nprtoseed+nprtoseed+1)).astype(int)
    #blop = length(parti.x)+1:length(parti.x)+1+nprtoseed
    seeding['id'] = blop
    
    xtemp=[]
    ytemp=[]
    ztemp=[]
    wsinktemp=[]
    for theclass in range( particle['numofclasses'] ): # careful with zero indexing
        # Particle seeding location
        [x,y,z] = np.meshgrid( np.arange(particle['istart'],particle['istart']+particle['irange'],particle['irez']),\
            np.arange(particle['jstart'],particle['jstart']+particle['jrange'],particle['jrez']),\
            np.arange(particle['kstart'],particle['kstart']+particle['krange'],particle['krez']))
        
        xtemp = np.append(xtemp,x)
        ytemp = np.append(ytemp,y)
        ztemp = np.append(ztemp,z)
    
        # Prescribe sinking velocity (in m/s)
        if theclass == 0:
            wsinktemp = np.append(wsinktemp,np.zeros((nprtoseed//particle['numofclasses'],1)))
        elif theclass == 1:
            wsinktemp = np.append(wsinktemp,-1/86400*np.ones((nprtoseed//particle['numofclasses'],1)))
        elif theclass == 2:
            wsinktemp= np.append(wsinktemp,-10/86400*np.ones((nprtoseed//particle['numofclasses'],1)))
        elif theclass == 3:
            wsinktemp = np.append(wsinktemp,-50/86400*np.ones((nprtoseed//particle['numofclasses'],1)))
        else:
            print('blop') 
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
    seeding['wtotal'] = (seeding.w+seeding.wsink)*0
    
    parti = pd.concat([parti,seeding])
    
    if model['verbose']:
        print('DONE')
        print('#################')
        print(' ')

    return parti

def getvelocityfield_test(model,t,mode='uniform'):
    """
    This routine makes up a velocity field for testing purposes.
    
    Parameters
    ----------
    model : xarray Dataset
        Contains model fields and meta data    
    tt : int
        Day of year at which velocity field is needed.   
    
    Returns
    ----------
    model : xarray Dataset 
        Updated structure now including model velocity.
    """
    
    import sys
    from glob import glob
    import numpy as np
    import xarray as xr
    import pandas as pd
    
    if model['verbose']:
        print('MAKING TEST VELOCITY FIELDS...')
       
    dx = 1000
    dy = 1000
    dz = 20
    model['x'] = np.arange(0,100000,dx)
    model['y'] = np.arange(0,100000,dy)
    model['z'] = np.arange(-1000,1,dz)
    
    model['xf'] = model['x'][:-1]+ 0.5*dx
    model['yf'] = model['y'][:-1]+ 0.5*dy  
    model['zf'] = model['z'][:-1]+ 0.5*dy  
    
    if mode == 'random':
        model['uf'] = np.ones((len(model['xf']),len(model['yf']),len(model['zf']) ))
        model['vf'] = np.ones((len(model['xf']),len(model['yf']),len(model['zf']) ))
        model['wf'] = np.ones((len(model['xf']),len(model['yf']),len(model['zf']) ))

        model['u'] = np.ones((len(model['x']),len(model['y']),len(model['z']) ))
        model['v'] = np.ones((len(model['x']),len(model['y']),len(model['z']) ))
        model['w'] = np.ones((len(model['x']),len(model['y']),len(model['z']) ))
        
    elif mode == 'uniform':
        model['uf'] = np.ones((len(model['xf']),len(model['yf']),len(model['zf']) ))
        model['vf'] = np.ones((len(model['xf']),len(model['yf']),len(model['zf']) ))
        model['wf'] = np.ones((len(model['xf']),len(model['yf']),len(model['zf']) ))

        model['u'] = np.ones((len(model['x']),len(model['y']),len(model['z']) ))
        model['v'] = np.ones((len(model['x']),len(model['y']),len(model['z']) ))
        model['w'] = np.ones((len(model['x']),len(model['y']),len(model['z']) ))
    
    if model['verbose']:
        print('DONE')
        print('#################')
        print(' ')

    return model

def modelinterpolants_test(tt,model):
    '''
    Take model fields and create 3D interpolants
    These will be used in "getpartivar" to interpolate model
    data on particle positions. 
    
    previously in getpartivar but removed for parallelization.
    
    '''
    from scipy.interpolate import RegularGridInterpolator,interp2d
    import xarray as xr
    import numpy as np
    
    
    # velocity
    Fu = RegularGridInterpolator((model['x'],model['y'],model['z']),model['u'])
    Fv = RegularGridInterpolator((model['x'],model['y'],model['z']),model['v'])
    Fw = RegularGridInterpolator((model['x'],model['y'],model['z']),model['w'])
    
    
    interpolants={'Fu':Fu,'Fv':Fv,'Fw':Fw}
    
    return interpolants