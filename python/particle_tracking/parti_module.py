import dask

def saveparti(parti,tt,particle,model):
    """
    This routine saves the particle data either in a CSV file, or in an SQLite database.
    If saved in a CSV file, one file is created for every recorded timestep.
    If saved as sqlite, data is added to the database at every recorded timestep.
    
    """
    import os
    from sqlalchemy import create_engine
    import pandas as pd
    
    parti2 = parti.copy()
    parti2.columns = ['DOY', 'ID','x','y','z','u','v','w','wsink','wtotal','salinity','temperature','density','PV','vorticity','mld']
    
    # change data typed to reduce filesize
    df_float = parti2.select_dtypes(include=['float'])
    converted_float = df_float.apply(pd.to_numeric,downcast='float')
    parti2[converted_float.columns] = converted_float
    
    df_int = parti2.select_dtypes(include=['int'])
    converted_int = df_int.apply(pd.to_numeric,downcast='unsigned')
    parti2[converted_int.columns] = converted_int
    
    ## CSV file format
    if particle.outputformat=='csv':
        print('SAVE PARTICLE DATA in CSV-file... ')
        
        # Define the variable to print in output file
        # exclude velocity field at previous timestep
        #vartoprint = [1,2,6:16]
        
        # Define filename
        #FILENAME = os.path.join(particle.outputdir,particle.outputfilename+'_doy'+str(tt)+'.csv')
        FILENAME = os.path.join(particle.outputdir,particle.outputfilename+'.csv')
        
        # test if file already exist
        if tt==particle.initime:
            if os.path.isfile(FILENAME):
                input("File exists, press Enter to continue...")
        
        if tt==particle.initime:
            parti2.to_csv(FILENAME, index=False, mode='w')
        else:
            parti2.to_csv(FILENAME, index=False, mode='a',header=False)
            
    
        ## SQLite database
    elif particle.outputformat=='sqlite':
        print('SAVE PARTICLE DATA in SQLite... ')
        
        # Set database filename
        DBFILENAME = os.path.join(particle.outputdir,particle.outputfilename+'_'+particle.direction+'.db')
        
        if tt==particle.initime:
            if os.path.isfile(DBFILENAME):
                input("File exists, press Enter to continue...")
            
        # Set database table name
        tablename = 'particles'
        
        # Creates the SQLite database if does not exist
        engine =  create_engine('sqlite:///'+DBFILENAME, echo=False)
        # Creates the table in the database
        
        parti2.to_sql(tablename, con=engine,if_exists='append',index=False)
    
        # Close connection
        engine.dispose()
        if model['verbose']:
            print('Disonnected to SQLite database')
        
    else:
        raise ValueError('Cannot recognize particle output format.')
    
    if model['verbose']:
        print('DONE')
        print('#################')
        print(' ')
    return

def iniparticles(parti,tt,particle,model):
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
    seeding['doy'] = tt*np.ones((nprtoseed,1))[0] if nprtoseed==1 else tt*np.ones((nprtoseed,1))
    
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
            raise ValueError('Sinking Class not recognized.')
            
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
    
    # Get particle characteristics
    # Must be zero when seeded
    seeding['s'] = np.zeros((nprtoseed,1))
    seeding['t'] = np.zeros((nprtoseed,1))
    seeding['rho'] = np.zeros((nprtoseed,1))
    seeding['pv'] = np.zeros((nprtoseed,1))
    seeding['vor'] = np.zeros((nprtoseed,1))
    
    parti = pd.concat([parti,seeding])
    
    if model['verbose']:
        print('DONE')
        print('#################')
        print(' ')

    return parti

def getpartivar(parti,tt,model,interpolants):
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
    
    ## Hydrography
    parti['t'] = interpolants['FT']((parti.x,parti.y,parti.z)) 
    parti['s'] = interpolants['FS']((parti.x,parti.y,parti.z)) 
    parti['rho'] = interpolants['FRHO']((parti.x,parti.y,parti.z)) 
    
    # mixed layer depth
    parti['mld'] = interpolants['Fmld'](parti.x,parti.y)
        
    ## Potential-vorticity
    parti['pv'] = interpolants['FPV']((parti.x,parti.y,parti.z)) 
        
    ## Relative vorticity
    parti['vor'] = interpolants['Fvor']((parti.x,parti.y,parti.z)) 
    parti['wtotal'] = parti.wsink+parti.w
    
    ##
    if model['verbose']:
        print('DONE')
        print('#################')
        print(' ')

    return parti

def advectparticles(parti,thetimestep,thedirection,model):

    """
    This routine projects the particle forward.
    
    This is modified for parallelization. It moves one particle at a call which is not optimal. Vectorization tricky though.
    
    STEPS
    ---------
           - Get the velocity fluxes at the faces surrounding the particle's
               location
           - Compute the variables necessary for the differential equation
               (Beta and delta)
           - Compute the shortest time it would take for the particle to reach
               one of the cell face (in x, y, or z)
               if the time is shorter than the particle tracking timestep,
               then divide the timestep into subtimesteps
           - Get new particle position
    
    """
    import numpy as np
    
    timestep_leftover = thetimestep*86400
    while timestep_leftover!=0:
        # Find indices of cell faces surrounding the particle
        # position in x
        face_im1 = np.where(model['xf']<=parti.x)[0][-1]
        face_i = face_im1+1
        face_jm1 = np.where(model['yf']<=parti.y)[0][-1]
        face_j = face_jm1+1
        face_km1 = np.where(model['zf']<=parti.z)[0][-1]
        face_k = face_km1+1
        
        # Compute some useful variables
        model_im1 = model['xf'][face_im1]
        model_i = model['xf'][face_i]
        model_jm1 = model['yf'][face_jm1]
        model_j = model['yf'][face_j]
        model_km1 = model['zf'][face_km1]
        model_k = model['zf'][face_k]

        # These are the conditional indexing used for every face velocity.
        # conditional indexing is time consuming so indices are identified
        # separately and once.
        i1 = np.where((model['x'][1:-1]>=model_im1) & \
                      (model['x'][1:-1]<=model_i))[0][0]            
        i2 = np.where((model['y'][1:-1]>=model_jm1) & \
                      (model['y'][1:-1]<=model_j))[0][0]
        i3 = np.where((model['z'][1:-1]>=model_km1) & \
                      (model['z'][1:-1]<=model_k))[0][0]
        
        # Correct indices if particle on a cell face. For instance, if a
        # particle is on a cell face in x, the cell to be considered for
        # advection is the one to the right of the particle if the face
        # velocity is +ve, and to the left of the particle is face velocity
        # is -ve. The cell to the right is considered by default (see
        # above), and this piece of code modifies this if the face velocity
        # is negative.
        
        # define direction of the tracking
        if thedirection=='forward':
            direction = 1
        elif thedirection=='backward':
            direction = -1
        else:
            raise ValueError('Direction must be forward or backward.')
        
        #A
        if direction*(model['uf'][face_im1,i2,i3])<0:
            if parti.x-model_im1 == 0:
            
            # Introduce periodicity at western boundary
                if face_im1 == 0 and model['periodic_ew'] ==1:
                # Select the last cell instead of the first one
                    face_im1 = len(model['xf'])-2
                    face_i = len(model['xf'])-1
                    parti['x'] = model_i
                
                # If not at model boundary
                else:
                    face_im1 = face_im1-1
                    face_i = face_i-1

            # NOTE: Periodicity at eastern boundary is not necessary, as
            # the particle location in x is modular of the model domain
            # (see mod function below). so when eastern most cell face is
            # reached, particle location is changed western most cell
            
        #B
        if direction*model['vf'][i1,face_jm1,i3]<0: 
            if parti.y-model_jm1 == 0:
            
            # If the particle is located at the very 1st cell, indices
            # cannot be switched to the previous cell (southward).
                if face_jm1 == 0:
                # If model is periodic in y, then periodicity is introduced
                    if model['periodic_ns'] ==1:
                    # Select the last cell instead of the first one
                        face_jm1 = len(model['yf'])-2
                        face_j = len(model['yf'])-1
                        parti['y'] = model_j

                # If not at model boundary
                else:
                    face_jm1 = face_jm1-1
                    face_j = face_j-1
            # NOTE: Periodicity at northern boundary is not necessary, as the
            # particle location in y is modular of the model domain (see mod
            # function below). so when northern most cell face is reached, particle
            # location is changed to southern most cell

        #C
        if direction*model['wf'][i1,i2,face_km1]<0:
            if parti.z-model_km1 == 0:
            
                face_km1 = face_km1-1
                face_k = face_k-1
        
        # Update the useful variables
        model_im1 = model['xf'][face_im1]
        model_i = model['xf'][face_i]
        model_jm1 = model['yf'][face_jm1]
        model_j = model['yf'][face_j]
        model_km1 = model['zf'][face_km1]
        model_k = model['zf'][face_k]
        
        i1 = np.where((model['x'][1:-1]>=model_im1) & \
                      (model['x'][1:-1]<=model_i))[0][0]            
        i2 = np.where((model['y'][1:-1]>=model_jm1) & \
                      (model['y'][1:-1]<=model_j))[0][0]
        i3 = np.where((model['z'][1:-1]>=model_km1) & \
                      (model['z'][1:-1]<=model_k))[0][0]
        
        # Compute Cell spacings in x, y, and z
        # Compute Dx,Dy,Dz
        Dx = model_i-model_im1
        Dy = model_j-model_jm1
        Dz = model_k-model_km1
        
        # Compute some variables in x
        rx0 = parti.x/Dx
        rxim1 = model_im1/Dx
        rxi = model_i/Dx
        betax = model['uf'][face_im1,i2,i3] - model['uf'][face_i,i2,i3]
        deltax = -model['uf'][face_im1,i2,i3]- betax*rxim1
        
        # Compute some variables in y
        ry0 = parti.y/Dy
        ryjm1 = model_jm1/Dy
        ryj = model_j/Dy
        betay = model['vf'][i1,face_jm1,i3] - model['vf'][i1,face_j,i3]
        deltay = -model['vf'][i1,face_jm1,i3]-betay*ryjm1
        
        # Compute some variables in z
        rz0 = parti.z/Dz
        rzkm1 = model_km1/Dz
        rzk = model_k/Dz
        betaz = model['wf'][i1,i2,face_km1] - model['wf'][i1,i2,face_k]
        deltaz = -model['wf'][i1,i2,face_km1] - betaz*rzkm1
        
        ## Compute the shortest time it would take the particle to reach one of the
        # cell faces
        Dtmaxtemp = np.zeros(6,dtype='complex64')

        # If no velocity gradient (i.e. beta == 0), the solution to the
        # differential equation changes, hence the if-loop
        # note that all the invuf, infvf, infwf can be infinite!
        if betax == 0:
            invuf = 1/model['uf'][face_im1,i2,i3]
            # time to reach i-1 th face            
            Dtmaxtemp[0] = (rxim1-rx0)*invuf*Dx*Dy*Dz
            # time to reach i th face
            Dtmaxtemp[1] = (rxi-rx0)*invuf*Dx*Dy*Dz
        else:
            invbetax = 1/betax
            # time to reach i-1 th face
            Dtmaxtemp[0] = -invbetax*np.log((rxim1 + deltax*invbetax)/(rx0 + deltax*invbetax))*Dx*Dy*Dz
            # time to reach i th face
            Dtmaxtemp[1] = -invbetax*np.log((rxi + deltax*invbetax)/(rx0 + deltax*invbetax))*Dx*Dy*Dz
            
        if betay == 0:
            invvf = 1/model['vf'][i1,face_jm1,i3]
            # time to reach j-1 th face
            Dtmaxtemp[2] = (ryjm1-ry0)*invvf*Dx*Dy*Dz
            # time to reach j th face
            Dtmaxtemp[3] = (ryj-ry0)*invvf*Dx*Dy*Dz
        else:
            # time to reach j-1 th face
            invbetay = 1/betay
            Dtmaxtemp[2] = -invbetay*np.log((ryjm1 + deltay*invbetay)/(ry0 + deltay*invbetay))*Dx*Dy*Dz
            # time to reach j th face
            Dtmaxtemp[3] = -invbetay*np.log((ryj + deltay*invbetay)/(ry0 + deltay*invbetay))*Dx*Dy*Dz

        if betaz == 0:
            invwf = 1/model['wf'][i1,i2,face_km1]
            # time to reach k-1 th face
            Dtmaxtemp[4] = (rzkm1-rz0)*invwf*Dx*Dy*Dz
            # time to reach j th face
            Dtmaxtemp[5] = (rzk-rz0)*invwf*Dx*Dy*Dz
        else:
            invbetaz = 1/betaz
            # time to reach k-1 th face
            Dtmaxtemp[4] = -invbetaz*np.log((rzkm1 + deltaz*invbetaz)/(rz0 + deltaz*invbetaz))*Dx*Dy*Dz
            # time to reach j th face
            Dtmaxtemp[5] = -invbetaz*np.log((rzk + deltaz*invbetaz)/(rz0 + deltaz*invbetaz))*Dx*Dy*Dz
                    
        # Imaginary times occur when the face considered has a velocity in
        # the opposite direction of the velocity at the particle's
        # location. Hence the imaginary number: does not matter how long,
        # the particle will never reach that face.
#%           
        if thedirection=='forward':
            Dtmax = Dtmaxtemp[ (Dtmaxtemp>0) & (Dtmaxtemp.imag==0)].min().real
            # If the shortest time it would take the particle to reach a cell face
            # is shorter than the particle tracking timestep, then the integration
            # timestep needs to be shortened
            if not Dtmax.any():
                intermediate_timestep = timestep_leftover
            else:
                intermediate_timestep = np.min([Dtmax,timestep_leftover])
                
            timestep_leftover = timestep_leftover-intermediate_timestep
            ds = intermediate_timestep/(Dx*Dy*Dz)
            
            
        elif thedirection=='backward':
            Dtmax = Dtmaxtemp[ (Dtmaxtemp>0) & (Dtmaxtemp.imag==0)].min().real
            
            # If the shortest time it would take the particle to reach a cell face
            # is shorter than the particle tracking timestep, then the integration
            # timestep needs to be shortened
            if not Dtmax.any():
                intermediate_timestep = timestep_leftover
            else:
                intermediate_timestep = np.max([Dtmax,timestep_leftover])
            
            timestep_leftover = timestep_leftover-intermediate_timestep
            ds = intermediate_timestep/(Dx*Dy*Dz)
            
        #=================
        #=== Assign x-position to particle.
        #=================
        # If no velocity gradient (i.e. beta == 0), the solution to the
        # differential equation changes, hence the if-loop
        if betax == 0:
            rx1 = rx0 + model['uf'][face_im1,i2,i3]*ds;
        else:
            rx1 = (rx0 + deltax*invbetax)*np.exp(-betax*ds)-(deltax*invbetax)

        if abs(rx1-rxim1)<1e-11:
            parti['x'] = model_im1
        elif abs(rx1-rxi)<1e-11:
            parti['x'] = model_i
        else:
            parti['x'] = rx1*Dx
        
        #=================
        #=== Assign y-position to particle.
        #=================
        # If no velocity gradient (i.e. beta == 0), the solution to the
        # differential equation changes, hence the if-loop
        if betay == 0:
            ry1 = ry0 + model['vf'][i1,face_jm1,i3]*ds;
        else:
            ry1 = (ry0 + deltay*invbetay)*np.exp(-betay*ds)-(deltay*invbetay)

        if abs(ry1-ryjm1)<1e-11:
            parti['y'] = model_jm1
        elif abs(ry1-ryj)<1e-11:
            parti['y'] = model_j
        else:
            parti['y'] = ry1*Dy
        #=================
        #=== Assign z-position to particle.
        #=================
        # If no velocity gradient (i.e. beta == 0), the solution to the
        # differential equation changes, hence the if-loop
        if betaz == 0:
            rz1 = rz0 + model['wf'][i1,i2,face_km1]*ds;
        else:
            rz1 = (rz0 + deltaz*invbetaz)*np.exp(-betaz*ds)-(deltaz*invbetaz)

        if abs(rz1-rzkm1)<1e-11:
            parti['z'] = model_km1
        elif abs(rz1-rzk)<1e-11:
            parti['z'] = model_k
        else:
            parti['z'] = rz1*Dz + parti['wsink']*intermediate_timestep
        
        if model['periodic_ew'] == 1:
            # Wrap the domain around if periodicity
            parti['x'] = parti.x % model['xf'][-1] #maximum of model.xf assumes xf to be sorted
        
        if model['periodic_ns'] == 1:
            # Wrap the domain around if periodicity
            parti['y'] = parti.y % model['yf'][-1] #maximum of model.xf assumes xf to be sorted
        
        # Particles can't go airborn (z>0), and stop tracking if too deep
        if parti.z>0:
            parti['z'] = 0
        elif parti.z< model['zf'][0]: #minimum of zf
            parti['z'] = np.nan
            timestep_leftover=0
        
        # If hits a solid boundary, kill the particle
        if parti.z>=model['yf'][-1] or parti.y<=model['yf'][0]:
            parti['y'] = np.nan
            timestep_leftover=0


    return parti