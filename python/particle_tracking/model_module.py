def saveconfiguration(model,particle):
    '''
    This routine prints the particle-tracking experiment configuration to a text file. 
    The filename includes the date and time of the simulation.
    
    '''
    if model['verbose']:
        print('SAVE EXPERIMENT CONFIGURATION... ')
    import datetime

    # Get current time with datetime
    nowtime = datetime.datetime.now().strftime(' %d%m%y-%H%M')

    # Open the file
    file = open('%s/%s_%s_configuration.out' %(particle['outputdir'],nowtime,particle['outputfilename']),'w')

    # Write the data
    file.write('******** OCEAN MODEL PARAMETERS ********\n')
    file.write('Start day =  %f \n' %model['start_day'])
    file.write('Ocean model timestep =  %f seconds \n' %model['timestep'])
    file.write('Directory of ocean model outputs = %s \n' %model['path'])
    file.write('Ocean model resolution used = every  %f file(s) \n' %model['outputskip'])
    file.write('Ocean model periodicities: E-W =  %f  N-S =  %f \n' %(model['periodic_ew'],model['periodic_ns']))
    file.write('\n')
    file.write('\n')
    file.write('******** PARTICLE TRACKING PARAMETERS ********\n')
    file.write('\n')
    file.write('#----------------\n')
    file.write('SEEDING\n')
    file.write('#----------------\n')
    file.write('Particle initialization (1st seeding) = Day  %f \n' %particle['initime'])
    file.write('Particle seeding frequency = every  %f day(s) \n' %particle['inifreq'])
    file.write('Number of seeding events =  %f \n' %particle['ininumber'])
    file.write('Type of seeding (static vs dynamic) = %s \n' %particle['initype'])
    file.write('\n')
    file.write('Number of particle classes =  %f \n' %particle['numofclasses'])
    file.write('\n')
    file.write('Particle seeding strategy in x:\n')
    file.write('istart =  %f \n' %particle['istart'])
    file.write('irange =  %f \n' %particle['irange'])
    file.write('irez =  %f \n' %particle['irez'])
    file.write('Particle seeding strategy in y:\n')
    file.write('jstart =  %f \n' %particle['jstart'])
    file.write('jrange =  %f \n' %particle['jrange'])
    file.write('jrez =  %f \n' %particle['jrez'])
    file.write('Particle seeding strategy in z:\n')
    file.write('kstart =  %f \n' %particle['kstart'])
    file.write('krange =  %f \n' %particle['krange'])
    file.write('krez =  %f \n' %particle['krez'])
    file.write('\n')
    file.write('#----------------\n')
    file.write('ADVECTION\n')
    file.write('#----------------\n')
    file.write('Timestep for particle tracking =  %f day(s)\n' %particle['timestep'])
    file.write('Length of particle tracking simulation =  %f day(s)\n' %particle['length'])
    file.write('Particle tracking direction: %s \n' %particle['direction'])
    file.write('\n')
    file.write('#----------------\n')
    file.write('OUTPUT\n')
    file.write('#----------------\n')
    file.write('Particle output frequency =  %f day(s)\n' %particle['outfreq'])
    file.write('Particle output format = %s \n' %particle['outputformat'])
    file.write('Particle output directory =%s \n' %particle['outputdir'])
    file.write('Particle output filename = %s \n' %particle['outputfilename'])

    # Close the file
    file.close()

    if model['verbose']:
        print('DONE\n')
        print('#################\n')
        print('\n')
    return

def getvelocityfield(model,tt,particle):
    """
    This routine identifies the model output files surrounding the particle
    tracking timestep to extract the velocity fields. It then linearly
    interpolates in time the velocity field to get the velocity field for the
    timestamp considered
    
    
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
    from scipy.io import netcdf_file
    
    if model['verbose']:
        print('IMPORTING VELOCITY FIELDS...')
    
    # List of model outputs present in the folder path
    face_files = np.sort( glob(model['path']+'face_*.cdf'))
    full_files = np.sort( glob(model['path']+'full_*.cdf'))
    
    # Delete some outputfile filenames to match the desired model output
    # temporal resolution to be used in the track
    face_files = face_files[1:-1:model['outputskip']]
    full_files = full_files[1:-1:model['outputskip']]
    
    # Rename the considered day of year
    doy = tt
    
    # find the two files to be imported
    timestepind = []
    for ii in range(len(face_files)-1):
        tmin = np.int(face_files[ii].split('_')[-1].split('.')[0])/(86400/model['timestep']) + model['start_day']
        tmax = np.int(face_files[ii+1].split('_')[-1].split('.')[0])/(86400/model['timestep']) + model['start_day']
        if tmin <= doy and tmax >= doy:
            timestepind = [ii, ii+1]
            timestepdoy = [tmin,tmax]
            break     
    assert len(timestepind)>0, 'No output files were found'
    
    def read_face_files(files,model):
        ''' Reads FACE model output files and puts them into numpy array,
        these are stored in the dictionary "model" 
        
        !CAREFUL check units and scaling factors of model output, pre-Mat model output NOT SI units CAREFUL!
        '''
        uf=[]
        vf=[]
        wf=[]
        for file in files:
            with netcdf_file(file,'r',) as f: # scipy
                uf.append( f.variables['uf'][:].astype('float')*1e3*1e5 )
                vf.append( f.variables['vf'][:].astype('float')*1e3*1e5 )
                wf.append( f.variables['wf'][:].astype('float')*1e3*1e5 )   
                
        model['uf'] = np.array(uf)
        model['vf'] = np.array(vf)
        model['wf'] = np.array(wf)
        return model
    
    def read_full_files(files,model):
        ''' Reads FULL model output files and puts them into numpy array,
        these are stored in the dictionary "model" 
        
        !CAREFUL check units and scaling factors of model output, pre-Mat model output NOT SI units CAREFUL!
        '''
        x=[];y=[];z=[]
        u=[];v=[];w=[]
        t=[];s=[];rho=[]
        pv=[];vor=[]
        for file in files:
            with netcdf_file(file,'r',) as f: # scipy
                # careful with scaling factors!!
                x.append( f.variables['xc'][:].astype('float')*1e3 )
                y.append( f.variables['yc'][:].astype('float')*1e3 )
                z.append( f.variables['zc'][:].astype('float')*1e3 )
                
                u.append( f.variables['u'][:].astype('float'))
                v.append( f.variables['v'][:].astype('float') )
                w.append( f.variables['w'][:].astype('float')*1e-3 )
                
                t.append( f.variables['temp'][:].astype('float') )
                s.append( f.variables['s'][:].astype('float') )
                rho.append( f.variables['rho'][:].astype('float') )
                pv.append( f.variables['pv'][:].astype('float') )
                vor.append( f.variables['vor'][:].astype('float') )
                
        model['x'] = np.array(x)[0,:]
        model['y'] = np.array(y)[0,:]
        model['z'] = np.array(z)[0,:,0,0].squeeze()
        model['u'] = np.array(u); model['v'] = np.array(v); 
        
        model['w'] = np.array(w)
        model['T'] = np.array(t); model['S'] = np.array(s); model['RHO'] = np.array(rho)
        model['PV'] = np.array(pv); model['vor'] = np.array(vor)
        
        return model

    # reads all variables in face_*
    model = read_face_files(face_files[timestepind],model)
     # reads all variables in face_*
    model = read_full_files(full_files[timestepind],model)
    
    def interp_doy(dataarray,timestepdoy=timestepdoy,tt=tt):
        '''Interpolate linearly in time the velocity field'''
        weight = (tt-timestepdoy[0])/(timestepdoy[1]-timestepdoy[0])
        return dataarray[0]+weight*(dataarray[1]-dataarray[0])
    
    def interpolate_between_doy(model):
        ''' Interpolate all variable in list'''
        varlist = ['u','v','w','uf','vf','wf','T','S','RHO','PV','vor']
        
        for var in varlist:
            model[var] = interp_doy(model[var]).T
        return model
    
    # interpolate linearly in time
    model = interpolate_between_doy(model)
    
    if particle['isobaric']:
        model['wf']=0*model['wf']
    
    # define coordinates (not in PSOM output files yet...)
    # if delta x and y are constant:
    if model['delta_x']>0 and model['delta_y']>0:
        model['xf'] = model['x'][:-1]+ 0.5*model['delta_x'] 
        model['yf'] = model['y'][:-1]+ 0.5*model['delta_y']  
    else:
        model['xf'] = model['x'][:-1]+ 0.5*np.diff(model['x'])
        model['yf'] = model['y'][:-1]+ 0.5*np.diff(model['y'])
    
    def read_zgrid(path):
        '''Read zgrid from model grid and store in "model" dictionary'''
        try:
            zgrid = pd.read_csv(model['path']+'zgrid.out',skipinitialspace=True,sep=' ',header=None)
            zgrid=zgrid[1].values
        except IOError:
            print('File not found.')
        return ( zgrid[ np.int( np.where( zgrid=='face' )[0])+2: -1] ).astype('float')
    
    # Import zf from zgrid.out
    model['zf'] = read_zgrid(model['path'])
    
    if model['verbose']:
        print('DONE')
        print('#################')
        print(' ')

    return model

def modelinterpolants(tt,model):
    '''
    Take model fields and create 3D interpolants
    These will be used in "getpartivar" to interpolate model
    data on particle positions. 
    
    previously in getpartivar but removed for parallelization.
    
    '''
    from scipy.interpolate import RegularGridInterpolator,interp2d
    import xarray as xr
    import numpy as np
    
    def interpolate_mld(model):
        ''' Build 2d interpolant for the mixed layer depth for later use in "getpartivar"
        
        criterion: 0.03kg/m3 difference from 10m density.
        
        finds depths of cells below and above (rho_10+0.03) and interpolates linearly between them. that's faster than 
        doing it for each profile.
        '''
        foo = xr.DataArray(model['RHO'], coords=[model['x'], model['y'], model['z']], dims=['x', 'y', 'z'])
        dss = xr.Dataset({'rho':foo})

        # find density at 10m
        a = dss.rho.interp(x=dss.x,y=dss.y,z=-10) 
        
        # find density at 10m + 0.03
        avalues = np.repeat(a.values[:,:,np.newaxis]+0.03,50,axis=2)
        
        # prepare arrays
        zz,yy,xx = np.meshgrid(dss.z,dss.y,dss.x,indexing='ij')   
        
        # find the two rho values large and smaller than rho10+0.03
        condition_lt = dss.rho < avalues
        condition_gt = dss.rho > avalues
        rho1 = np.max(np.ma.masked_where(condition_lt,dss.rho), axis=2)
        rho2 = np.min(np.ma.masked_where(condition_gt ,dss.rho), axis=2)
        z1 = np.max(np.ma.masked_where(condition_lt,zz.T), axis=2)
        z2 = np.min(np.ma.masked_where(condition_gt ,zz.T), axis=2)
        
        # find mixed layer depth surface
        w1=(rho1- a.values+0.03)/(rho1-rho2) 
        w2= 1-w1
        mld = z1*w2 + z2*w1

        # make interpolant
        f = interp2d(dss.x,dss.y,mld.T)

        return f
    
    # mixed layer depth
    Fmld = interpolate_mld(model)
    
    # velocity
    Fu = RegularGridInterpolator((model['x'],model['y'],model['z']),model['u'])
    Fv = RegularGridInterpolator((model['x'],model['y'],model['z']),model['v'])
    Fw = RegularGridInterpolator((model['x'],model['y'],model['z']),model['w'])
    
    ## Hydrography
    FT = RegularGridInterpolator((model['x'],model['y'],model['z']),model['T'])
    FS = RegularGridInterpolator((model['x'],model['y'],model['z']),model['S'])
    FRHO = RegularGridInterpolator((model['x'],model['y'],model['z']),model['RHO'])

    ## Potential-vorticity
    FPV = RegularGridInterpolator((model['x'],model['y'],model['z']),model['PV'])
        
    ## Relative vorticity
    Fvor = RegularGridInterpolator((model['x'],model['y'],model['z']),model['vor'])
    
    interpolants={'Fu':Fu,'Fv':Fv,'Fw':Fw,'FT':FT,'FS':FS,'FRHO':FRHO,
                  'FPV':FPV,'Fvor':Fvor,'Fmld':Fmld}
    
    return interpolants