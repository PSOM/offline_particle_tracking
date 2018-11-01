def getvelocityfield(model,tt):
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
    import sys,glob
    import numpy as np
    import xarray as xr
    import pandas as pd
    
    print('IMPORTING VELOCITY FIELDS...')
    
    # List of model outputs present in the folder path
    files = np.sort( glob.glob(model['path']+'face_*.cdf'))
    files2 = np.sort( glob.glob(model['path']+'full_*.cdf'))
    
    # Delete some outputfile filenames to match the desired model output
    # temporal resolution to be used in the track
    
    files = files[1:-1:model['outputskip']]
    files2 = files2[1:-1:model['outputskip']]
    
    # Rename the considered day of year
    doy = tt
    
    # find the two files to be imported
    timestepind = []
    for ii in range(len(files)-1):
        tmin = np.int(files[ii].split('_')[-1].split('.')[0])/(86400/model['timestep']) + model['start_day']
        tmax = np.int(files[ii+1].split('_')[-1].split('.')[0])/(86400/model['timestep']) + model['start_day']
        if tmin <= doy and tmax >= doy:
            timestepind = [ii, ii+1]
            timestepdoy = [tmin,tmax]
            break     
    del ii
    
    if not timestepind:
        print('No output files were found')
        sys.exit(0)
    
    # reads all variables in face_*
    ds_face = xr.open_mfdataset(files[timestepind],concat_dim='doy',chunks=None)
    # scale stuff to SI units
    #ds_face['uf'] =  ds_face['uf']*1e3*1e5   
    #ds_face['vf'] =  ds_face['vf']*1e3*1e5   
    #ds_face['wf'] =  ds_face['wf']*1e3*1e5   
    
    # reads all variables in full_*
    ds_full = xr.open_mfdataset(files2[timestepind],concat_dim='doy',chunks=None)
    # scale stuff to SI units
    model['x'] = ds_full['xc'][0,:].values*1e3
    model['y'] = ds_full['yc'][0,:].values*1e3
    #ds_full['w'] = ds_full['w']*1e-3
    ds_full['z'] = ds_full['zc'] # watch out maybe needs to be scaled!
    model['z'] = ds_full['z'][0,:,0,0].squeeze().values
    
    # Interpolate linearly in time the velocity field
    def interp_doy(dataarray,timestepdoy=timestepdoy,tt=tt):
        return dataarray.sel(doy=0)+((tt-timestepdoy[0])/(timestepdoy[1]-timestepdoy[0])*(dataarray.sel(doy=1)-dataarray.sel(doy=0)))
    
    model['uf'] = interp_doy(ds_face.uf).T.values
    model['vf']= interp_doy(ds_face.vf).T.values
    model['wf']= interp_doy(ds_face.wf).T.values
    
    # define coordinates (not in PSOM output files yet...)
    model['xf'] = model['x'][:-1]+ 0.5*np.diff(model['x'])
    model['yf'] = model['y'][:-1]+ 0.5*np.diff(model['y'])
    
    # Import zf from zgrid.out
    zgrid = pd.read_csv(model['path']+'zgrid.out',skipinitialspace=True,sep=' ',header=None)
    zgrid=zgrid[1].values
    model['zf'] = ( zgrid[ np.int( np.where( zgrid=='face' )[0])+2: -1] ).astype(np.float)
    
    model['u'] = interp_doy(ds_full.u).T.values
    model['v'] = interp_doy(ds_full.v).T.values
    model['w'] = interp_doy(ds_full.w).T.values
    model['T'] = interp_doy(ds_full.temp).T.values
    model['S'] = interp_doy(ds_full.s).T.values
    model['RHO'] = interp_doy(ds_full.rho).T.values
    model['PV']  = interp_doy(ds_full.pv).T.values
    model['vor']  = interp_doy(ds_full.vor).T.values
    
    # Grid the 3D coordinates
    #[model['x2'],model['y2'],model['z2']] = np.meshgrid(ds_full.x,ds_full.y,ds_full.z)
    #[model['ufx'],model['ufy'],model['ufz']] = np.meshgrid(model['xf'],ds_full.y[1:-1],ds_full.z[1:-1])
    #[model['vfx'],model['vfy'],model['vfz']] = np.meshgrid(ds_full.x[1:-1],model['yf'],ds_full.z[1:-1])
    #[model['wfx'],model['wfy'],model['wfz']] = np.meshgrid(ds_full.x[1:-1],ds_full.y[1:-1],model['zf'])
    
    #model['x2'] = ds_full.x
    #model['y2'] = ds_full.y
    #model['z2'] = model['z']
    
    print('DONE')
    print('#################')
    print(' ')

    return model