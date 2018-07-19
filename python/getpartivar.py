#!/usr/bin/env python3
# -*- coding: utf-8 -*-

def getpartivar(tt,parti,model):
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
    from scipy.interpolate import RegularGridInterpolator
    
    print('INTERPOLATE VARIABLES ONTO PARTICLES... ')
    
    parti['doy'] = tt
    
    ## Velocity field
    
    # If velocity field extracted from face velocity-fluxes (divergence free)
    
    # F = griddedInterpolant(model.ufx,model.ufy,model.ufz,model.uf)
    # parti.u = F(parti.x,parti.y,parti.z)
    # F = griddedInterpolant(model.vfx,model.vfy,model.vfz,model.vf)
    # parti.v = F(parti.x,parti.y,parti.z)
    # F = griddedInterpolant(model.wfx,model.wfy,model.wfz,model.wf)
    # parti.w = F(parti.x,parti.y,parti.z)
    
    # If velocity field extracted from centered velocities (NOT divergence
    # free)
    F = RegularGridInterpolator((model['x'],model['y'],model['z']),model['u'])
    parti['u'] = F((parti.x,parti.y,parti.z))
    F = RegularGridInterpolator((model['x'],model['y'],model['z']),model['v'])
    parti['v'] = F((parti.x,parti.y,parti.z)) 
    F = RegularGridInterpolator((model['x'],model['y'],model['z']),model['w'])
    parti['w'] = F((parti.x,parti.y,parti.z)) 
    
    ## Hydrography
    F = RegularGridInterpolator((model['x'],model['y'],model['z']),model['T'])
    parti['t'] = F((parti.x,parti.y,parti.z)) 
    F = RegularGridInterpolator((model['x'],model['y'],model['z']),model['S'])
    parti['s'] = F((parti.x,parti.y,parti.z)) 
    F = RegularGridInterpolator((model['x'],model['y'],model['z']),model['RHO'])
    parti['rho'] = F((parti.x,parti.y,parti.z)) 
        
    ## Potential-vorticity
    F = RegularGridInterpolator((model['x'],model['y'],model['z']),model['PV'])
    parti['pv'] = F((parti.x,parti.y,parti.z)) 
        
    ## Relative vorticity
    F = RegularGridInterpolator((model['x'],model['y'],model['z']),model['vor'])
    parti['vor'] = F((parti.x,parti.y,parti.z)) 
    parti['wtotal'] = parti.wsink+parti.w
    
    ##
    print('DONE')
    print('#################')
    print(' ')

    return parti
