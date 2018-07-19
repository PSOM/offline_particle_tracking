#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#import line_profiler

#@profile
def advectparticles(tt,particle,model,parti):
 #%%
    """
    This routine projects the particle forward.
    
    
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
    
    print('PROJECT PARTICLE DATA... ')    
  
    for ii in range(0,len(parti.x)):        
        timestep_leftover = particle.timestep*86400
        while timestep_leftover!=0:
            #print(timestep_leftover)
            #%% THIS SECTION IDENTIFIES THE PROPER MODEL CELL TO BE USED FOR ADVECTION
            
            # Find indices of cell faces surrounding the particle
            # position in x
            face_im1 = np.where(model['xf']<=parti.at[ii,'x'])[0][-1]
            face_i = face_im1+1
            face_jm1 = np.where(model['yf']<=parti.at[ii,'y'])[0][-1]
            face_j = face_jm1+1
            face_km1 = np.where(model['zf']<=parti.at[ii,'z'])[0][-1]
            face_k = face_km1+1
        
            # Compute some useful variables
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
            
            # Correct indices if particle on a cell face. For instance, if a
            # particle is on a cell face in x, the cell to be considered for
            # advection is the one to the right of the particle if the face
            # velocity is +ve, and to the left of the particle is face velocity
            # is -ve. The cell to the right is considered by default (see
            # above), and this piece of code modifies this if the face velocity
            # is negative.
            
            # define direction of the tracking
            if particle.direction=='forward':
                direction = 1
            elif particle.direction=='backward':
                direction = -1
            else:
                print('Direction must be forward or backward.')
                exit
            
            #A
            if direction*(model['uf'][face_im1,i2,i3])<0:
                if parti.at[ii,'x']-model_im1 == 0:
                
                    print('block A') 
                # Introduce periodicity at western boundary
                    if face_im1 == 0 and model.periodic_ew ==1:
                    # Select the last cell instead of the first one
                        face_im1 = len(model['xf'])-2
                        face_i = len(model['xf'])-1
                        parti.at[ii,'x'] = model_i
                    
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
                if parti.at[ii,'y']-model_jm1 == 0:
                
                    print('block B') 
                # If the particle is located at the very 1st cell, indices
                # cannot be switched to the previous cell (southward).
                    if face_jm1 == 0:
                    # If model is periodic in y, then periodicity is introduced
                        if model.periodic_ns ==1:
                        # Select the last cell instead of the first one
                            face_jm1 = len(model['yf'])-2
                            face_j = len(model['yf'])-1
                            parti.at[ii,'y'] = model_j
    
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
                if parti.at[ii,'z']-model_km1 == 0:
                
                    print('block C') 
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
            
            #%%
            # Compute Cell spacings in x, y, and z
            # Compute Dx,Dy,Dz
            Dx = model_i-model_im1
            Dy = model_j-model_jm1
            Dz = model_k-model_km1
            
            # Compute some variables in x
            rx0 = parti.at[ii,'x']/Dx
            rxim1 = model_im1/Dx
            rxi = model_i/Dx
            betax = model['uf'][face_im1,i2,i3] - model['uf'][face_i,i2,i3]
            deltax = -model['uf'][face_im1,i2,i3]- betax*rxim1
            
            # Compute some variables in y
            ry0 = parti.at[ii,'y']/Dy
            ryjm1 = model_jm1/Dy
            ryj = model_j/Dy
            betay = model['vf'][i1,face_jm1,i3] - model['vf'][i1,face_j,i3]
            deltay = -model['vf'][i1,face_jm1,i3]-betay*ryjm1
            
            # Compute some variables in z
            rz0 = parti.at[ii,'z']/Dz
            rzkm1 = model_km1/Dz
            rzk = model_k/Dz
            betaz = model['wf'][i1,i2,face_km1] - model['wf'][i1,i2,face_k]
            deltaz = -model['wf'][i1,i2,face_km1] - betaz*rzkm1
            
            ## Compute the shortest time it would take the particle to reach one of the
            # cell faces
            
            # time to reach i-1 th face
            Dtmaxtemp = np.zeros(6,dtype='complex64')
            invbetax = 1/betax
            invbetay = 1/betay
            invbetaz = 1/betaz
            
            Dtmaxtemp[0] = -invbetax*np.log((rxim1 + deltax*invbetax)/(rx0 + deltax*invbetax))*Dx*Dy*Dz
            # time to reach i th face
            Dtmaxtemp[1] = -invbetax*np.log((rxi + deltax*invbetax)/(rx0 + deltax*invbetax))*Dx*Dy*Dz
            # time to reach j-1 th face
            Dtmaxtemp[2] = -invbetay*np.log((ryjm1 + deltay*invbetay)/(ry0 + deltay*invbetay))*Dx*Dy*Dz
            # time to reach j th face
            Dtmaxtemp[3] = -invbetay*np.log((ryj + deltay*invbetay)/(ry0 + deltay*invbetay))*Dx*Dy*Dz
            # time to reach k-1 th face
            Dtmaxtemp[4] = -invbetaz*np.log((rzkm1 + deltaz*invbetaz)/(rz0 + deltaz*invbetaz))*Dx*Dy*Dz
            # time to reach j th face
            Dtmaxtemp[5] = -invbetaz*np.log((rzk + deltaz*invbetaz)/(rz0 + deltaz*invbetaz))*Dx*Dy*Dz
            
            # Imaginary times occur when the face considered has a velocity in
            # the opposite direction of the velocity at the particle's
            # location. Hence the imaginary number: does not matter how long,
            # the particle will never reach that face.
    #%           
            if particle.direction=='forward':
                Dtmax = Dtmaxtemp[ (Dtmaxtemp>0) & (Dtmaxtemp.imag==0)].min().real
                # If the shortest time it would take the particle to reach a cell face
                # is shorter than the particle tracking timestep, then the integration
                # timestep needs to be shortened
                if not Dtmax:
                    intermediate_timestep = timestep_leftover
                else:
                    intermediate_timestep = np.min([Dtmax,timestep_leftover])
                    
                timestep_leftover = timestep_leftover-intermediate_timestep
                ds = intermediate_timestep/(Dx*Dy*Dz)
                
                
            elif particle.direction=='backward':
                Dtmax = Dtmaxtemp[ (Dtmaxtemp>0) & (Dtmaxtemp.imag==0)].min().real
                Dtmax = Dtmaxtemp[ (Dtmaxtemp<0) & (Dtmaxtemp.imag==0)].min().real
                
                # If the shortest time it would take the particle to reach a cell face
                # is shorter than the particle tracking timestep, then the integration
                # timestep needs to be shortened
                if not Dtmax:
                    intermediate_timestep = timestep_leftover
                else:
                    intermediate_timestep = np.max([Dtmax,timestep_leftover])
                
                timestep_leftover = timestep_leftover-intermediate_timestep
                ds = intermediate_timestep/(Dx*Dy*Dz)
                
            #=================
            #=== Assign x-position to particle.
            #=================
            rx1 = (rx0 + deltax*invbetax)*np.exp(-betax*ds)-(deltax*invbetax)
            if abs(rx1-rxim1)<1e-11:
                parti.at[ii,'x'] = model_im1
            elif abs(rx1-rxi)<1e-11:
                parti.at[ii,'x'] = model_i
            else:
                parti.at[ii,'x'] = rx1*Dx
            
            #=================
            #=== Assign y-position to particle.
            #=================
            ry1 = (ry0 + deltay*invbetay)*np.exp(-betay*ds)-(deltay*invbetay)
            if abs(ry1-ryjm1)<1e-11:
                parti.at[ii,'y'] = model_jm1
            elif abs(ry1-ryj)<1e-11:
                parti.at[ii,'y'] = model_j
            else:
                #parti.loc[ii,'y'] = ry1*Dy
                parti.at[ii,'y'] = ry1*Dy
            #=================
            #=== Assign z-position to particle.
            #=================
            rz1 = (rz0 + deltaz*invbetaz)*np.exp(-betaz*ds)-(deltaz*invbetaz)
            if abs(rz1-rzkm1)<1e-11:
                parti.at[ii,'z'] = model_km1
            elif abs(rz1-rzk)<1e-11:
                parti.at[ii,'z'] = model_k
            else:
                parti.at[ii,'z'] = rz1*Dz + parti.at[ii,'wsink']*intermediate_timestep
            
            if model['periodic_ew'] == 1:
                # Wrap the domain around if periodicity
                parti.at[ii,'x'] = parti.at[ii,'x'] % model['xf'].max()
            
            if model['periodic_ns'] == 1:
                # Wrap the domain around if periodicity
                parti.at[ii,'y'] = parti.at[ii,'y'] % model['yf'].max()
            
            # Particles can't go airborn (z>0), and stop tracking if too deep
            #if parti.at[ii,'z']>0:
            #    parti.at[ii,'z'] = 0
            #elif parti.at[ii,'z']<model['zf'].min():
            #    parti.drop(parti.index[ii],inplace=True)
            #    counter_remover = counter_remover + 1
            #    break
            
            # If hits a solid boundary, kill the particle
            #if parti.at[ii,'y']>model['yf'][-1] or parti.at[ii,'y']<model['yf'][0]:
            #    parti.drop(parti.index[ii],inplace=True)
            #    counter_remover = counter_remover + 1
            #    break

    return parti