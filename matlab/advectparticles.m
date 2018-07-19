% This routine projects the particle forward

disp('PROJECT PARTICLE DATA... ')

%% STEPS
%       - Get the velocity fluxes at the faces surrounding the particle's
%           location
%       - Compute the variables necessary for the differential equation
%           (Beta and delta)
%       - Compute the shortest time it would take for the particle to reach
%           one of the cell face (in x, y, or z)
%           if the time is shorter than the particle tracking timestep,
%           then divide the timestep into subtimesteps
%       - Get new particle position

%h = waitbar(0,'Advection of particles, please wait...');

counter_remover = 0;
end_loop = length(parti.x);
for xx = 1:end_loop
    ii = xx-counter_remover;
    
    %waitbar(ii / length(parti.x),h,['Advection of particles, please wait... [',num2str(ii),'/',num2str(length(parti.x)),']'])
    
    % ii can be larger than length(parti.x) only because some values are
    % removed when particles are outside of the domain. so parti.x possibly
    % changes length within the loop
    %if ii > length(parti.x)+counter_remover
    %    break
    %end
    
    timestep_leftover = particle.timestep*86400;
    while timestep_leftover~=0
        % Find indices of cell faces surrounding the particle
        % position in x
        face_im1 = find(model.xf<=parti.x(ii),1,'last');
        face_i = face_im1+1;
        face_jm1 = find(model.yf<=parti.y(ii),1,'last');
        face_j = face_jm1+1;
        face_km1 = find(model.zf<=parti.z(ii),1,'last');
        face_k = face_km1+1;
        
        % Correct indices if particle on a cell face. For instance, if a
        % particle is on a cell face in x, the cell to be considered for
        % advection is the one to the right of the particle if the face
        % velocity is +ve, and to the left of the particle is face velocity
        % is -ve. The cell to the right is considered by default (see
        % above), and this piece of code modifies this if the face velocity
        % is negative. This reasoning is reversed if simulation is ran
        % backwards
        
        % define direction of the tracking
        if strcmp(particle.direction,'forward')
            direction = 1;
        else
            direction = -1;
        end
        
        % These are the conditional indexing used for every face velocity.
        % conditional indexing is time consuming so indices are identified
        % separately and once.
        indxf = find(model.x(2:end-1)>=model.xf(face_im1) & model.x(2:end-1)<=model.xf(face_i));
        indyf = find(model.y(2:end-1)>=model.yf(face_jm1) & model.y(2:end-1)<=model.yf(face_j));
        indzf = find(model.z(2:end-1)>=model.zf(face_km1) & model.z(2:end-1)<=model.zf(face_k));
        
        if parti.x(ii)-model.xf(face_im1) == 0 &&...
                direction*model.uf(face_im1,indyf,indzf)<0
            
            % Introduce periodicity at western boundary
            if face_im1 == 1 && model.periodic_ew ==1
                % Select the last cell instead of the first one
                face_im1 = length(model.xf)-1;
                face_i = length(model.xf);
                parti.x(ii) = model.xf(face_i);
                
                % If not at model boundary
            else
                face_im1 = face_im1-1;
                face_i = face_i-1;
            end
            
            % NOTE: Periodicity at eastern boundary is not necessary, as
            % the particle location in x is modular of the model domain
            % (see mod function below). so when eastern most cell face is
            % reached, particle location is changed western most cell
            
        end
        
        if parti.y(ii)-model.yf(face_jm1) == 0 &&...
                direction*model.vf(indxf,face_jm1,indzf)<0
            
            % If the particle is located at the very 1st cell, indices
            % cannot be switched to the previous cell (southward).
            if face_jm1 == 1
                
                % If model is periodic in y, then periodicity is introduced
                if model.periodic_ns ==1
                    % Select the last cell instead of the first one
                    face_jm1 = length(model.yf)-1;
                    face_j = length(model.yf);
                    parti.y(ii) = model.yf(face_j);
                else
                end
                
                
                % If not at model boundary
            else
                face_jm1 = face_jm1-1;
                face_j = face_j-1;
            end
            
            % NOTE: Periodicity at northern boundary is not necessary, as the
            % particle location in y is modular of the model domain (see mod
            % function below). so when northern most cell face is reached, particle
            % location is changed to southern most cell
            
        end
        
        if parti.z(ii)-model.zf(face_km1) == 0 &&...
                direction*model.wf(indxf,indyf,face_km1)<0
            
            face_km1 = face_km1-1;
            face_k = face_k-1;
        end
        
        % Update useful variables
        indxf = find(model.x(2:end-1)>=model.xf(face_im1) & model.x(2:end-1)<=model.xf(face_i));
        indyf = find(model.y(2:end-1)>=model.yf(face_jm1) & model.y(2:end-1)<=model.yf(face_j));
        indzf = find(model.z(2:end-1)>=model.zf(face_km1) & model.z(2:end-1)<=model.zf(face_k));
        
        % Compute Cell spacings in x, y, and z
        % Compute Dx
        Dx = model.xf(face_i)-model.xf(face_im1);
        % Compute Dy
        Dy = model.yf(face_j)-model.yf(face_jm1);
        % Compute Dz
        Dz = model.zf(face_k)-model.zf(face_km1);
        
        %%
        % Compute some variables in x
        rx0 = parti.x(ii)/Dx;
        rxim1 = model.xf(face_im1)/Dx;
        rxi = model.xf(face_i)/Dx;
        betax = model.uf(face_im1,indyf,indzf) - ...
            model.uf(face_i,indyf,indzf);
        deltax = -model.uf(face_im1,indyf,indzf) -...
            betax*rxim1;
        
        % Compute some variables in y
        ry0 = parti.y(ii)/Dy;
        ryjm1 = model.yf(face_jm1)/Dy;
        ryj = model.yf(face_j)/Dy;
        betay = model.vf(indxf,face_jm1,indzf) - ...
            model.vf(indxf,face_j,indzf);
        deltay = -model.vf(indxf,face_jm1,indzf) -...
            betay*ryjm1;
        
        % Compute some variables in z
        rz0 = parti.z(ii)/Dz;
        rzkm1 = model.zf(face_km1)/Dz;
        rzk = model.zf(face_k)/Dz;
        betaz = model.wf(indxf,indyf,face_km1) - ...
            model.wf(indxf,indyf,face_k);
        deltaz = -model.wf(indxf,indyf,face_km1) -...
            betaz*rzkm1;
        
        %% Compute the shortest time it would take the particle to reach one of the
        % cell faces
        
        % time to reach i-1 th face
        Dtmaxtemp(1) = -1/betax*log((rxim1 + deltax/betax)/(rx0 + deltax/betax))*Dx*Dy*Dz;
        % time to reach i th face
        Dtmaxtemp(2) = -1/betax*log((rxi + deltax/betax)/(rx0 + deltax/betax))*Dx*Dy*Dz;
        % time to reach j-1 th face
        Dtmaxtemp(3) = -1/betay*log((ryjm1 + deltay/betay)/(ry0 + deltay/betay))*Dx*Dy*Dz;
        % time to reach j th face
        Dtmaxtemp(4) = -1/betay*log((ryj + deltay/betay)/(ry0 + deltay/betay))*Dx*Dy*Dz;
        % time to reach k-1 th face
        Dtmaxtemp(5) = -1/betaz*log((rzkm1 + deltaz/betaz)/(rz0 + deltaz/betaz))*Dx*Dy*Dz;
        % time to reach j th face
        Dtmaxtemp(6) = -1/betaz*log((rzk + deltaz/betaz)/(rz0 + deltaz/betaz))*Dx*Dy*Dz;
        
        
        % Imaginary times occur when the face considered has a velocity in
        % the opposite direction of the velocity at the particle's
        % location. Hence the imaginary number: does not matter how long,
        % the particle will never reach that face.
        
        if strcmp(particle.direction,'forward')
            [Dtmax,~] = min(Dtmaxtemp(Dtmaxtemp>0 & imag(Dtmaxtemp)==0));
            % If the shortest time it would take the particle to reach a cell face
            % is shorter than the particle tracking timestep, then the integration
            % timestep needs to be shortened
            if isempty(Dtmax)
                intermediate_timestep = timestep_leftover;
            else
                intermediate_timestep = min(Dtmax,timestep_leftover);
            end
            timestep_leftover = timestep_leftover-intermediate_timestep;
            ds = intermediate_timestep/(Dx*Dy*Dz);
            
            
        elseif strcmp(particle.direction,'backward')
            Dtmax = max(Dtmaxtemp(Dtmaxtemp<0 & imag(Dtmaxtemp)==0));
            
            % If the shortest time it would take the particle to reach a cell face
            % is shorter than the particle tracking timestep, then the integration
            % timestep needs to be shortened
            if isempty(Dtmax)
                intermediate_timestep = timestep_leftover;
            else
                intermediate_timestep = max(Dtmax,timestep_leftover);
            end
            timestep_leftover = timestep_leftover-intermediate_timestep;
            ds = intermediate_timestep/(Dx*Dy*Dz);
            
        end
        
        %=================
        %=== Assign x-position to particle.
        %=================
        rx1 = (rx0 + deltax/betax)*exp(-betax*ds)-(deltax/betax);
        % If particle is really close to boundary, then particle is
        % assigned the boundary value.
        if abs(rx1-rxim1)<1e-11
            %rx1 = rxim1;
            parti.x(ii) = model.xf(face_im1);
            
        elseif abs(rx1-rxi)<1e-11
            %rx1 = rxi;
            parti.x(ii) = model.xf(face_i);
        else
            parti.x(ii) = rx1*Dx;
        end
        
        %=================
        %=== Assign y-position to particle.
        %=================
        ry1 = (ry0 + deltay/betay)*exp(-betay*ds)-(deltay/betay);
        % If particle is really close to boundary, then particle is
        % assigned the boundary value.
        if abs(ry1-ryjm1)<1e-11
            %ry1 = ryjm1;
            parti.y(ii) = model.yf(face_jm1);
            
        elseif abs(ry1-ryj)<1e-11
            %ry1 = ryj;
            parti.y(ii) = model.yf(face_j);
        else
            parti.y(ii) = ry1*Dy;
        end
        
        %=================
        %=== Assign z-position to particle.
        %=================
        rz1 = (rz0 + deltaz/betaz)*exp(-betaz*ds)-(deltaz/betaz);
        
        % If particle is really close to boundary, then particle is
        % assigned the boundary value.
        if abs(rz1-rzkm1)<1e-11
            %rz1 = rzkm1;
            parti.z(ii) = model.zf(face_km1);
        elseif abs(rz1-rzk)<1e-11
            %rz1 = rzk;
            parti.z(ii) = model.zf(face_k);
        else
            parti.z(ii) = rz1*Dz + parti.wsink(ii)*intermediate_timestep;
        end
        
        if model.periodic_ew == 1
            % Wrap the domain around if periodicity
            parti.x(ii) = mod(parti.x(ii),max(model.xf));
        end
        if model.periodic_ns == 1
            % Wrap the domain around if periodicity
            parti.y(ii) = mod(parti.y(ii),max(model.yf));
        end
        
        % Particles can't go airborn (z>0), and stop tracking if too deep
        if parti.z(ii)>0
            parti.z(ii) = 0;
        elseif parti.z(ii)<min(model.zf)
            thefields = fieldnames(parti);
            for jj = 1:length(thefields)
                eval(['parti.',thefields{jj},'(ii) = [];'])
            end; clear jj
            counter_remover = counter_remover + 1;
            break
        end
        
        % If hits a solid boundary, kill the particle
        if parti.y(ii)>model.yf(end) || parti.y(ii)<model.yf(1)
            warning('Particle removed at boundary!')
            thefields = fieldnames(parti);
            for jj = 1:length(thefields)
                eval(['parti.',thefields{jj},'(ii) = [];'])
            end; clear jj
            counter_remover = counter_remover + 1;
            break
        end
        
    end
end; clear ii

% remove useless fields
clear thefields beta* delta* center_* face_* direction ds Dtmax Dx Dy Dz intermediate_timestep rx0 rx1 rxi rxim1 ry0 ry1 ryj ryjm1 rz0 rz1 rzk rzkm1 Dtmaxtemp indxf indyf indzf

disp(['Removed ',num2str(counter_remover),' particles'])
%close(h)
clear timestep_leftover counter_remover h