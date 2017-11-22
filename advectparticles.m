% This routine projects the particle forward using a 2nd order
% Adams?Bashforth method
% https://en.wikipedia.org/wiki/Linear_multistep_method

disp('PROJECT PARTICLE DATA... ')

%% STEP 1
%       - Get the velocity fluxes at the faces surrounding the particle's
%           location
%       - Compute the variables necessary for the differential equation
%           (Beta and delta)
%       - Compute the shortest time it would take for the particle to reach
%           one of the cell face (in x, y, or z)
%           if the time is shorter than the particle tracking timestep,
%           then divide the timestep into subtimesteps

counter_remover = 0;
for ii = 1:length(parti.x)
    
    % ii can be larger than length(parti.x) only because some values are
    % removed when particles are outside of the domain. so parti.x possibly
    % changes length within the loop
    if ii > length(parti.x)
        break
    end
    
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
        
        % Compute Dx
        Dx = model.xf(face_i)-model.xf(face_im1);
        
        % Compute Dy
        Dy = model.yf(face_j)-model.yf(face_jm1);
        
        % Compute Dz
        Dz = model.zf(face_k)-model.zf(face_km1);
        
        % Correct indices if particle on a cell face. For instance, if a particle
        % is on a cell face in x, the cell to be considered for advection is the
        % one to the right of th particle if the face velocity is +ve, and to the
        % left of the particle is face velocity is -ve.
        % The cell to the right is considered by default, and this piece of code
        % modifies this if the face velocity is negative.
        
        % define direction of the tracking
        if strcmp(particle.direction,'forward')
            direction = 1;
        else
            direction = -1;
        end
        
        if parti.x(ii)-model.xf(face_im1) == 0 &&...
                direction*model.uf(face_im1,...
                model.y(2:end-1)>=model.yf(face_jm1) & model.y(2:end-1)<=model.yf(face_j),...
                model.z(2:end-1)>=model.zf(face_km1) & model.z(2:end-1)<=model.zf(face_k))<0
            
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
                direction*model.vf(model.x(2:end-1)>=model.xf(face_im1) & model.x(2:end-1)<=model.xf(face_i),...
                face_jm1,...
                model.z(2:end-1)>=model.zf(face_km1) & model.z(2:end-1)<=model.zf(face_k))<0
            
            % Introduce periodicity at southern boundary
            if face_jm1 == 1 && model.periodic_ns ==1
                % Select the last cell instead of the first one
                face_jm1 = length(model.yf)-1;
                face_j = length(model.yf);
                parti.y(ii) = model.yf(face_j);
                
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
                direction*model.wf(model.x(2:end-1)>=model.xf(face_im1) & model.x(2:end-1)<=model.xf(face_i),...
                model.y(2:end-1)>=model.yf(face_jm1) & model.y(2:end-1)<=model.yf(face_j),...
                face_km1)<0
            
            face_km1 = face_km1-1;
            face_k = face_k-1;
        end
        
        
        %%
        % Compute some variables in x
        rx0 = parti.x(ii)/Dx;
        rxim1 = model.xf(face_im1)/Dx;
        rxi = model.xf(face_i)/Dx;
        betax = model.uf(face_im1,...
            model.y(2:end-1)>=model.yf(face_jm1) & model.y(2:end-1)<=model.yf(face_j),...
            model.z(2:end-1)>=model.zf(face_km1) & model.z(2:end-1)<=model.zf(face_k)) - ...
            model.uf(face_i,...
            model.y(2:end-1)>=model.yf(face_jm1) & model.y(2:end-1)<=model.yf(face_j),...
            model.z(2:end-1)>=model.zf(face_km1) & model.z(2:end-1)<=model.zf(face_k));
        deltax = -model.uf(face_im1,...
            model.y(2:end-1)>=model.yf(face_jm1) & model.y(2:end-1)<=model.yf(face_j),...
            model.z(2:end-1)>=model.zf(face_km1) & model.z(2:end-1)<=model.zf(face_k)) -...
            betax*rxim1;
        
        % Compute some variables in y
        ry0 = parti.y(ii)/Dy;
        ryjm1 = model.yf(face_jm1)/Dy;
        ryj = model.yf(face_j)/Dy;
        betay = model.vf(model.x(2:end-1)>=model.xf(face_im1) & model.x(2:end-1)<=model.xf(face_i),...
            face_jm1,...
            model.z(2:end-1)>=model.zf(face_km1) & model.z(2:end-1)<=model.zf(face_k)) - ...
            model.vf(model.x(2:end-1)>=model.xf(face_im1) & model.x(2:end-1)<=model.xf(face_i),...
            face_j,...
            model.z(2:end-1)>=model.zf(face_km1) & model.z(2:end-1)<=model.zf(face_k));
        deltay = -model.vf(model.x(2:end-1)>=model.xf(face_im1) & model.x(2:end-1)<=model.xf(face_i),...
            face_jm1,...
            model.z(2:end-1)>=model.zf(face_km1) & model.z(2:end-1)<=model.zf(face_k)) -...
            betay*ryjm1;
        
        % Compute some variables in z
        rz0 = parti.z(ii)/Dz;
        rzkm1 = model.zf(face_km1)/Dz;
        rzk = model.zf(face_k)/Dz;
        betaz = model.wf(model.x(2:end-1)>=model.xf(face_im1) & model.x(2:end-1)<=model.xf(face_i),...
            model.y(2:end-1)>=model.yf(face_jm1) & model.y(2:end-1)<=model.yf(face_j),...
            face_km1) - ...
            model.wf(model.x(2:end-1)>=model.xf(face_im1) & model.x(2:end-1)<=model.xf(face_i),...
            model.y(2:end-1)>=model.yf(face_jm1) & model.y(2:end-1)<=model.yf(face_j),...
            face_k);
        deltaz = -model.wf(model.x(2:end-1)>=model.xf(face_im1) & model.x(2:end-1)<=model.xf(face_i),...
            model.y(2:end-1)>=model.yf(face_jm1) & model.y(2:end-1)<=model.yf(face_j),...
            face_km1) -...
            betaz*rzkm1;
%         
%         rz0 = -parti.z(ii)/Dz;
%         rzkm1 = -model.zf(face_km1)/Dz;
%         rzk = -model.zf(face_k)/Dz;
%         betaz = model.wf(model.x(2:end-1)>=model.xf(face_im1) & model.x(2:end-1)<=model.xf(face_i),...
%             model.y(2:end-1)>=model.yf(face_jm1) & model.y(2:end-1)<=model.yf(face_j),...
%             face_km1) - ...
%             model.wf(model.x(2:end-1)>=model.xf(face_im1) & model.x(2:end-1)<=model.xf(face_i),...
%             model.y(2:end-1)>=model.yf(face_jm1) & model.y(2:end-1)<=model.yf(face_j),...
%             face_k);
%         betaz=-betaz;
%         deltaz = model.wf(model.x(2:end-1)>=model.xf(face_im1) & model.x(2:end-1)<=model.xf(face_i),...
%             model.y(2:end-1)>=model.yf(face_jm1) & model.y(2:end-1)<=model.yf(face_j),...
%             face_km1) -...
%             betaz*rzkm1;
        
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
        
        
        % This bit gets rid of imaginary numbers. Imaginary times occur
        % when the the face considered has a velocity in the opposite
        % direction of the velocity at the particle's location. Hence the
        % imaginary number: does not matter how long, the particle will
        % never reach that face.
%         if isreal(Dtmaxtemp)==0
%             warning('imaginary time...')
%             error('blop')
%         end
        
        if strcmp(particle.direction,'forward')
            % = NaN;
            [Dtmax,~] = min(Dtmaxtemp(Dtmaxtemp>0));
            
            % If the shortest time it would take the particle to reach a cell face
            % is shorter than the particle tracking timestep, then the integration
            % timestep needs to be shortened
            intermediate_timestep = min(Dtmax,particle.timestep*86400);
            timestep_leftover = particle.timestep*86400-intermediate_timestep;
            ds = intermediate_timestep/(Dx*Dy*Dz);
            
            
        elseif strcmp(particle.direction,'backward')
            Dtmax = max(Dtmaxtemp(Dtmaxtemp<0));
            
            % If the shortest time it would take the particle to reach a cell face
            % is shorter than the particle tracking timestep, then the integration
            % timestep needs to be shortened
            intermediate_timestep = max(Dtmax,particle.timestep*86400);
            timestep_leftover = particle.timestep*86400-intermediate_timestep;
            ds = intermediate_timestep/(Dx*Dy*Dz);
            
        end
        
        %=================
        %=== Assign x-position to particle.
        %=================
        rx1 = (rx0 + deltax/betax)*exp(-betax*ds)-(deltax/betax);
         if abs(rx1-rxim1)<1e-11
             rx1 = rxim1;
         end
         if abs(rx1-rxi)<1e-11
             rx1 = rxi;
         end
        parti.x(ii) = rx1*Dx;
        
        %=================
        %=== Assign y-position to particle.
        %=================
        ry1 = (ry0 + deltay/betay)*exp(-betay*ds)-(deltay/betay);
        if abs(ry1-ryjm1)<1e-11
            ry1 = ryjm1;
        elseif abs(ry1-ryj)<1e-11
            ry1 = ryj;
        end
        parti.y(ii) = ry1*Dy;
        
        %=================
        %=== Assign z-position to particle.
        %=================
        rz1 = (rz0 + deltaz/betaz)*exp(-betaz*ds)-(deltaz/betaz);
         if abs(rz1-rzkm1)<1e-11
             rz1 = rzkm1;
         elseif abs(rz1-rzk)<1e-11
             rz1 = rzk;
         end
        parti.z(ii) = rz1*Dz;
        
        if model.periodic_ew == 1
            % Wrap the domain around if periodicity
            parti.x(ii) = mod(parti.x(ii),max(model.xf));
        end
        if model.periodic_ns == 1
            % Wrap the domain around if periodicity
            parti.x(ii) = mod(parti.x(ii),max(model.xf));
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
        
        clear thefields beta* delta* center_* face_* direction ds Dtmax Dx Dy Dz intermediate_timestep rx0 rx1 rxi rxim1 ry0 ry1 ryj ryjm1 rz0 rz1 rzk rzkm1 Dtmaxtemp
    end
end; clear ii

disp(['Removed ',num2str(counter_remover),' particles'])

clear timestep_leftover counter_remover