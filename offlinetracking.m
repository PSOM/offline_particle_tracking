clear

%% This routine conducts offline particle tracking based on imported 3D
% velocity field

%% PARAMETERS

%%%%%%%%%%%%%%%%%%%%%
% PSOM parameters
%%%%%%%%%%%%%%%%%%%%%
% Day of year where simulation starts
model.start_day = 90;
% Model timestep in second (needed to convert from timestamp to DOY
model.timestep = 432;
% Path of folder containing model outputs
model.path = '/Volumes/garuda/Mathieu/output_Papa1km/particle_tracking';
% constant to play with model output temporal resolution
% "2" means every other model output file
% "4" means every 4 model output file
% "1" means every file
model.outputskip = 1;

% Define periodicities in model domain
model.periodic_ew = 1;
model.periodic_ns = 0;

%%%%%%%%%%%%%%%%%%%%%
% Particle parameters
%%%%%%%%%%%%%%%%%%%%%
% day of year for particle initialization (i.e., first seeding)
particle.initime = 183;
% frequency of particles seeding (in days)
particle.inifreq = 100;
% timestep of particle tracking in days
particle.timestep = 2^-5;
% Number of seeding events
particle.ininumber = 12;
% length of particle tracking experiment in days
particle.length = 60;
% frequency of particles output (in days)
particle.outfreq = 2^-2;
% Format of output ('csv' or 'sqlite')
particle.outputformat = 'sqlite';
particle.outputdir = '/Users/mathieudever/Documents/EXPORTS/Particles/offline_particle_tracking/Papa_1km';
particle.outputfilename = 'offlineparticles_6hrsPSOM_full_relax3days';

% define particle-tracking direction ('forward' or 'backward')
particle.direction = 'forward';

% parameters to customize the initialization of the particles (in meters)
particle.istart = 10000;
particle.irange = 1000;
particle.irez = 1000;
particle.jstart = 100000;
particle.jrange = 100000;
particle.jrez = 1000;
particle.kstart = -25;
particle.krange = 1;
particle.krez = 1;
particle.numofclasses = 1;

%%%%%%%%%%%%%%%%%%%%%
% CORE CODE
%%%%%%%%%%%%%%%%%%%%%
mytimer = tic;

% Change a few things if backward-tracking
if strcmp(particle.direction,'backward')
    particle.timestep = -particle.timestep;
    particle.length = -particle.length;
end

% Print the model configuration to a text-file
run saveconfiguration

for tt = particle.initime:particle.timestep:particle.initime+particle.length
    
    disp('')
    disp([num2str((tt-particle.initime)/(particle.length)*100),'% DONE (doy ',num2str(tt),')'])
    
    % Start timer
    tstart = tic;
    
    %% Get the velocity field relevant to the time step considered, and 
    % already interpolated onto the time step

    run getvelocityfield
    
    %% Advect particles backward in time if backward particle tracking
    
    % Skip a time step in the case of backward particle tracking. this
    % is because the forward case uses the velocity field at time t to
    % compute trajectories from t to t+Dt, so backward tracking should
    % use the same velocity field to compute the trajectory from t+Dt
    % to t.
    if strcmp(particle.direction,'backward') && tt == particle.initime
    % Skip a time step
    elseif strcmp(particle.direction,'backward')
        run advectparticles
    end
    
    %% Seed particle when required
    
    % If day of year matches the seeding frequency
    if  mod(tt-particle.initime,particle.inifreq)<=1e-10
    % If day of year matches the seeding frequency (threshold is to account
    % for rouding error is inifreq is not a power of 2).
    if  mod(tt-particle.initime,particle.inifreq)<=1e-10 && particle.ininumber ~= 0
        % Seed particles
        run iniparticles
        particle.ininumber = particle.ininumber -1;
    end
    
    % Display total number of particle in the experiment
    disp(['Total number of particles is ',num2str(length(parti.x))])
    
    %% Get the model variables interpolated onto the particle's position
    % Interpolate model variables onto particle's position in 4D
    run getpartivar
    
    %% Save particle position when required

    if mod(tt-particle.initime,particle.outfreq)<=1e-10
        run saveparti
    end
    
    %% Advect particles forward in time if forward particle tracking
    
    if strcmp(particle.direction,'forward')
        run advectparticles
    end
    
    %% Displays some info to user
    
    % End timer
    telapsed = toc(tstart);
    disp(['This step took ',num2str(telapsed/60),' minutes'])
    disp(['At this rate, you have ', num2str(length(tt:particle.timestep:particle.initime+particle.length)*telapsed/60),' minutes left']);
    
end; clear tt
totaltime = toc(mytimer);
disp(['Entire run took ', num2str(totaltime/60),' minutes']);
