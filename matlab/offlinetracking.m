clear

%% This routine conducts offline particle tracking based on imported 3D
% velocity field

%% PARAMETERS

%%%%%%%%%%%%%%%%%%%%%
%% PSOM parameters
%%%%%%%%%%%%%%%%%%%%%
% Day of year where simulation starts
model.start_day = 0;
% Model timestep in second (needed to convert from timestamp to DOY
model.timestep = 216;
% Path of folder containing model outputs
model.path = '/Volumes/garuda-1/Mathieu/output_Papa_Dx=500_from_1km/K1/highrez/';
% constant to play with model output temporal resolution
% "2" means every other model output file
% "4" means every 4 model output file
% "1" means every file
model.outputskip = 1;

% Define periodicities in model domain
model.periodic_ew = 1;
model.periodic_ns = 0;

%%%%%%%%%%%%%%%%%%%%%
%% Particle parameters
%%%%%%%%%%%%%%%%%%%%%

%----------------
% SEEDING
%----------------
% day of year for particle initialization (i.e., first seeding)
particle.initime = 60.5;
% frequency of particles seeding (in days)
particle.inifreq = 2^-2;
% Number of seeding events
particle.ininumber = 1;
% Type of particle seeding ('dynamic' or 'static')
particle.initype = 'static';

% Number of particle classes to be released. Please refer to iniparticles
% to customize the sinking velocity of each particle class
particle.numofclasses = 1;

% Parameters to customize the initialization of the particles (in meters)
% Only applied if seeding type is 'static'
particle.istart = 50000;
particle.irange = 1000;
particle.irez = 2000;

particle.jstart = 150000;
particle.jrange = 1000;
particle.jrez = 2000;

particle.kstart = -10;
particle.krange = 1;
particle.krez = 3;

%----------------
% ADVECTION
%----------------
% timestep of particle tracking in days
particle.timestep = 2^-4;
% length of particle tracking experiment in days
particle.length = 5;
% Define particle-tracking direction ('forward' or 'backward')
particle.direction = 'forward';

%----------------
% OUTPUT
%----------------
% frequency of particles output (in days)
particle.outfreq = 2^-3;
% Format of output ('csv' or 'sqlite')
particle.outputformat = 'sqlite';
particle.outputdir = '/Users/mathieudever/Desktop/particle_test/MATLAB';
particle.outputfilename = 'matlab_test_1parti_nopara';

%%%%%%%%%%%%%%%%%%%%%
%% CORE CODE
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
    
    % If day of year matches the seeding frequency (threshold is to account
    % for rouding error is inifreq is not a power of 2).
    if  mod(tt-particle.initime,particle.inifreq)<=1e-10 && particle.ininumber ~= 0

        % If seeding is set to 'dynamic', the iniparticle routine uses the
        % current particle position of the 1st particle to seed (at a 
        % constant depth)
        if strcmp(particle.initype,'dynamic') && tt ~= particle.initime
            particle.istart = parti.x(parti.id==1);
            particle.irange = 1;
            particle.irez = 1;
            
            particle.jstart = parti.y(parti.id==1);
            particle.jrange = 1;
            particle.jrez = 1;
            
            particle.kstart = -50;
            particle.krange = 1;
            particle.krez = 1;
        end
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
        %run advectparticles_parallel
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
