%% This routine prints the particle-tracking experiment configuration to a
% text file. The filename includes the date and time of the simulation.

disp('SAVE EXPERIMENT CONFIGURATION... ')

% Get current time
nowtime = clock;
nowtime = datestr(nowtime,'ddmmmyy-HHMM');

% Open the file
fileID = fopen([particle.outputdir,'/',nowtime,'particle_tracking_configuration','.out'],'w');

% Write the data
fprintf(fileID,'%s\n','******** OCEAN MODEL PARAMETERS ********');
fprintf(fileID,'%s\n',['Start day = ',num2str(model.start_day)]);
fprintf(fileID,'%s\n',['Ocean model timestep = ',num2str(model.timestep),' seconds']);
fprintf(fileID,'%s\n',['Directory of ocean model outputs = ',model.path]);
fprintf(fileID,'%s\n',['Ocean model resolution used = every ',num2str(model.outputskip),' file(s)']);
fprintf(fileID,'%s\n',['Ocean model periodicities: E-W = ',num2str(model.periodic_ew),', N-S = ',num2str(model.periodic_ns)]);
fprintf(fileID,'%s\n',[]);
fprintf(fileID,'%s\n',[]);
fprintf(fileID,'%s\n','******** PARTICLE TRACKING PARAMETERS ********');
fprintf(fileID,'%s\n',['Particle initialization (1st seeding) = Day ',num2str(particle.initime)]);
fprintf(fileID,'%s\n',['Particle seeding frequency = every ',num2str(particle.inifreq),' day(s)']);
fprintf(fileID,'%s\n',['Timestep for particle tracking = ',num2str(particle.timestep),' day(s)']);
fprintf(fileID,'%s\n',['Length of particle tracking simulation = ',num2str(particle.length),' day(s)']);
fprintf(fileID,'%s\n',['Particle output frequency = ',num2str(particle.outfreq),' day(s)']);
fprintf(fileID,'%s\n',['Particle output format = ',particle.outputformat]);
fprintf(fileID,'%s\n',['Particle output directory = ',particle.outputdir]);
fprintf(fileID,'%s\n',['Particle output filename = ',particle.outputfilename]);
fprintf(fileID,'%s\n',[]);
fprintf(fileID,'%s\n',['Particle tracking direction: ',particle.direction]);
fprintf(fileID,'%s\n',[]);
fprintf(fileID,'%s\n','Particle seeding strategy in x:');
fprintf(fileID,'%s\n',['istart = ',num2str(particle.istart)]);
fprintf(fileID,'%s\n',['irange = ',num2str(particle.irange)]);
fprintf(fileID,'%s\n',['irez = ',num2str(particle.irez)]);
fprintf(fileID,'%s\n','Particle seeding strategy in y:');
fprintf(fileID,'%s\n',['jstart = ',num2str(particle.jstart)]);
fprintf(fileID,'%s\n',['jrange = ',num2str(particle.jrange)]);
fprintf(fileID,'%s\n',['jrez = ',num2str(particle.jrez)]);
fprintf(fileID,'%s\n','Particle seeding strategy in z:');
fprintf(fileID,'%s\n',['kstart = ',num2str(particle.kstart)]);
fprintf(fileID,'%s\n',['krange = ',num2str(particle.krange)]);
fprintf(fileID,'%s\n',['krez = ',num2str(particle.krez)]);
fprintf(fileID,'%s\n',['Number of particle classes = ',num2str(particle.numofclasses)]);

% Close the file
fclose(fileID);
clear fileID nowtime ans

disp('DONE')
disp('%%%%%%%%%%%%%%%%%')
disp(' ')