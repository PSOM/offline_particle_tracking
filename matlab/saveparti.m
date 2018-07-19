%% This routine saves the particle data either in a CSV file, or in an
% SQLite database.
%       If saved in a CSV file, one file is created for every recorded
%       timestep
%       If saved as sqlite, data is added to the database at every recorded
%       timestep

%% CSV file format
if strcmp(particle.outputformat,'csv')
    disp('SAVE PARTICLE DATA in CSV-file... ')
    
    % Define the variable to print in output file
    % exclude velocity field at previous timestep
    vartoprint = [1,2,6:16];
    
    % Define header
    thefields = fieldnames(parti);
    for ii = 1:length(vartoprint)
        if ii == 1
            HEADER = thefields{vartoprint(ii)};
        else
            HEADER = [HEADER,', '];
            HEADER = [HEADER,thefields{vartoprint(ii)}];
        end
    end; clear ii
    
    % Define filename
    FILENAME = fullfile(particle.outputdir,[particle.outputfilename,'_doy',num2str(tt,'%03.2f'),'.csv']);
    
    % test if file already exist
    if exist(FILENAME,'file')~=0
        error('File already exist - This code is set to never overwrite a file');
    end
    
    % Create file
    fileID = fopen(FILENAME,'w');
    
    % Prints Header
    fprintf(fileID,'%s\n',HEADER);
    fprintf(fileID,'%03.2f, %u, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f\n',...
        transpose(cat(2,parti.doy, parti.id, parti.x, parti.y, parti.z, parti.u, parti.v, ...
        parti.w, parti.wsink, parti.wtotal, parti.s, parti.t, parti.rho, parti.pv, parti.vor)));
    fclose(fileID);
    clear vartoprint thefields HEADER FILENAME fileID
    
    %% SQLite database
elseif strcmp(particle.outputformat,'sqlite')
    disp('SAVE PARTICLE DATA in SQLite... ')
    
    % Set database filename
    dbfile = [particle.outputfilename,'_',particle.direction,'.db'];
    % Set database table name
    tablename = 'particles';
    
    % Test if database file exists
    if exist(fullfile(particle.outputdir,dbfile),'file')==0
        % Creates the SQLite database if does not exist
        conn = sqlite(fullfile(particle.outputdir,dbfile),'create');
        
        % Creates the table in the database
        sqltable = ['CREATE TABLE ',tablename,' (DOY REAL, ID INTEGER, x REAL, y REAL, z REAL, u REAL, v REAL, w REAL, wsink REAL, wtotal REAL, salinity REAL, temperature REAL, density REAL, PV REAL, vorticity REAL);'];
        exec(conn,sqltable);
    else
        % test if file already exist at the beginning of the run
        if tt == particle.initime
            error('File already exist - This code is set to never overwrite a file');
        end
        % connects to database if file exists
        conn = sqlite(fullfile(particle.outputdir,dbfile),'connect');
        % Test connection
        if conn.IsOpen
            disp('Connected to SQLite database');
        else
            error('Failed connection to SQLite database');
        end
    end
    
    % Export particle data into the database
    colnames = {'DOY', 'ID','x','y','z','u','v','w','wsink','wtotal','salinity','temperature','density','PV','vorticity'};
    if isempty(parti.x)==0
        insert(conn,tablename,colnames,cat(2,...
            parti.doy,...
            parti.id,...
            parti.x,...
            parti.y,...
            parti.z,...
            parti.u,...
            parti.v,...
            parti.w,...
            parti.wsink,...
            parti.wtotal,...
            parti.s,...
            parti.t,...
            parti.rho,...
            parti.pv,...
            parti.vor))
    end
    % Close connection
    close(conn)
    disp('Disonnected to SQLite database');
    
    clear dbfile tablename conn sqltable colnames
else
    error('Cannot recognize particle output format');
end

disp('DONE')
disp('%%%%%%%%%%%%%%%%%')
disp(' ')