%% This routine identifies the model output files surrounding the particle
% tracking timestep to extract the velocity fields. It then linearly
% interpolates in time the velocity field to get the velocity field for the
% timestamp considered

disp('IMPORTING VELOCITY FIELDS...')

% List of model outputs present in the folder path
files = dir([model.path,'face_*.cdf']);
files2 = dir([model.path,'full_*.cdf']);

% Delete some outputfile filenames to match the desired model output
% temporal resolution to be used in the tracking
files = files(1:model.outputskip:length(files));
files2 = files2(1:model.outputskip:length(files2));

% Rename the considered day of year
doy = tt;

% find the two files to be imported
timestepind = [];
for ii = 1:length(files)-1
    if str2double(files(ii).name(end-8:end-4))/(86400/model.timestep) + model.start_day <= doy &&...
            str2double(files(ii+1).name(end-8:end-4))/(86400/model.timestep) + model.start_day >= doy
        timestepind = [ii, ii+1];
        break
    end
end; clear ii

if isempty(timestepind)
    error('No output files were found')
end

% reset variables
model.uf = [];
model.vf = [];
model.wf = [];
model.x = [];
model.y = [];
model.z = [];
model.RHO = [];
model.T = [];
model.S = [];
model.doy = [];

for ii = 1:length(timestepind)
    
    % Define filename
    filename = [files(timestepind(ii)).folder,'/',files(timestepind(ii)).name];
    
    % Import velocity fluxes
    model.uf(:,:,:,ii) = ncread(filename,'uf')*1e3*1e5;
    model.vf(:,:,:,ii) = ncread(filename,'vf')*1e3*1e5;
    model.wf(:,:,:,ii) = ncread(filename,'wf')*1e5*1e5/1000;
    
    % extract coordinate from full_* files
    filename2 = [files2(timestepind(ii)).folder,'/',files2(timestepind(ii)).name];
    model.x = ncread(filename2,'xc')*1000;
    model.y = ncread(filename2,'yc')*1000;
    model.z = ncread(filename2,'zc')*1000;
    model.z = squeeze(model.z(1,1,:));
    
    % Import variables
    model.RHO(:,:,:,ii) = ncread(filename2,'rho');
    model.T(:,:,:,ii) = ncread(filename2,'temp');
    model.S(:,:,:,ii) = ncread(filename2,'s');
    
    % Extract DOY from file name
    model.doy(ii) = str2double(filename(end-8:end-4))*432/86400 + model.start_day;
    
end; clear ii filename files timestepind doy files2 filename2

% Interpolate linearly in time the velocity field
if size(model.uf,4)>1
    model.uf = model.uf(:,:,:,1)+((tt-model.doy(1))/(model.doy(2)-model.doy(1))*(model.uf(:,:,:,2)-model.uf(:,:,:,1)));
    model.vf = model.vf(:,:,:,1)+((tt-model.doy(1))/(model.doy(2)-model.doy(1))*(model.vf(:,:,:,2)-model.vf(:,:,:,1)));
    model.wf = model.wf(:,:,:,1)+((tt-model.doy(1))/(model.doy(2)-model.doy(1))*(model.wf(:,:,:,2)-model.wf(:,:,:,1)));
end

% define coordinates (not in PSOM output files yet...)
model.xf = model.x(1:end-1)+ 0.5*(model.x(2:end)-model.x(1:end-1));
model.yf = model.y(1:end-1)+ 0.5*(model.y(2:end)-model.y(1:end-1));

% Import zf from zgrid.out
filename = [model.path,'zgrid.out'];
fileID = fopen(filename,'r');
dataArray = textscan(fileID, '%*s%f%*s%[^\n\r]', 'Delimiter', ' ', 'MultipleDelimsAsOne', true, 'TextType', 'string', 'HeaderLines' ,36, 'ReturnOnError', false, 'EndOfLine', '\r\n');
fclose(fileID);
model.zf = dataArray{1:end-1};
clearvars filename dataArray ans fileID;
model.zf(1) = [];
model.zf(end) = [];

% Grid the 3D coordinates
[model.x2,model.y2,model.z2,model.doy2] = ndgrid(model.x,model.y,model.z,model.doy);
[model.ufx,model.ufy,model.ufz] = ndgrid(model.xf,model.y(2:end-1),model.z(2:end-1));
[model.vfx,model.vfy,model.vfz] = ndgrid(model.x(2:end-1),model.yf,model.z(2:end-1));
[model.wfx,model.wfy,model.wfz] = ndgrid(model.x(2:end-1),model.y(2:end-1),model.zf);

disp('DONE')
disp('%%%%%%%%%%%%%%%%%')
disp(' ')

