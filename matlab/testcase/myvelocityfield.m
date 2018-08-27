%% This piece of code overwrites the velocity field extracted from the model 
% with a user-specified velocity field. This is intended to be used to
% validate tyhe particle-tracking code by imposing a known flow.

%% Face velocities (those matter for advection)

% create face area matrix to go from velocities to velocity fluxes
Dy = diff(model.yf);
Dy = repmat(Dy',size(model.uf,1),1,size(model.uf,3));
Dz = diff(model.zf);
Dz = repmat(Dz,1,size(model.uf,1),size(model.uf,2));
Dz = permute(Dz,[2 3 1]);

model.uf = ones(size(model.uf)).*Dy.*Dz;


Dx = diff(model.xf);
Dx = repmat(Dx,1,size(model.vf,2),size(model.vf,3));
Dz = diff(model.zf);
Dz = repmat(Dz,1,size(model.vf,1),size(model.vf,2));
Dz = permute(Dz,[2 3 1]);

model.vf = ones(size(model.vf)).*Dx.*Dz;

Dy = diff(model.yf);
Dy = repmat(Dy',size(model.wf,1),1,size(model.wf,3));
Dx = diff(model.xf);
Dx = repmat(Dx,1,size(model.wf,2),size(model.wf,3));
model.wf = 0*ones(size(model.wf)).*Dx.*Dy;

clear Dx Dy Dz

%% Center velocities (NOT used for advection, just recorded.)
model.u = ones(size(model.u));
model.v = ones(size(model.v));
model.w = ones(size(model.w));
