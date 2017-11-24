%% This routine interpolate the velocty field in time and space onto the
% particles location using a 4-D linear interpolant. The velocity field has
% been interpolated linearly in time already in "getvelocityfield.m", so
% only the hydrographic variables are interpolated in time as well as in
% space.

disp('INTERPOLATE VARIABLES ONTO PARTICLES... ')

parti.doy(:) = tt;

%% Velocity field

% If velocity field extracted from face velocity-fluxes (divergence free)

% F = griddedInterpolant(model.ufx,model.ufy,model.ufz,model.uf);
% parti.u = F(parti.x,parti.y,parti.z);
% F = griddedInterpolant(model.vfx,model.vfy,model.vfz,model.vf);
% parti.v = F(parti.x,parti.y,parti.z);
% F = griddedInterpolant(model.wfx,model.wfy,model.wfz,model.wf);
% parti.w = F(parti.x,parti.y,parti.z);

% If velocity field extracted from centered velocities (NOT divergence
% free)
F = griddedInterpolant(model.x2,model.y2,model.z2,model.doy2,model.u);
parti.u = F(parti.x,parti.y,parti.z,parti.doy);
F = griddedInterpolant(model.x2,model.y2,model.z2,model.doy2,model.v);
parti.v = F(parti.x,parti.y,parti.z,parti.doy);
F = griddedInterpolant(model.x2,model.y2,model.z2,model.doy2,model.w);
parti.w = F(parti.x,parti.y,parti.z,parti.doy);

%% Hydrography
F = griddedInterpolant(model.x2,model.y2,model.z2,model.doy2,model.T);
parti.t = F(parti.x,parti.y,parti.z,parti.doy);
F = griddedInterpolant(model.x2,model.y2,model.z2,model.doy2,model.S);
parti.s = F(parti.x,parti.y,parti.z,parti.doy);
F = griddedInterpolant(model.x2,model.y2,model.z2,model.doy2,model.RHO);
parti.rho = F(parti.x,parti.y,parti.z,parti.doy);

%% Potential-vorticity
F = griddedInterpolant(model.x2,model.y2,model.z2,model.doy2,model.pv);
parti.pv = F(parti.x,parti.y,parti.z,parti.doy);

%% Relative vorticity
F = griddedInterpolant(model.x2,model.y2,model.z2,model.doy2,model.vor);
parti.vor = F(parti.x,parti.y,parti.z,parti.doy);
clear F
parti.wtotal = parti.wsink+parti.w;

%%
disp('DONE')
disp('%%%%%%%%%%%%%%%%%')
disp(' ')

