%GENERAL PARAMETERS
iter_end = 100;       % The iteration after which execution ends
startfromdumpfile='no';  % if 'yes' the distribution is read from the
                         % dumps directory
dump_period_dump=1000;   % number of iterations between dumps one can start from
dump_period_particles=100000; % number of iterations between particle dumps
dump_period_fields=500;  % number of iterations between field dumps
write_period_screen=100; % number of iterations between showing onscreen
                         % life signs.

Nspecies = 4;        % number of species included in the simulation
dt   =  4.0d-7;      % timestep
iprocs = 8;          % processes in x dimension, should be multiple of 2
jprocs = 8;          % processes in y dimension, should be multiple of 2
kprocs = 8;          % processes in z dimension, should be multiple of 2
Nx_local = 30;       % should be equal to INTEGER*8-2, e.g., 6, 30, 54, 78, 102
Ny_local = 30;       % should be equal to INTEGER*8-2, e.g., 6, 30, 54, 78, 102
Nz_local = 30;       % should be equal to INTEGER*8-2, e.g., 6, 30, 54, 78, 102

xmin = -960.0d0;     % physical domain
ymin = -960.0d0;     % physical domain
zmin = -960.0d0;     % physical domain
xmax =  960.0d0;     % physical domain
ymax =  960.0d0;     % physical domain
zmax =  960.0d0;     % physical domain

nucleusradius = 5.0e0;   % radius of the nucleus [m], 2.0d3/400.0d0
flatradius = 2.5e1;      % radius of flat density part [m], 1.0d4/400.0d0
Galandradius = 8.0e2;    % outer radius of Galand model
fadeoutradius = 8.8e2;   % radius at which the cometary ion density is zero

B0x = 0.0d-9;                % Constant magnetic field in the simulation region.
B0y = 1.9846592926771838d-6; % 5nT is in the range observed (Goetz et al. 2017)
B0z = 0.0d-9;                % Scaled thus: 400.0d0*4.9616482316929589d-9

Nprobes = 5;                 % Number of probes defined below
%END

%SPECIES 1 - electrons. Species 1 must be electrons.
ppc = 64; % particles per cell, will be recomputed for SW electrons
mass = 9.10938356d-31;     % mass [kg]
charge = -1.602176634d-19; % charge [C]
upstreamdensity = 7.78d5;  % upstream (SW) density [m^-3]. 
upstreamkelvin = 7.21d4;   % solar wind temperature for this species
v0x =-430d3;               % solar wind velocity x component [m/s]
v0y = 0.0d3;               % solar wind velocity y component [m/s]
v0z = 0.0d3;               % solar wind velocity z component [m/s]
vn  = 0.0d0;               % neutral gas radial velocity component
Qn  = 0.0d0;               % neutral gas production rate [s^-1]
Ek  = 0.0d0;               % excess energy in photo-ionisation [eV]
nu_i = 0.0d0;              % ionisation frequency [s^-1]
cometion='no';             % If 'yes', a Galand model density profile is created
                           % It should be 'no' for electrons.
productspecies = 0;        % species number for the electrons that quasi-
                           % neutralise this species. This must be 0
                           % for electrons
%END

%SPECIES 2 - solar wind protons. Species 2 must be solar wind protons.
ppc = 64; % particles per cell
mass = 1.672621898d-27;    % mass [kg]
charge = 1.602176634d-19;  % charge [C]
upstreamdensity = 7.78d5;  % upstream (SW) density [m^-3]
upstreamkelvin = 3.85d4;   % solar wind temperature for this species
v0x =-430d3;               % solar wind velocity x component [m/s]
v0y = 0.0d3;               % solar wind velocity y component [m/s]
v0z = 0.0d3;               % solar wind velocity z component [m/s]
vn  = 0.0d0;               % neutral gas radial velocity component
Qn  = 0.0d0;               % neutral gas production rate [s^-1]
Ek  = 0.0d0;               % excess energy in photo-ionisation [eV]
nu_i = 0.0d0;              % ionisation frequency [s^-1]
cometion='no';             % If 'yes', a Galand model density profile is created
                           % It should be 'no' for electrons.
productspecies = 1;        % species number for the electrons that quasi-
                           % neutralise this species. This must be 0
                           % for electrons
%END

%SPECIES 3 - cometary H2O+. 
ppc = 16;                  % particles per cell relative the upstream
                           % electron (species 1) value
mass = 2.9914162400592d-26;% mass [kg]
charge = 1.602176634d-19;  % charge [C]
upstreamdensity = 0.0d0;   % upstream (SW) density [m^-3]
upstreamkelvin = 0.0d0;    % solar wind temperature for this species
v0x = 0.0d3;               % solar wind velocity x component [m/s]
v0y = 0.0d3;               % solar wind velocity y component [m/s]
v0z = 0.0d3;               % solar wind velocity z component [m/s]
vn  = 6.045d2;             % neutral gas radial velocity component
Qn  = 2.1865d23;           % neutral gas production rate
                           % (here, 8.746e25/400) [s^-1]
Ek  = 15.0d0;              % excess energy in photo-ionisation [eV]
nu_i = 3.68d-8;            % ionisation frequency (here, 3.31d-7 at 1 AU,
                           % low UV) [s^-1] 
cometion='yes';            % If 'yes', a Galand model density profile is created
                           % It should be 'no' for electrons.
productspecies = 4;        % Species number for the electron produced
                           % when ionising this species. This must be 0
                           % for electrons
%END

%SPECIES 4 - electrons. These shall be electrons created in ionisations
ppc = 16;                  % particles per cell relative the upstream
                           % electron (species 1) value
mass = 9.10938356d-31;     % mass [kg]
charge = -1.602176634d-19; % charge [C]
upstreamdensity = 0.0d0;   % will be recomputed for SW electrons [m^-3]
upstreamkelvin = 0.0d0;    % solar wind temperature for this species
v0x = 0.0d3;               % solar wind velocity x component [m/s]
v0y = 0.0d3;               % solar wind velocity y component [m/s]
v0z = 0.0d3;               % solar wind velocity z component [m/s]
vn  = 0.0d0;               % neutral gas radial velocity component
Qn  = 0.0d0;               % neutral gas production rate [s^-1]
Ek  = 0.0d0;               % excess energy in photo-ionisation [eV]
nu_i = 0.0d0;              % ionisation frequency [s^-1]
cometion='no';             % If 'yes', a Galand model density profile is created
                           % It should be 'no' for electrons.
productspecies = 0;        % Species number for the electron produced
                           % when ionising this species. This must be 0
                           % for electrons
%END

% Here follow the probe definitions. Particles are saved at intervals
% given by dump_period_pprobes set in the general parameters section. The
% electric field is saved every timestep. Try to put (xc,yc,zc) at the 
% centre of a grid cell. The program finds the nearest centre, but if you
% have (xc,yc,zc) at a cell boundary, it is not clear where the nearest
% centre is.

%PROBE 1
xc = 4.0d0;                % x coordinate of probe
yc = 4.0d0;                % y coordinate of probe
zc = 4.0d0;                % z coordinate of probe
rprobe = 4.0e0;            % Probe radius for particle saves.
                           % If rprobe<=0, no particles are saved.
%END

%PROBE 2
xc =-4.0d0;                % x coordinate of probe
yc =-4.0d0;                % y coordinate of probe
zc =-4.0d0;                % z coordinate of probe
rprobe =-1.0e0;            % Probe radius for particle saves.
                           % If rprobe<=0, no particles are saved.
%END

%PROBE 3
xc = 4.0d0;                % x coordinate of probe
yc = 12.0d0;               % y coordinate of probe
zc = 4.0d0;                % z coordinate of probe
rprobe = 4.0e0;            % Probe radius for particle saves.
                           % If rprobe<=0, no particles are saved.
%END

%PROBE 4
xc = 100.0d0;              % x coordinate of probe
yc =-4.0d0;                % y coordinate of probe
zc = 4.0d0;                % z coordinate of probe
rprobe = 4.0e0;            % Probe radius for particle saves.
                           % If rprobe<=0, no particles are saved.
%END

%PROBE 5
xc = -100.0d0;             % x coordinate of probe
yc =-4.0d0;                % y coordinate of probe
zc = 4.0d0;                % z coordinate of probe
rprobe = 0.0e0;            % Probe radius for particle saves.
                           % If rprobe<=0, no particles are saved.
%END
