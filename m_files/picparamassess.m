% picparamassess is used to assess parameters and whether they be
% suitable for picard simulations.
%
% HG 2018-11-01

disp('=================')
R = 3.0 % heliocentric distance [AU] (Jesper used 2.56 AU)
sf = 400; % squeeze factor that artificially increases B and decreases Q

% The resolutions
dt = 4e-7 % Imagine that we use this timestep
dx = 8    % Imagine that we use this grid cell size
Lx = 720  % Imagine that we use this max(|x|)


% solar wind parameters
nsw = 7e6*R.^(-2)
Tesw = 15e4*R^(-2/3) % This seems to be the way Jesper scaled it. I don't
                     % know why, but the result is not unreasonable.
kTesw = kB*Tesw/elc
Tpsw = 8e4*R^(-2/3)
kTpsw = kB*Tpsw/elc

% Comet
Q = sf^(-3)*2.59d28 *R^(-5.18)
u = (-55.5*R+771)*(1+0.171*exp(-(R-1.24)/0.13))
nui = 8.28d-7/R^2 % ionisation frequency, high EUV
%nui = 3.31d-7/R^2 % ionisation frequency, low EUV


n=200e6*(nui/1.26e-7)*(Q/3.11e18)*min(3,(20/dx)^2) % plasma density
                                                   % [m^-3] (scaling 
                                  % Jesper's result of 200 cm^-3). The
                                  % last factor which includes the min
                                  % function is due to ion production
                                  % capping at 2 comet radii. The capping
                                  % is not yet implemented, but I will
                                  % make it so. 
kTe = 15;    % electron temperature [eV]
B = 4.9616482316929589e-9*sf; % B field [T] (Jesper had 1.78e-9*400=7.12e-7)

fpe = sqrt(n*elc^2/me/eps0)/(2*pi)
fce = (elc*B/me)/(2*pi)

vte = sqrt(2*kTe*elc/me)

rge = vte/(2*pi*fce)

lD = sqrt(eps0*kTe/(elc*n))

% Reasonability
timestepsPerCyclotronPeriod = (1/fce)/dt

timestepsPerPlasmaPeriod = (1/fpe)/dt

CellsPerDebyeLength = lD/dx

CellsPerGyroRadius = rge/dx

% Computational cost
cells = (2*Lx/dx)^3 % Our number of cells
Jesper.dt=4e-6; % Jesper's dt
Jesper.dx=20; % Jesper's dx
Jesper.Lx=1800; % Jesper's Lx
Jesper.cells=(2*Jesper.Lx/Jesper.dx)^3; % Jesper's number of cells
Jesper.job.time = 88 % Time it took to run Jesper's job [h]

temporalfactor = Jesper.dt/dt

spatialfactor = cells/Jesper.cells

totaltime = temporalfactor*spatialfactor*Jesper.job.time
