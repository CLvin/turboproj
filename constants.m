CZ = 0; % Axial Velocity
CR = 0; % Radial Velocity
CM = 0; % Meridional Vel.
CU = 0; % Tangential Vel.
VU = 0; %Relative Tangential Vel.
C = 0; %Absolute Velocity
V = 0; %Relative Velocity

USPEED = 6000; %Rot. Speed
ALPHA = 0; %Abs. angle in rad
OMEGA = (2*pi*USPEED)/60; %Wheel Speed
BETA = 0; %Rel. angle
Z = 0; %Z-coord.

TSTATIC = 0; %Static Temp.
TTOTAL = 288; %Total Temp.
TOREL = 0; %Total Rel. Temp.
HOREL = 0; %Total Rel. Enth.
HTOTAL = 0; %Stag. Enth.
DENSITY = 1.5; %Density
ENTROPY = 0; %Entropy
ROTHALPY = 0; %Rothalpy
PSTATIC = 0; %Static Pressure
PTOTAL = 0; %Total Pressure
PSI = 0; %Stream function

CP = 1005; %Cp
GAMMA = 1.4; %Gamma
GAMAM = 3.5;
RGAS = 287.058; %Gas const.
XMASS = 30.4; %Mass
RHS = 0; %Right hand side

RHUB = 0.45; %radius of hub
RSHROUD = 0.50; %radius of shroud
G = 9.81; %gravity m/s^2
LOSS_COEFF = 0.03; %omeglos


NSTATN = 51; %number of vertical grid lines
NSTREAM = 11; %number of horizontal grid lines
NSTRM = 11;


PSI = zeros(NSTREAM, NSTATN);
RHS = zeros(NSTREAM, NSTATN);
%% Calculating radial velocity at inlet
CR = zeros(NSTREAM, NSTATN);

%% Calculating axial velocity at inlet
CZ = zeros(NSTREAM, NSTATN);
CZ(:,1) = XMASS/(DENSITY*pi*(RSHROUD^2-RHUB^2));
%% Caluculating radius at the grid points
RADIUS = zeros(11,51);
for i=1:11
    for j=1:51
        RADIUS(i,j) = RSHROUD - (i-1)*(RSHROUD - RHUB)/10;
    end
end
%% HTOTAL at inlet
HTOTAL = zeros(NSTREAM, NSTATN);
for i=1:11
    for j=1:51
        HTOTAL(i,1)=OMEGA*RADIUS(i,j)*(117.8-39.3);
    end
end
%% RCU values at inlet
RCU = zeros(NSTREAM, NSTATN);
RCU(:, 1:21) = 39.3; %Swirl Before rotor
%% PTOTAL at inlet
PTOTAL = zeros(NSTREAM, NSTATN);


DENSITY = zeros(NSTREAM, NSTATN);
ENTROPY = zeros(NSTREAM, NSTATN);
BETA = zeros(NSTREAM, NSTATN);

