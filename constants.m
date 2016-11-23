CZ = 0; % Axial Velocity
CR = 0; % Radial Velocity
CM = 0; % Meridional Vel.
CU = 0; % Tangential Vel.
RCU = 0; %Swirl
VU = 0; %Relative Tangential Vel.
C = 0; %Absolute Velocity
V = 0; %Relative Velocity

USPEED = 0; %Rot. Speed
ALPHA = 0; %Abs. angle in rad
OMEGA = 0; %Wheel Speed
BETA = 0; %Rel. angle
Z = 0; %Z-coord.
RADIUS = 0; %R-coord.

TSTATIC = 0; %Static Temp.
TTOTAL = 288; %Total Temp.
TOREL = 0; %Total Rel. Temp.
HOREL = 0; %Total Rel. Enth.
HTOTAL = 0; %Stag. Enth.
DENSITY = 0; %Density
ENTROPY = 0; %Entropy
ROTHALPY = 0; %Rothalpy
PSTATIC = 0; %Static Pressure
PTOTAL = 0; %Total Pressure
PSI = 0; %Stream function

CP = 1005; %Cp
GAMMA = 1.4; %Gamma
GAMAM = 3.5;
RGAS = 287.058; %Gas const.
XMASS = 0; %Mass
RHS = 0; %Right hand side

RHUB = 0.45; %radius of hub
RSHROUD = 0.50; %radius of shroud
G = 9.81; %gravity m/s^2
LOSS_COEFF = 0; %omeglos

NSTATN = 10; %number of vertical grid lines
NSTREAM = 10; %number of horizontal grid lines
NSTRM = 10;

PSI = zeros(NSTATN, NSTREAM);
RHS = zeros(NSTATN, NSTREAM);
CZ = zeros(NSTATN, NSTREAM);
CR = zeros(NSTATN, NSTREAM);
RADIUS = zeros(NSTATN, NSTREAM);
HTOTAL = zeros(NSTATN, NSTREAM);
RCU = zeros(NSTATN, NSTREAM);
PTOTAL = zeros(NSTATN, NSTREAM);
DENSITY = zeros(NSTATN, NSTREAM);
ENTROPY = zeros(NSTATN, NSTREAM);
BETA = zeros(NSTATN, NSTREAM);

