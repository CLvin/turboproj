RHUB = 0.45; %radius of hub
RSHROUD = 0.50; %radius of shroud
RADIUS = zeros(11,51); %linearly interpolating the radius
H= (RSHROUD-RHUB)/10; %Distance between stations
NLE = 21;
NTE = 30;

for i=1:11
    for j=1:51
        RADIUS(i,j) = RSHROUD - (i-1)*(RSHROUD - RHUB)/10;
    end
end

DENSITY(1:11, 1:51) = 1.5; %Setting the initial density

PSI = zeros(11,51); %Setting initial value of PSI
for i=1:11
    for j=1:51
        PSI(i,j) = (RADIUS(i,j)^2-RHUB^2)/(RSHROUD^2 - RHUB^2);
    end
end
%PSI(11,:)= 0;
%PSI(1,:) = 1;

CZ = zeros(11,51); %Axial velocity
CR = zeros(11,51); %Radial velocity
XMASS = 30.4; %Mass flow rate
for i=2:10 %calculating Axial velocity and Radial velocity
    for j=2:50
        CZ(i,j) = (XMASS/(2*pi*DENSITY(i,j)*RADIUS(i,j)))*((PSI(i,j+1)-PSI(i,j-1))/(2*H));
        CR(i,j) = -(XMASS/(2*pi*DENSITY(i,j)*RADIUS(i,j)))*((PSI(i+1,j)-PSI(i-1,j))/(2*H));
    end
end

USPEED = 6000;
OMEGA = 2*pi*USPEED/60;
RCU = zeros(11,51);
RCU(:,1:21) = 39.3;
RCU(:,30:51) = 117.8;
TTOTAL(1:11,1) = 288;
CP = 1005;
HTOTAL = zeros(11,51);

for i = 1:11;
    C = CZ(i,1)^2 +CR(i,1)^2 +(RCU(i,1)/RADIUS(i,1))^2;
    TSTATIC1 = TTOTAL(i,1) - C^2/(2*CP);
    HSTATIC1 = CP * TSTATIC1;
    HTOTAL(i,1:21) = HSTATIC1 + C^2/2;
end


