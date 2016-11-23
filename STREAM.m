%% TRACING THERMODYNAMIC VARIABLES STATION-BY-STATION
function STREAM(PSI,CZ,CR,RCU,DENSITY,RADIUS,HTOTAL,KOUNT,TOTAL,RHS,NSTRM,NSTATN,ERRRHS,ERRDENS,OMEGLOS,ENTROPY)
constants;
NLE = ???;
NTE = ???;
H= ???;
%% UPDATING THE DENSITY ON THE INLET

for j=1:NSTRM
    CSQ = (RCU(1,J)/RADIUS(1,J))^2 + CZ(1,J)^2 + CR(1,J)^2;
    HSTATIC = HTOTAL(1,J) - CSQ/2;
    PSTATIC = PTOTAL(1,J)*(HSTATIC/HTOTAL(1,J))^3.5;
    DENSITY(1,J) = PSTATIC/(RGAS*HSTATIC/CP);
end

%% SWEEPING THE PLANES FROM PLANE 2 TO NSTAN

for I=2:NSTATN
    LEFT = I - 1;
    if (I > NLE && I < NTE)
        LEFT = NLE;
    end
    NSTART = 1;
    PSIDN = PSI(LEFT,NSTART);
    PSIUP = PSI(LEFT,NSTART+1);
    
    for J=1:NSTRM
        DESIRED = PSI(I,J);
        if (~(DESIRED<=PSIUP) && ~(DESIRED>=PSIDN))
            NSTART = NSTART + 1;
            PSIDN = PSI(LEFT,NSTART);
            PSIUP = PSI(LEFT,NSTART+1);
        else
            DELTA = (DESIRED - PSIDN)/(PSIUP - PSIDN); %Streamline origin 
                                                       %has been located
        end
        ROTATE = 0;
        if (I>NLE && I<NTE)
            ROTATE = OMEGA;
        end
        % QUANTITIES AT "1" NEEDED FOR CONSERVATION PRINCIPLE
        RCU1 = DELTA*(RCU(LEFT,NSTART+1) - RCU(LEFT,NSTART)) + RCU(LEFT,NSTART);
        HTOTAL1 = DELTA*(HTOTAL(LEFT,NSTART+1) - HTOTAL(LEFT,NSTART)) + HTOTAL(LEFT,NSTART);
        PTOTAL1 = DELTA*(PTOTAL(LEFT,NSTART+1) - PTOTAL(LEFT,NSTART)) + PTOTAL(LEFT,NSTART);    
        RAD1 = DELTA * (RADIUS(LEFT,NSTART+1) - RADIUS(LEFT,NSTART)) + RADIUS(LEFT,NSTART);
        % QUANTITIES AT "1" NEEDED FOR LOSS CALCULATION
        DENS1 = DELTA*(DENSITY(LEFT,NSTART+1) - DENSITY(LEFT,NSTART)) + DENSITY(LEFT,NSTART);
        CZ1 = DELTA*(CZ(LEFT,NSTART+1) - CZ(LEFT,NSTART)) + CZ(LEFT,NSTART);
        CR1 = DELTA*(CR(LEFT,NSTART+1) - CR(LEFT,NSTART)) + CR(LEFT,NSTART);
        ENTROP1 = DELTA*(ENTROPY(LEFT,NSTART+1) - ENTROPY(LEFT,NSTART)) + ENTROPY(LEFT,NSTART);
        C1SQ = CZ1^2 + CR1^2 + (RCU1/RAD1)^2;
        HSTAT1 = HTOTAL - C1SQ/2;
        PSTAT1 = PTOTAL1*(HSTAT/HTOTAL1)^3.5;
        
        %ROTATING AND NON-ROTATING QUANTITIES AT REFERENCE STATION
        
        ROTALP1 = HTOTAL1 - ROTATE * RCU1;
        HOR2 = ROTALP1 + (ROTATE * RADIUS(I,J))^2/2;
        HOR1 = ROTALP1 + (ROTATE * RAD1)^2/2;
        POR1 = PTOTAL1 * (HOR1/HTOTAL1)^GAMAM;
        POR2IDL = POR1 *(HOR2/HOR1)^GAMAM;
        if ((I > NLE) && (I < NTE))
            OMEGLOS = 0.03*(I - NLE)/(NTE-NLE);
            PLOSS = OMEGLOS*(POR1 - PSTAT1);
        else
            OMEGLOS = 0;
            PLOSS = 0;
        end
        POR2 = POR2IDL - PLOSS;
        
        %For Project 1,2
        RCU(I,J) = RCU1;
        if (I == NLE + 1)
            RCU(I,J) = 117.85;
        end
        if (I == NLE + 2)
            RCU(I,J) = 157.1;
        end
        if (I == NLE + 3)
            RCU(I,J) = 196.35;
        end
        if (I == NTE)
            RCU(I,J) = 235.6;
        end
        
        HTOTAL(I,J) = HTOTAL1 + ROTATE * (RCU(I,J) - RCU1);
        PTOTAL(I,J) = POR2 * (HTOTAL(I,J)/HOR2).^(GAMAM);
        
        %For Project 3.
        %PTOTAL(I,J) = PTOTAL1 * PressureRatio;
        %HTOTAL(I,J) = HOR2 * (PTOTAL(I,J)/POR2).^GAMAM;
        %RCU(I,J) = RCU1 + (HTOTAL(I,J) - HTOTAL1)/ ROTATE;
        
        %Calculation Block for Common Vars.
        CU = RCU(I,J)/RADIUS(I,J);
        VU = CU - ROTATE*RADIUS(I,J);
        V2SQ = VU^2 + CZ(I,J)^2 + CR(I,J)^2;
        C2SQ = CU^2 + CZ(I,J)^2 + CR(I,J)^2;
        HSTATIC = HOR2 - V2SQ/2;
        HTOTAL(I,J) = HSTATIC + C2SQ/2;
        PTOTAL(I,J) = POR2 * (HTOTAL(I,J)/HOR2)^GAMAM;
        DENSOLD = DENSITY(I,J);
        DENSITY(I,J) = PSTATIC/(RGAS*HSTATIC/CP);
        CHANGE = abs(DENSOLD - DENSITY(I,J));
        ERRDENS = max(CHANGE, ERRDENS);
        ENTROPY(I,J) = CP * log(HTOTAL(I,J)/HTOTAL1) - RGAS*log(PTOTAL(I,J)/PTOTAL1) + ENTROP1;
        %Assemble RHS for Unknown
        ERRRHS = 0;
        CHANGE = 0;
        for i = 1:NSTATN
            for j = 2: NSTRM-1
                RHSOLD = RHS(i,j);
                CU = RCU(i,j)/RADIUS(i,j);
                TSTATIC = (HTOTAL(i,j)-(CZ(i,j)^2+CR(i,j)^2 + CU^2)/2)/CP;
                RHS(i,j) = -1/CZ(i,j) * 2*pi/XMASS*(CU/RADIUS(i,j)*RCU(i,j+1) - RCU(i,j-1))/(2*H) + TSTATIC*(ENTROPY(i,j+1)-ENTROPY(i,j-1))/(2*H);
                CHANGE = abs(RHS(i,j) - RHSOLD);
                ERRRHS = max(ERRRHS, CHANGE);
            end
        end
    end
    
end
end
        