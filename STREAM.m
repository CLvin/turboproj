%% TRACING THERMODYNAMIC VARIABLES STATION-BY-STATION
constants;

NLE = 21;
NTE = 30;
H= (RSHROUD-RHUB)/11;



%% UPDATING THE DENSITY ON THE INLET

for J=1:NSTRM
    CSQ = (RCU(J,1)/RADIUS(J,1))^2 + CZ(J,1)^2 + CR(J,1)^2;
    HSTATIC = HTOTAL(J,1) - CSQ/2;
    PSTATIC = PTOTAL(J,1)*(HSTATIC/HTOTAL(J,1))^3.5;
    DENSITY(J,1) = PSTATIC/(RGAS*HSTATIC/CP);
end

%% SWEEPING THE PLANES FROM PLANE 2 TO NSTAN
ERRDENS = 0;
CHANGE = 0;
for I=2:NSTATN
    LEFT = I - 1;
    if (I > NLE && I < NTE)
        LEFT = NLE;
    end
    NSTART = 1;
    PSIDN = PSI(NSTART,LEFT);
    PSIUP = PSI(NSTART+1,LEFT);
    
    for J=1:NSTRM
        DESIRED = PSI(J, I);
        if (~(DESIRED<=PSIUP) && ~(DESIRED>=PSIDN))
            NSTART = NSTART + 1;
            PSIDN = PSI(NSTART,LEFT);
            PSIUP = PSI(NSTART+1,LEFT);
        else
            DELTA = (DESIRED - PSIDN)/(PSIUP - PSIDN); %Streamline origin 
                                                       %has been located
        end
        ROTATE = 0;
        if (I>NLE && I<NTE)
            ROTATE = OMEGA;
        end
        % QUANTITIES AT "1" NEEDED FOR CONSERVATION PRINCIPLE
        RCU1 = DELTA*(RCU(NSTART+1,LEFT) - RCU(NSTART,LEFT)) + RCU(NSTART,LEFT);
        HTOTAL1 = DELTA*(HTOTAL(NSTART+1,LEFT) - HTOTAL(NSTART,LEFT)) + HTOTAL(NSTART,LEFT);
        PTOTAL1 = DELTA*(PTOTAL(NSTART+1,LEFT) - PTOTAL(NSTART,LEFT)) + PTOTAL(NSTART,LEFT);    
        RAD1 = DELTA * (RADIUS(NSTART+1,LEFT) - RADIUS(NSTART,LEFT)) + RADIUS(NSTART,LEFT);
        % QUANTITIES AT "1" NEEDED FOR LOSS CALCULATION
        DENS1 = DELTA*(DENSITY(NSTART+1,LEFT) - DENSITY(NSTART,LEFT)) + DENSITY(NSTART,LEFT);
        CZ1 = DELTA*(CZ(NSTART+1,LEFT) - CZ(NSTART,LEFT)) + CZ(NSTART,LEFT);
        CR1 = DELTA*(CR(NSTART+1,LEFT) - CR(NSTART,LEFT)) + CR(NSTART,LEFT);
        ENTROP1 = DELTA*(ENTROPY(NSTART+1,LEFT) - ENTROPY(NSTART,LEFT)) + ENTROPY(NSTART,LEFT);
        C1SQ = CZ1^2 + CR1^2 + (RCU1/RAD1)^2;
        HSTAT1 = HTOTAL1 - C1SQ/2;
        PSTAT1 = PTOTAL1*(HSTAT1/HTOTAL1)^3.5;
        
        %ROTATING AND NON-ROTATING QUANTITIES AT REFERENCE STATION
        
        ROTALP1 = HTOTAL1 - ROTATE * RCU1;
        HOR2 = ROTALP1 + (ROTATE * RADIUS(J, I))^2/2;
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
        
        
        RCU(J, I) = RCU1;
        if (I == NLE + 1)
            RCU(J, I) = 117.85;
        end
        if (I == NLE + 2)
            RCU(J, I) = 157.1;
        end
        if (I == NLE + 3)
            RCU(J, I) = 196.35;
        end
        if (I == NTE)
            RCU(J, I) = 235.6;
        end
        
        HTOTAL(J, I) = HTOTAL1 + ROTATE * (RCU(J, I) - RCU1);
        PTOTAL(J, I) = POR2 * (HTOTAL(J, I)/HOR2).^(GAMAM);
        
                
        %Calculation Block for Common Vars.
        CU = RCU(J, I)/RADIUS(J, I);
        VU = CU - ROTATE*RADIUS(J, I);
        V2SQ = VU^2 + CZ(J, I)^2 + CR(J, I)^2;
        C2SQ = CU^2 + CZ(J, I)^2 + CR(J, I)^2;
        HSTATIC = HOR2 - V2SQ/2;
        HTOTAL(J, I) = HSTATIC + C2SQ/2;
        PTOTAL(J, I) = POR2 * (HTOTAL(J, I)/HOR2)^GAMAM;
        DENSOLD = DENSITY(J, I);
        DENSITY(J, I) = PSTATIC/(RGAS*HSTATIC/CP);
        CHANGE = abs(DENSOLD - DENSITY(J, I));
        ERRDENS = max(CHANGE, ERRDENS);
        ENTROPY(J, I) = CP * log(HTOTAL(J, I)/HTOTAL1) - RGAS*log(PTOTAL(J, I)/PTOTAL1) + ENTROP1;
        %Assemble RHS for Unknown
        ERRRHS = 0;
        CHANGE = 0;
        for I = 1:NSTATN
            for J = 2: NSTRM-1
                RHSOLD = RHS(J, I);
                CU = RCU(J, I)/RADIUS(J, I);
                TSTATIC = (HTOTAL(J, I)-(CZ(J, I)^2+CR(J, I)^2 + CU^2)/2)/CP;
                RHS(J, I) = -1/CZ(J, I) * 2*pi/XMASS*(CU/RADIUS(J, I)*RCU(J+1, I) - RCU(J-1, I))/(2*H) + TSTATIC*(ENTROPY(J+1,I)-ENTROPY(J-1,I))/(2*H);
                CHANGE = abs(RHS(J, I) - RHSOLD);
                ERRRHS = max(ERRRHS, CHANGE);
            end
        end
    end
    
end
%% Solving the equation

%[X, ITER] = GaussSolver(A,RHS);

        