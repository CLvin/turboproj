function [] = gaussbns(A, X, RHS, MBAND, MM)
%**************************************************************** 
%PROF. W.G. HABASHI, McGILL UNIVERSITY  
%PROGRAM TO SOLVE N SIMULTANEOUS EQUATIONS BY GAUSS-SEIDEL
%FOR BANDED NON-SYMMETRIC MATRICIES
%
%Ported to MATLAB by Calvin Liao, Nov. 21, 2016
%
%A = matrix of coefficients (N*2*(MBAND-1)+1)
%RHS = right hand side
%X = vector of unknowns, initial guess made when function called
%TOLER = max acceptable tol.
%NTRY = number of iterations
%ERROR = half band of matrix, including diagonal
%MM = width of matrix, (2*(MBAND-1)+1)
%
%****************************************************************
[N, MM] = size(A);
TOLER = 1.0E-05;
NTRY = 50;
A = (2.*N).*(MBAND-1) + 1;

for i = 1:NTRY
    ERROR = 0.0;
    
    %row sweep
    for j = 1:N
        XOLD = X(j);
        SUM = RHS(j);
    end
    %column sweeep
    JSTART = MAX(1,i-
end
end
