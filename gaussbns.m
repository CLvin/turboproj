function [err, tol, iter, xnew] = gaussbns(A, X, RHS, MBAND, MM)
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
%or, MBAND = ((MM - 1)/2)+1
%****************************************************************
[N, MM] = size(A);
MBAND2 = ((MM-1)/2) + 1;
TOLER = 1.0E-05;
NTRY = 50;
A = (2.*N).*(MBAND2-1) + 1;

for i = 1:NTRY
    ERROR = 0.0;
    
    %row sweep
    for j = 1:N
        XOLD = X(j);
        SUM = RHS(j);
    %column sweeep
    JSTART = MAX(1,i - MBAND2 + 1);
    JFINIS = MIN(N, i + MBAND2 - 1);
        for k = JSTART:JFINIS
            II = i;
            JJ = j-i+1;
            SUM = SUM - A(II,JJ)*X(j);
        end
    %newly computed values at j'th iteration
    SUM = SUM + A(i, MBAND2)*X(i);
    X(i) = SUM / A(i, MBAND2);
    ERROR = max([ERROR (abs(XOLD-X(i)))]);
    end
    %check if max error is within tol.
    if (ERROR < TOLER)
        err = ERROR;
        tol = TOLER;
        iter = i;
        xnew = X;
        break
    else
        err = ERROR;
        tol = TOLER;
        xnew = X;
        break
    end  
end
end
