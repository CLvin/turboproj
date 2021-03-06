function [err, tol, iter, xnew] = gaussbns(A, X, RHS)
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
[lower , upper] = bandwidth(A);
q = max(lower,upper);
MBAND = q+1;
[N , MM] = size(A);
TOLER = 1.0E-05;
NTRY = 50;
[L, U] = lu(A);

for ITERAT = 1:NTRY
  ERROR = 0.0;
        %row sweep
        for I = 2:N
            XOLD = X(I);
            SUM = RHS(I);
        %column sweeep
        JSTART = max(1, I - MBAND + 1)
        JFINIS = min(N, I + MBAND - 1)
            for J = JSTART:JFINIS
                II = I
                J
                JJ = J-I+1
                SUM = SUM - (A(II,JJ)*X(J));
            end
        %newly computed values at j'th iteration
        SUM = SUM + A(I, MBAND)*X(I);
        X(I) = SUM / A(I, MBAND);
        ERROR = max([ERROR ,abs(XOLD-X(I))]);
        end
  %check if max error is within tol.
end
  if (ERROR < TOLER)
      err = ERROR;
      tol = TOLER;
      iter = I;
      xnew = X;
  else
      fprintf('Iteration fails to converge');
  end
end
