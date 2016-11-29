function [X, ITER] = GaussSolver(A,RHS)
n = length(RHS);
X = zeros(n,1);
TOLER = ones(n,1);
iteration = 0;
while max(TOLER) > 1*10^-5
    iteration = iteration + 1;
    Z = X;  % save current values to calculate error later
    for i = 1:n
        j = 1:n; % define an array of the coefficients' elements
        j(i) = [];  % eliminate the unknow's coefficient from the remaining coefficients
        Xtemp = X;  % copy the unknows to a new variable
        Xtemp(i) = [];  % eliminate the unknown under question from the set of values
        X(i) = (RHS(i) - sum(A(i,j) * Xtemp)) / A(i,i);
    end;
    Xsolution(:,iteration) = X;
    TOLER = sqrt((X - Z).^2);
end
[~, ITER] = size(Xsolution);
end
