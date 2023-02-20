function dN = computeQuadraticBasisFunctionsAndFirstDerivatives(xi)
%% COMPUTEQUADRATICBASISFUNCTIONSANDFIRSTDERIVATIVES Returns the quadratic basis functions and their derivatives at a given parametric location
%
% Returns the quadratic shape functions corresponding corresponding to an
% one-dimensional element with parametric space [-1, +1]
%
%    Input :
%       xi : The coordinate on the parametric domain [-1, +1]
%
%   Output :
%       dN : 2x2 array for which dN(:,1) contains the basis functions
%            themselves, dN(:,2) the derivatives of the basis functions 
%
%% Function Implementation

%% 0. Input validation
arguments
    xi (1, 1) double {mustBeInRange(xi, -1, +1)}
end

%% 1. Computation of the basis functions' values
dN = [(xi*(xi - 1))/2, xi - 1/2
      (xi*(xi + 1))/2, xi + 1/2
      1 - xi^2       , -2*xi];

end