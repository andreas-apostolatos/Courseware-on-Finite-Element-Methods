function dN = computeLinearBasisFunctionsAndFirstDerivatives(xi)
%% COMPUTELINEARBASISFUNCTIONSANDDERIVATIVES Returns the linear basis functions and their derivatives for an one-dimensional element
%
% Returns the linear (hat) shape functions and their parametric derivatives 
% corresponding to an one-dimensional element
%
%    Input :
%       xi : The coordinate on the parametric domain [-1, 1]
%
%   Output :
%       dN : 2x2 array for which dN(:,1) contains the basis functions
%            themselves, dN(:,2) the derivatives of the basis functions 
%
% Function layout :
%
% 0. Input validation
%
% 1. Computation of the basis functions' values
%
%% Function Implementation

%% 0. Input validation
arguments
    xi (1, 1) double {mustBeInRange(xi, -1, 1)}
end

%% 1. Computation of the basis functions' values
dN = [(1 - xi)/2, -1/2
      (1 + xi)/2, 1/2];

end