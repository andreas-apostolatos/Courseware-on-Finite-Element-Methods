function dN = computeConstantBasisFunctionAndFirstDerivatives(xi)
%% COMPUTECONSTANTBASISFUNCTIONSANDFIRSTDERIVATIVES Returns the constant basis functions and their derivatives for an one-dimensional element
%
% Returns the constant shape function and its derivative corresponding to 
% an one-dimensional element
%
%    Input :
%       xi : The coordinate on the parametric domain [-1, 1]
%
%   Output :
%       dN : 1x2 array for which dN(1, 1) contain the basis function
%            itself and dN(1, 2) its derivative
%
% Function layout :
%
% 0. Input validation
%
% 1. Compute the constant basis function and its parametric derivatives
%
%% Function main body

%% 0. Input validation
arguments
    xi (1, 1) double {mustBeInRange(xi, -1, 1)} = 0
end

%% 1. Compute the constant basis function and its parametric derivatives
dN = [1, 0];

end