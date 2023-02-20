function [dN, isInside] = ...
    computeBilinearBasisFunctionsAndFirstDerivatives(xi, eta)
%% COMPUTEBILINEARBASISFUNCTIONSANDFIRSTDERIVATIVES Returns the bilinear basis functions and their derivatives corresponding to a quadrilateral element
%
% Returns the bilinear basis functions and their derivatives corresponding 
% to a quadrilateral element.
%
%    Input :
%   xi,eta : The parameters on the parametric domain
%
%   Output :
%       dN : 4x3 array for which dN(:,1) contains the basis functions
%            themselves, dN(:,2) the derivatives of the basis functions 
%            with respect to xi-direction and dN(:,3) the derivatives of 
%            the basis functions with respect to eta-direction
% isInside : Flag on whether the point where the basis functions are to be
%            computed is inside or outside the quadrilateral
%
% Function layout :
%
% 0. Input validation
%
% 1. Read input
%
% 2. Compute the bilinear basis functions and their parametric derivatives
%
% 3. Check if the given point is inside the canonical quadrilateral
%
%% Function Implementation

%% 0. Input validation
arguments
    xi (1, 1) double {mustBeInRange(xi, -1, 1)}
    eta (1, 1) double {mustBeInRange(eta, -1, 1)}
end

%% 1. Read input

% Initialize output array
isInside = true;
    
%% 2. Compute the bilinear basis functions and their parametric derivatives

% Assign the values of the functions at the parametric locations
dN(:,1) = [(1 - xi)*(1 - eta)/4
           (1 - xi)*(1 + eta)/4
           (1 + xi)*(1 + eta)/4
           (1 + xi)*(1 - eta)/4];
           
% Derivatives dN/dxi
dN(:,2) = [-(1 - eta)/4
           -(1 + eta)/4
           (1 + eta)/4
           (1 - eta)/4];

% Derivatives dN/deta
dN(:,3) = [-(1 - xi)/4
            (1 - xi)/4
            (1 + xi)/4
           -(1 + xi)/4];
       
%% 3. Check if the given point is inside the canonical quadrilateral
if dN(1,1) < -1 || dN(1,1) > 1 || dN(2,1) < -1 || dN(2,1) > 1 || ...
        dN(3,1) < -1 || dN(3,1) > 1 || dN(4,1) < -1 || dN(4,1) > 1
    isInside = false;
end

end