function [xiEta, hasConverged] = ...
    computePointCoordinatesOnCanonicalBilinearQuadrilateral...
    (X, XCorner, propNewtonRaphson)
%% COMPUTEPOINTCOORDINATESONCANONICALBILINEARQUADRILATERAL Returns the parametric coordinates of a Cartesian point in the parent space of a bilinear quadrilateral
%
% Returns the coordinates of a point with coordinates x within the 2D 
% quadrilateral x1-x2-x3-x4 in an arbitrary Cartesian space to the
% canonical bilinear quadrilateral. The applied method for solving the
% non-linear system is Newton-Rapson.
%
%             Input :
%                 X : The coordinates of the point into the Cartesian 
%                     system
%           XCorner : The coordinates of the corner nodes of the
%                     quadrilateral as follows :
%                           xCorner = [x1 x2 x3 x4]
% propNewtonRaphson : Parameters for the Newton-Raphson iterations :
%                        .eps : Termination criterion for the 2-norm of the
%                               residual
%                      .maxIt : Maximum number of Newton-Rapshon iterations
%
%           Output :
%            xiEta : The coordinates of the point in the canonical space of 
%                    the bilinear quadrilateral xiEta = [xi eta]
%     hasConverged : Flag on the convergence of the Newton iterations
%
% Function layout :
%
% 0. Input validation
% 
% 1. Read input
%
% 2. Loop over all the Newton iterations
% ->
%
%    2i. Compute the bilinear basis functions and its derivatives at u
%
%   2ii. Compute the right-hand side residual
%
%  2iii. Compute the Jacobian
%
%   2iv. Solve the linear system
%
%    2v. Check condition for convergence
%
%   2vi. Update the iteration counter
% <-
%
% 3. Check if convergence has been achieved
%
%% Function Implementation

%% 0. Input validation
arguments
    X (2, 1) double
    XCorner (2, 4) double
    propNewtonRaphson (1, 1)
end
if ~isfield(propNewtonRaphson, 'eps')
    error("Struct input variable 'propNewtonRaphson' should contain a numeric " + ...
        "field 'eps'")
else
    if ~isnumeric(propNewtonRaphson.eps)
        error("Field 'propNewtonRaphson.eps' should be numeric")
    else
        if ~isscalar(propNewtonRaphson.eps)
            error("Field 'propNewtonRaphson.eps' should scalar")
        end
    end
end
if ~isfield(propNewtonRaphson, 'maxIt')
    error("Struct input variable 'propNewtonRaphson' should contain a numeric " + ...
        "field 'maxIt'")
else
    if ~isnumeric(propNewtonRaphson.maxIt)
        error("Field 'propNewtonRaphson.maxIt' should be numeric")
    else
        if ~isscalar(propNewtonRaphson.maxIt)
            error("Field 'propNewtonRaphson.maxIt' should scalar")
        else
            if propNewtonRaphson.maxIt ~= floor(propNewtonRaphson.maxIt) || ...
                propNewtonRaphson.maxIt <= 0
                error("Field 'propNewtonRaphson.maxIt' should be a strictly positive integer")
            end
        end
    end
end

%% 1. Read input

% Coordinates of point x into the canonical space. Initialize them into the
% center of the canonical bilinear quadrilateral
xiEta = zeros(2,1);

% Initialize counter
counter = -1;

% Initialize convergence flag to true
hasConverged = true;

%% 2. Loop over all the Newton iterations
while counter <= propNewtonRaphson.maxIt
    %% 2i. Compute the bilinear basis functions and its derivatives at xi
    dN = computeBilinearBasisFunctionsAndFirstDerivatives ...
        (xiEta(1,1), xiEta(2,1));
    
    %% 2ii. Compute the right-hand side residual
    
    % Compute the physical coordinates of the point with canonical
    % coordinates xi
    xPredicted = zeros(2,1);
    DxPredictedDxi = zeros(2,1);
    DxPredictedDeta = zeros(2,1); 
    
    % Loop over all the basis functions and add contributions
    for i = 1:length(dN(:, 1))
        % The coordinates of the point
        xPredicted = xPredicted + dN(i, 1)*XCorner(:, i);
        
        % The derivatives of the position vector w.r.t. tp xi
        DxPredictedDxi = DxPredictedDxi + dN(i, 2)*XCorner(:, i);
        
        % The derivatives of the position vector w.r.t. tp eta
        DxPredictedDeta = DxPredictedDeta + dN(i, 3)*XCorner(:, i);
    end
    
    % Compute the residual vector
    residual = X - xPredicted;
    
    %% 2iii. Compute the Jacobian
    J = [DxPredictedDxi(1, 1) DxPredictedDeta(1, 1)
         DxPredictedDxi(2, 1) DxPredictedDeta(2, 1)];
    
    %% 2iv. Solve the linear system
    delta = J\residual;
    
    % Update the parametric locations
    xiEta = xiEta + delta;

    %% Enforce parametric coordinates to stay within the parametric space [xi,eta] \in [-1,+1]x[-1,+1]
    xiEta(1, 1) = enforceParametricCoordinate(xiEta(1, 1), [-1, +1]);
    xiEta(2, 1) = enforceParametricCoordinate(xiEta(2, 1), [-1, +1]);
    
    %% 2v. Check condition for convergence
    if norm(residual) <= propNewtonRaphson.eps
        break;
    end
    
    %% 2vi. Update the iteration counter
    counter = counter + 1;
end

%% 3. Check if convergence has been achieved
if counter > propNewtonRaphson.maxIt
    hasConverged = false;
end

end

%% ENFORCEPARAMETRICCOORDINATE Enforces parametric coordinate to stay in the parametric space
function xi = enforceParametricCoordinate(xi, span)
    arguments
        xi (1, 1) double
        span (1, 2) double
    end
    if xi < span(1)
        xi = span(1);
    end
    if xi > span(2)
        xi = span(2);
    end
end