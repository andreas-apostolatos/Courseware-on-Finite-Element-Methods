function [idEl, xiTilde] = findElementByPhysicalCoordinate1d(mshL1, X)
%% FINDELEMENTBYPHYSICALCOORDINATE1D Returns the element id that contains a Cartesian point for an one-dimensional two-noded linear mesh
%
% Returns element ID that contains the given physical coordinate and its 
% xi-parametric coordinate with respect to the parametric space of this
% element
%
%       Input : 
%       mshL1 : One-dimensional Finite Element mesh
%                   .nodes : Nodal coordinates in the one-dimensional mesh,
%                            e.g. mshL1.nodes(2, :) is a row-vector
%                            containing the X-, Y- and Z-coordinates of the
%                            second node in the mesh
%                .elements : Nodal ID's per element in the one-dimensional
%                            mesh, e.g. mshL1(3, :) is a row-vector
%                            containing all the nodal ID's of the third
%                            element in the mesh
%           X : The coordinates of the query point
%
%      Output :
%        idEl : The ID of the element that contains node with coordinates
%               xCoord in the one-dimensional Finite Element mesh mshL1
%     xiTilde : The parametic location of the query point in the element
%               of the Finite Element mesh that contains it
%
% Function layout :
%
% 0. Input validation
%
% 1. Read input
%
% 2. Setup the symbolic expressions of the linear basis functions
%
% 3. Find the element containing the query point and its local (parametric coordinates)
% ->
%    3i. Construct the position vector corresponding to the current element
%
%   3ii. Construct the inverse position vector corresponding to the current element
%
%  3iii. Compute the parametric coordinate of the evaluation point in the current element using the element's inverse parametric representation
%
%   3iv. Check whether the computed parametric coordinate is within the closed space [-1, +1] in which case the element containing the node has been found
% <-
%
% 4. Check the validity of the element found containing the current point
%
%% Function main body

%% 0. Input validation
arguments
    mshL1 {mustHaveNodesAndElements(mshL1)}
    X (1,1) double
end

%% 1. Read input
idEl = nan;

%% 2. Setup the symbolic expressions of the linear basis functions
syms pos_vct_L1(xi) N1_L1(xi) N2_L1(xi)
N1_L1(xi) = (1 - xi)/2;
N2_L1(xi) = (1 + xi)/2;

%% 3. Find the element containing the query point and its local (parametric coordinates)
numEl_L1 = numel(mshL1.elements(:, 1));
for jj = 1:numEl_L1
    %% 3i. Construct the position vector corresponding to the current element
    pos_vct_L1(xi) = [N1_L1 N2_L1]*mshL1.nodes(mshL1.elements(jj, :));

    %% 3ii. Construct the inverse position vector corresponding to the current element
    pos_vct_inv = finverse(pos_vct_L1);
    
    %% 3iii. Compute the parametric coordinate of the evaluation point in the current element using the element's inverse parametric representation
    xiTilde = pos_vct_inv(X);
    
    %% 3iv. Check whether the computed parametric coordinate is within the closed space [-1, +1] in which case the element containing the node has been found
    if xiTilde >= -1 && xiTilde <= 1
        idEl = jj;
        break;
    end

end

%% 4. Check the validity of the element found containing the current point
if isnan(idEl)
    error("The point with coordinate %d was found in no element in the mesh", X);
end

end