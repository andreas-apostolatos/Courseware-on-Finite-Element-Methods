function [xs, ys, zs] = createSupports3D(nds, homDOFs, scalVer)
%% CREATESUPPORT3D Returns the coordinates for the visualization of the supports on a plate as planar triangles in 3d
%
% It returns the coordinates in X-, Y- and Z-directions for the creation of
% triagles in three-dimensions indicating the directions of supports in a
% plate configuration.
%
% - A triangle aligned with the X-Y plane indicates a fixed support with
% respect to the vertical deflection
%
% - A triangle aligned with the Y-Z plane indicates a clamped support with
% respect to the rotation around the X-axis
%
% - A triangle aligned with the X-Z plane indicates a clamped support with
% respect to the rotation around the Y-axis
%
%       Input : 
%         nds : The set of nodes in the Finite Element mesh
%     homDOFs : Global numbering of the suported DOFs
%     scalVer : Vertical scaling of the Z-dimension of the supports in case
%               a field is displayed along with the mesh in the Z-direction
%
%      Output :
%  xs, ys, zs : The coordinates of the triangles graphically representing
%               the boundary conditions
%
% Function layout :
%
% 0. Input validation
%
% 1. Initialize output
%
% 2. Axis alignment necessary for placing the triangles at the corresponding locations
%
% 3. Loop over all the DOFs along the Dirichlet boundary
% ->
%    3i. ID of node associated with the DOF
%   3ii. Create the triangle aligned with the corresponding plane (X-Y, Y-Z, X-Z) depending on the DOF direction
% <-
%
%% Function Implementation

%% 0. Input validation
arguments
    nds (:, 3) double
    homDOFs (1, :) double {mustBeInteger, mustBePositive}
    scalVer (1, 1) double {mustBeReal}
end

%% 1. Initialize output
numDOFs = numel(homDOFs);
xs = zeros(numDOFs, 4);
ys = zeros(numDOFs, 4);
zs = zeros(numDOFs, 4);

%% 2. Axis alignment necessary for placing the triangles at the corresponding locations
up = max(max(nds));
lo = min(min(nds));
axisAlignment = (up - lo)/50;

%% 3. Loop over all the DOFs along the Dirichlet boundary
for ii = 1:numDOFs
    %% 3i. ID of node associated with the DOF
    homDOF = homDOFs(ii);
    id = ceil(homDOF/3);

    %% 3ii. Create the triangle aligned with the corresponding plane (X-Y, Y-Z, X-Z) depending on the DOF direction
    if mod(homDOF, 3) == 0 % support on rotation around Y-axis
        xs(ii, 1) = nds(id, 1);
        xs(ii, 2) = nds(id, 1) -  1.732*axisAlignment;
        xs(ii, 3) = nds(id, 1) - 1.732*axisAlignment;
        xs(ii, 4) = xs(ii, 1);
        ys(ii, 1) = nds(id, 2);
        ys(ii, 2) = nds(id, 2) + axisAlignment;
        ys(ii, 3) = nds(id, 2) - axisAlignment;
        ys(ii, 4) = ys(ii, 1);
        zs(ii, 1:4) = nds(id, 3);
    elseif mod(homDOF, 3) == 2 % support on vertical deflection
        zs(ii, 1) = nds(id, 3);
        zs(ii, 2) = nds(id, 3) -  scalVer*1.732*axisAlignment;
        zs(ii, 3) = nds(id, 3) - scalVer*1.732*axisAlignment;
        zs(ii, 4) = zs(ii, 1);
        ys(ii, 1) = nds(id, 2);
        ys(ii, 2) = nds(id, 2) + axisAlignment;
        ys(ii, 3) = nds(id, 2) - axisAlignment;
        ys(ii, 4) = ys(ii, 1);
        xs(ii, 1:4) = nds(id, 1);
    elseif mod(homDOF, 3) == 1 % support on rotation around X-axis
        xs(ii, 1) = nds(id, 1);
        xs(ii, 2) = nds(id, 1) - axisAlignment;
        xs(ii, 3) = nds(id, 1) + axisAlignment;
        xs(ii, 4) = xs(ii, 1);
        ys(ii, 1) = nds(id, 2);
        ys(ii, 2) = nds(id, 2) - 1.732*axisAlignment;
        ys(ii, 3) = nds(id, 2) - 1.732*axisAlignment;
        ys(ii, 4) = ys(ii, 1);
        zs(ii, 1:4) = nds(id, 3);
    else
        error("Array 'homDOFs' should only contain integers corresponding to three" + ...
            "types of support: a vertical deflection and two cross-sectional" + ...
            "rotations. However, an extra unsupported type of support is found  " + ...
            "corresponding to DOF %i with value %d", ii, homDOF);
    end
end

end