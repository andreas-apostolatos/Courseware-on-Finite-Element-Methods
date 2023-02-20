function mustHaveNodesAndElements(msh)
%% MUSTCONTAINNODESANDELEMENTS Verifies that a mesh is defined appropriately
%
% Checks whether input MATLAB®-struct 'msh' contains the expected fields 
% 'nodes' and 'elements' and whether these have the expected types
% 
% Input :
% 
%   msh : Nodal coordinates of the end vertices of a finite element
%       .nodes : Array containing the coordinates of the nodes in the mesh
%    .elements : Array containing the IDs of the nodes of each element in 
%                the mesh
%
% Function layout :
%
% 1. Verification of the existence of the 'nodes'-field of the MATLAB-struct
%
% 2. Verification of the existence of the 'elements'-field of the MATLAB-struct
%
%% Function Implementation

%% 1. Verification of the existence of the 'nodes'-field of the MATLAB-struct
if ~isfield(msh, 'nodes')
    error("Input struct 'msh' should contain a field 'nodes'")
else
    if ~isnumeric(msh.nodes)
        error("Field 'msh.nodes' should be a numeric array")
    else
        if isscalar(msh.nodes)
            error("The array of nodes should contain at least two")
        end
    end
end

%% 2. Verification of the existence of the 'elements'-field of the MATLAB-struct
if ~isfield(msh, 'elements')
    error("Input struct 'msh' should contain a field 'elements'")
else
    if ~isnumeric(msh.nodes)
        error("Field 'msh.elements' should be a numeric array")
    else
        if isscalar(msh.nodes)
            error("The array of elements should contain at least one element with more than one nodes")
        end
    end
end

end