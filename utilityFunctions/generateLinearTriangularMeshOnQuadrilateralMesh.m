%[text] %[text:anchor:T_9805BF00] # **Generation of a linear triangular mesh using the Delaunay triangulation**
%[text] It uses an existing set of nodes to create a triangular mesh by means of the Delaunay triangulation
%[text:tableOfContents]{"heading":"**Table of Contents**"}
function mshT1 = generateLinearTriangularMeshOnQuadrilateralMesh(nodesX, nodesY)
%[text] %[text:anchor:H_930C6DF3] ## Function description
%[text] Returns a linear triangular mesh for the following rectangular geometry using nodes that are gene
%[text]  **Input :**
%[text]  `nodesX,` `nodesY` : Two-dimensional array `[numNodesX, numNodesY]` containing the $X$- and $Y$-coordinates of the nodes in a provided quadrilateral mesh, respectively
%[text]  **Output :**
%[text]  `mshT1` : MATLAB-struct containing the nodes and elements in the generated linear triangular mesh. If the Community Toolbox cannot be detected, then this is returned as an empty array.
%[text]                                         .`nodes` : Array `[numNodes, 3]` containing the nodes in the linear triangular mesh
%[text]  `.elements` : Array `[numEl, 3]` containing the node connectivity by means of the node numbering for each element in the linear triangular mesh
%[text] %[text:anchor:H_2F5E0E0A] ## **Function implementation**
%[text] %[text:anchor:H_3877BB59] ### **Read input**
    arguments
        nodesX (:, 1) {mustBeNumeric, mustBeReal}
        nodesY (:, 1) {mustBeNumeric, mustBeReal, mustBeEqualSize(nodesY, nodesX)}
    end
    if numel(nodesX) ~= numel(nodesY)
        error("Arrays 'nodesX' and 'nodesY' must have the same number of elements")
    end
%[text] %[text:anchor:H_6F973E73] ### **Mesh generation**
    mshT1.nodes = [nodesX nodesY zeros(numel(nodesX), 1)];
    dt = delaunay(nodesX, nodesY);
    mshT1.elements = dt;
%[text] 
end

%[appendix]{"version":"1.0"}
%---
