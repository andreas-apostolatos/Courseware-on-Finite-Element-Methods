%[text] %[text:anchor:T_8AF2E205] # Plot a Two-Dimensional Finite Element Mesh
%[text] Plot a two-dimensional mesh including the numbering of the nodes and the numbering of the elements.
%[text:tableOfContents]{"heading":"**Table of Contents**"}
function plotMesh2d(msh, numVertices, isPlotNumberings)
%[text] %[text:anchor:H_586A5901] ## **Function Description**
%[text] Plot a two-dimensional mesh together together with the node and element numberings. The method works for all types of quadrilaterals (e.g. bilinear, biquadratic, etc.)
%[text]  **Input :**
%[text]  `msh` : Nodal coordinates of the end vertices of the finite element
%[text]                                          .`nodes` : Array containing the coordinates of the nodes in the mesh.
%[text]                                     .`elements` : Array containing the IDs of the nodes of each element in the mesh
%[text]  `numVertices` : The number of vertex nodes corresponding to the quadrilateral elements. For higher order meshes, the number of vertex nodes is smaller than the total number of nodes, as the total number of nodes contains also the interior edge and domain nodes
%[text]  `isPlotNumberings` : Flag to indicate whether the node and element numberings should be plotted
%[text] 
%[text] %[text:anchor:H_C15C19CD] ## Function Implementation
%[text] %[text:anchor:H_948A5DE3] ### Input validation
    arguments
        msh (1, 1) {mustHaveNodesAndElements}
        numVertices (1, 1) {mustBeInteger, mustBePositive}
        isPlotNumberings (1, 1) logical
    end
%[text] %[text:anchor:H_CAFB5CED] ### Read input
    numNodesEl = width(msh.elements);
    if numNodesEl == 3
        numVerticesEl = 3;
    else
        numVerticesEl = 4;
    end
%[text] %[text:anchor:H_F7FC26A7] ### Plot the elements in the mesh
    if numVerticesEl == 3
        trimesh(msh.elements(:, 1:numVerticesEl), msh.nodes(1:numVertices, 1), msh.nodes(1:numVertices, 2), msh.nodes(1:numVertices, 3), ...
                FaceColor=[217 218 219]/255, EdgeColor='k', FaceAlpha=0.5);
    elseif numVerticesEl == 4
        patch('Faces', msh.elements(:,1:4), ...
            'Vertices', msh.nodes(1:numVertices,1:3), ...
            'FaceColor', [217 218 219]/255, ...
            'EdgeColor', 'k', ...
            'FaceAlpha', 0.5);
    end
    view(3)
    hold on;
%[text] %[text:anchor:H_E4DC393E] ### Plot the nodes in the mesh
    plot3(msh.nodes(:, 1), msh.nodes(:, 2), msh.nodes(:, 3), 'ok');
%[text] %[text:anchor:H_CDA574AF] ### Display the node numbering on the mesh
    if isPlotNumberings
        for ii = 1:length(msh.nodes(:, 1))
            text(msh.nodes(ii, 1), msh.nodes(ii, 2), msh.nodes(ii, 3) + max(msh.nodes, [], "all", "includenan")/20, ...
                sprintf("X_{%d}", ii));
        end
    end
%[text] %[text:anchor:H_A05D9E7A] ### Display the element numbering on the mesh
    if isPlotNumberings
        for ii = 1:length(msh.elements(:, 1))
            text(mean(msh.nodes(msh.elements(ii, :), 1)), mean(msh.nodes(msh.elements(ii, :), 2)), ...
                max(msh.nodes, [], "all", "includenan")/30, strcat(strcat("e_{", num2str(ii)), "}"));
        end
    end
%[text] 
    hold off;
end

%[appendix]{"version":"1.0"}
%---
