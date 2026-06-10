%[text] %[text:anchor:T_4A77ED39] # Plot an One-Dimensional Mesh
%[text] Plot an one-dimensional mesh including the numbering of the nodes and the numbering of the elements.
%[text:tableOfContents]{"heading":"**Table of Contents**"}
function plotMesh1d(msh)
%[text] %[text:anchor:H_214B5634] ## Function Description
%[text] Plot an one-dimensional mesh together together with the node and element numberings
%[text] %[text:anchor:H_AEB0D076]  **Input :**
%[text]  `msh` : Nodal coordinates of the end vertices of the finite element
%[text]                                .`nodes` : Array containing the coordinates of the nodes in the mesh
%[text]                          .`elements` : Array containing the IDs of the nodes of each element in the mesh
%[text] %[text:anchor:H_BE68C1F8] ## Function Implementation
%[text] %[text:anchor:H_F940506C] ### Input validation
    arguments
        msh (1, 1) {mustHaveNodesAndElements}
    end
%[text] %[text:anchor:H_3529A059] ### Plot the one-dimensional mesh
    plot(msh.nodes, zeros(numel(msh.nodes), 1), '-o');
%[text] %[text:anchor:H_1DF4C8BD] ### Display the node numbering on the mesh
    for ii = 1:length(msh.nodes(:, 1))
        text(msh.nodes(ii), max(msh.nodes)/50, sprintf("X_{%d}", ii));
    end
%[text] %[text:anchor:H_538290C7] ### Display the element numbering on the mesh
    for ii = 1:length(msh.elements(:, 1))
        text(mean(msh.nodes(msh.elements(ii, :))), -max(msh.nodes)/70, strcat("e_{", strcat(num2str(ii), "}")));
    end
%[text] 
end

%[appendix]{"version":"1.0"}
%---
