%[text] %[text:anchor:H_12B57CA4] # Retrieval of the free degrees of freedom corresponding to the biquadratic quadrilateral mesh
%[text:tableOfContents]{"heading":"**Table of Contents**"}
%[text] It returns the global numbering of the free Degrees of Freedom (DOFs) corresponding to the bilinear quadrilateral mesh created in section [Nine-noded biquadratic quadrilateral mesh](file:..\main_Chapter3_Locking_ReissnerMindlin_SquarePlate.mlx:H_4E93D9EB) in Live Script [main\_Chapter3\_Locking\_ReissnerMindlin\_SquarePlate.mlx](file:..\main_Chapter3_Locking_ReissnerMindlin_SquarePlate.mlx). This function is only relevant for the setup detailed in the latter Live Script.
function [freeDOFsQ2, homDOFsQ2] = getFreeDOFsBiquadraticQuadrilateralMesh ...
    (numNodesQ1, numNodesQ2, numEly, homDOFsQ1X, idxX, biasX)
%[text] %[text:anchor:H_EBFEFD5E] ## Function description
%[text] It returns the global numbering of the DOFs which do not lie on the Dirichlet boundary corresponding to the biquadratic quadrilateral mesh for the problem at hand, namely, the two-sided clamped rectangular plate.
%[text] %[text:anchor:H_C2A6403E]  **Input :**
%[text]  `numNodesQ1` : Number of nodes corresponding to the bilinear quadrilateral mesh onto which the biquadratic quadrilateral mesh is based
%[text]  `numNodesQ2` : Number of nodes corresponding to the biquadratic quadrilateral mesh
%[text]  `numEly` : Number of elements along the $Y${"editStyle":"visual"}-direction
%[text]  `homDOFsQ1X` : Global numbering of the DOFs corresponding to the Dirichlet boundary along the $X${"editStyle":"visual"}-direction
%[text]  `idxX` : Array containing the IDs of the nodes that lie on the Dirichlet boundary along the $X${"editStyle":"visual"}-direction
%[text]  `biasX` : Number of nodes of the background bilinear quadrilateral mesh and the set of intermediate for each edge along the $Y${"editStyle":"visual"}-direction
%[text] 
%[text] %[text:anchor:H_CBAFD59D]  **Output :**
%[text]  `freeDOFsQ2` : The global numbering of the free DOFs corresponding to the biquadratic quadrilateral mesh
%[text]  `homDOFsQ2` : The global numbering of DOFs along the Dirichlet boundary corresponding to the biquadratic quadrilateral mesh
%[text] %[text:anchor:H_2728769D] ## Function implementation
%[text] %[text:anchor:H_D88B8E8E] ### Input validation
    arguments
        numNodesQ1 (1, 1) double {mustBeInteger, mustBePositive}
        numNodesQ2 (1, 1) double {mustBeInteger, mustBePositive}
        numEly (1, 1) double {mustBeInteger, mustBePositive}
        homDOFsQ1X (1, :) double
        idxX (1, :) double
        biasX (1, 1) double
    end
%[text] %[text:anchor:H_E9B3FF26] ### Dirichlet DOFs
freeDOFsQ2 = 1:3*numNodesQ2;
homDOFsQ2Y = horzcat(1:3*(numEly + 1), 3*numNodesQ1+1:3*numNodesQ1+3*numEly);
homDOFsQ2 = horzcat(homDOFsQ2Y, homDOFsQ1X);
homDOFsQ2 = unique(sort(homDOFsQ2));
idxXX = idxX(1:end - 1) + biasX;
homDOFsQ2XX = sort([3*idxXX-2 3*idxXX-1 3*idxXX]);
homDOFsQ2 = unique(sort(horzcat(homDOFsQ2, homDOFsQ2XX)));
[~, idx] = ismember(homDOFsQ2, freeDOFsQ2);
freeDOFsQ2(idx) = [];
%[text] 
end

%[appendix]{"version":"1.0"}
%---
