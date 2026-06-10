%[text] %[text:anchor:H_E94CADDE] # Retrieval of the free degrees of freedom corresponding to the bilinear quadrilateral mesh
%[text:tableOfContents]{"heading":"**Table of Contents**"}
%[text] It returns the global numbering of the free Degrees of Freedom (DOFs) corresponding to the bilinear quadrilateral mesh created in section [Four-noded bilinear quadrilateral mesh](file:..\main_Chapter3_Locking_ReissnerMindlin_SquarePlate.mlx:H_711714C3) in Live Script [main\_Chapter3\_Locking\_ReissnerMindlin\_SquarePlate.mlx](file:..\main_Chapter3_Locking_ReissnerMindlin_SquarePlate.mlx). This function is only relevant for the setup detailed in the latter Live Script.
function [freeDOFsQ1, homDOFsQ1, idxX, homDOFsQ1X] = getFreeDOFsBilinearQuadrilateralMesh ...
    (numNodesQ1, numElx, numEly)
%[text] %[text:anchor:H_D41A3DDC] ## **Function Description**
%[text] It returns the global numbering of the DOFs which do not lie on the Dirichlet boundary corresponding to the bilinear quadrilateral mesh for the problem at hand, namely, the two-sided clamped rectangular plate.
%[text] %[text:anchor:H_C2A6403E]  **Input :**
%[text]  `numNodesQ1` : Number of nodes corresponding to the bilinear quadrilateral mesh
%[text]  `numElx`, `numEly` : Number of elements along the x- and the y-directions
%[text] %[text:anchor:H_CBAFD59D]  **Output :**
%[text]  `freeDOFsQ1` : The global numbering of free DOFs corresponding to the bilinear quadrilateral mesh
%[text]  `homDOFsQ1` : The global numbering of DOFs along the Dirichlet boundary corresponding to the bilinear quadrilateral mesh
%[text]  `idxX` : Array containing the IDs of the nodes that lie on the Dirichlet boundary along the x-direction
%[text]  `homDOFsQ1X` : Array containing the global numbering of the DOFs which are constrained along the x-direction
%[text] %[text:anchor:H_D5ED9987] ## **Function Implementation**
%[text] %[text:anchor:H_D88B8E8E] ### Input validation
    arguments
        numNodesQ1 (1, 1) double {mustBeInteger, mustBePositive}
        numElx (1, 1) double {mustBeInteger, mustBePositive}
        numEly (1, 1) double {mustBeInteger, mustBePositive}
    end
%[text] %[text:anchor:H_E9B3FF26] ### Dirichlet DOFs
    freeDOFsQ1 = 1:3*numNodesQ1;
    homDOFsQ1Y = 1:3*(numEly + 1);
    idxX = numEly+1:numEly+1:(numElx + 1)*(numEly + 1);
    homDOFsQ1X = sort([3*idxX-2 3*idxX-1 3*idxX]);
    homDOFsQ1 = horzcat(homDOFsQ1Y, homDOFsQ1X);
    homDOFsQ1 = unique(sort(homDOFsQ1));
    [~, idx] = ismember(homDOFsQ1, freeDOFsQ1);
    freeDOFsQ1(idx) = [];
%[text] 
end

%[appendix]{"version":"1.0"}
%---
