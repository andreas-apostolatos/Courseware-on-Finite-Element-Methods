%[text] %[text:anchor:T_83E91B43] # Retrieval of the free degrees of freedom corresponding to the bilinear quadrilateral mesh
%[text:tableOfContents]{"heading":"**Table of Contents**"}
%[text] It returns the global numbering of the free Degrees of Freedom (DOFs) corresponding to the bilinear quadrilateral mesh created in section [Four-noded bilinear quadrilateral mesh](file:..\..\2_ReissnerMindlinPlate\main_Chapter2_AShearDeformablePlateElement.mlx:H_23BD88B7) in Live Script [`main_Chapter2_AShearDeformablePlateElement.mlx`](file:..\..\2_ReissnerMindlinPlate\main_Chapter2_AShearDeformablePlateElement.mlx) and in section [Four-noded bilinear quadrilateral mesh](file:C:\Interactive_Teaching_Material\FiniteElementMethodsCourseware_GitLab\3_TransverseShearLocking\3_TransverseShearLocking_ReissnerMindlinCantileverPlate\main_Chapter3_Locking_ReissnerMindlinCantileverPlate.mlx:H_4A4EF0A6) in Live Script [`main_Chapter3_Locking_ReissnerMindlinCantileverPlate.mlx`](file:C:\Interactive_Teaching_Material\FiniteElementMethodsCourseware_GitLab\3_TransverseShearLocking\3_TransverseShearLocking_ReissnerMindlinCantileverPlate\main_Chapter3_Locking_ReissnerMindlinCantileverPlate.mlx)[](file:..\..\2_ReissnerMindlinPlate\main_Chapter2_AShearDeformablePlateElement.mlx). This function is only relevant for the setup detailed in the latter Live Scripts.
function [freeDOFsQ1, homDOFsQ1] = ...
    getFreeDOFsBilinearQuadrilateralMesh(numNodesQ1, numEly)
%[text] %[text:anchor:H_D41A3DDC] ## **Function Description**
%[text] It returns the global numbering of the DOFs which do not lie on the Dirichlet boundary corresponding to the bilinear quadrilateral mesh for the cantilever plate problem.
%[text] %[text:anchor:H_C2A6403E]  **Input :**
%[text]  `numNodesQ1` : Number of nodes corresponding to the bilinear quadrilateral mesh
%[text]  `numElx`, `numEly` : Number of elements along the x- and the y-directions
%[text] %[text:anchor:H_CBAFD59D]  **Output :**
%[text]  `freeDOFsQ1` : The global numbering of free DOFs corresponding to the bilinear quadrilateral mesh
%[text]  `homDOFsQ1` : The global numbering of DOFs along the Dirichlet boundary corresponding to the bilinear quadrilateral mesh
%[text] %[text:anchor:H_D5ED9987] ## **Function Implementation**
%[text] %[text:anchor:H_D88B8E8E] ### Input validation
    arguments
        numNodesQ1 (1, 1) double {mustBeInteger, mustBePositive}
        numEly (1, 1) double {mustBeInteger, mustBePositive}
    end
%[text] %[text:anchor:H_E9B3FF26] ### Dirichlet DOFs
    freeDOFsQ1 = 1:3*numNodesQ1;
    homDOFsQ1 = 1:3*(numEly + 1);
    [~, idx] = ismember(homDOFsQ1, freeDOFsQ1);
    freeDOFsQ1(idx) = [];
%[text] 
end

%[appendix]{"version":"1.0"}
%---
