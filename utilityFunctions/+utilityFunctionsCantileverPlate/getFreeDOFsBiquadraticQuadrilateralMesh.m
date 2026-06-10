%[text] %[text:anchor:H_12B57CA4] # Retrieval of the free degrees of freedom corresponding to the biquadratic quadrilateral mesh
%[text:tableOfContents]{"heading":"**Table of Contents**"}
%[text] It returns the global numbering of the free Degrees of Freedom (DOFs) corresponding to the bilinear quadrilateral mesh created in section [Nine-noded biquadratic quadrilateral mesh](file:..\main_Chapter3_Locking_ReissnerMindlinCantileverPlate.mlx:H_4E93D9EB) in Live Script [main\_Chapter3\_Locking\_ReissnerMindlinCantileverPlate.mlx](file:C:\Interactive_Teaching_Material\FiniteElementMethodsCourseware_GitLab\3_TransverseShearLocking\3_TransverseShearLocking_ReissnerMindlinCantileverPlate\main_Chapter3_Locking_ReissnerMindlinCantileverPlate.mlx). This function is only relevant for the setup detailed in the latter Live Script.
function [freeDOFsQ2, homDOFsQ2] = ...
    getFreeDOFsBiquadraticQuadrilateralMesh(numNodesQ1, numNodesQ2, numEly)
%[text] %[text:anchor:H_D41A3DDC] ## Function Description
%[text] It returns the global numbering of the DOFs which do not lie on the Dirichlet boundary corresponding to the biquadratic quadrilateral mesh for the problem at hand, namely, the cantilever plate.
%[text] %[text:anchor:H_C2A6403E]  **Input :**
%[text]  `numNodesQ1` : Number of nodes corresponding to the bilinear quadrilateral mesh onto which the biquadratic quadrilateral mesh is based
%[text]  `numEly` : Number of elements along the y-direction
%[text] %[text:anchor:H_CBAFD59D]  **Output :**
%[text]  `freeDOFsQ2` : The global numbering of the free DOFs corresponding to the biquadratic quadrilateral mesh
%[text]  `homDOFsQ2` : The global numbering of the constrained DOFs corresponding to the biquadratic quadrilateral mesh
%[text] %[text:anchor:H_D5ED9987] ## Function Implementation
%[text] %[text:anchor:H_D88B8E8E] ### Input validation
    arguments
        numNodesQ1 (1, 1) double {mustBeInteger, mustBePositive}
        numNodesQ2 (1, 1) double {mustBeInteger, mustBePositive}
        numEly (1, 1) double {mustBeInteger, mustBePositive}
    end
%[text] %[text:anchor:H_E9B3FF26] ### Dirichlet DOFs
freeDOFsQ2 = 1:3*numNodesQ2;
homDOFsQ2 = horzcat(1:3*(numEly + 1), 3*numNodesQ1+1:3*numNodesQ1+3*numEly);
[~, idx] = ismember(homDOFsQ2, freeDOFsQ2);
freeDOFsQ2(idx) = [];
%[text] 
end

%[appendix]{"version":"1.0"}
%---
