%[text] %[text:anchor:H_4C23294B] # Computation of the Penalty contribution matrix for the imposition of the weak Dirichlet boundary conditions
%[text:tableOfContents]{"heading":"**Table of Contents**"}
%[text] Considering the Penalty method for the application of weak Dirichlet boundary conditions for the Timoshenko beam problem, the following function returns the discrete Penalty matrix that results from the two Penalty terms, see Eqs. (4.1) and (4.2) [here](file:..\main_Chapter4_WeakBoundaryConditions_TimoshenkoBeam.mlx:M_BCB4701D), for the displacement and the rotation field, respectively.
function KPenalty = computePenaltyMtxWeakDirichletBoundaryConditions ...
    (msh, alphaU, alphaBeta, homDOFs, computeBasisFunctionsAndDerivs)
%[text] %[text:anchor:H_6E7F9CA0] ## Function description
%[text] This function returns the Penalty contribution to the system where the Dirichlet boundary conditions are imposed weakly. This contribution stems from the discretization of the following term,
%[text] $\\int\_0^L \\left\[ \\begin{array}{cc} \\delta w & \\delta \\beta \\end{array} \\right\] \\left\[  \\begin{array}{cc}  \\alpha\_w & 0 \\\\ 0 & \\alpha\_\\beta\\end{array} \\right\] \\left\[ \\begin{array}{c} w \\\\ \\beta \\end{array} \\right\] \\; \\text{d} \\Gamma$,
%[text] where $\\alpha\_w$ and $\\alpha\_\\beta$ stand for the Penalty parameters regarding the imposition of the Dirichlet boundary conditions for the vertical deflection and cross-sectional rotation, respectively.
%[text]  **Input :**
%[text]  `msh` : Struct-variable representing the Finite Element mesh containing the following information:
%[text]                                                                          .nodes : Two-dimensional array containing the Cartesian coordinates of each node in the mesh
%[text]                                                                     .elements : Two-dimensional array containing the global numbering of the nodes that belong to each element
%[text]  `alphaU` : The Penalty parameter $\\alpha\_w$ for the imposition of the Dirichlet boundary condition on the vertical deflections
%[text]  `alphaBeta` : The Penalty parameter $\\alpha\_\\beta$ for the imposition of the Dirichlet boundary condition on the cross-sectional rotations
%[text]  `homDOFs` : Global numbering of the DOFs that lie on the Dirichlet boundary
%[text] `computeBasisFunctionsAndDerivs` : Function handle for the computation of the basis functions and their derivatives
%[text] 
%[text]  **Output :**
%[text]  `Kpenalty` : The Penalty contribution to the stiffness matrix
%[text] %[text:anchor:H_AFEE50D5] ## Function implementation
%[text] %[text:anchor:H_F1A03792] ### Input validation
    arguments
        msh (1, 1) {mustHaveNodesAndElements}
        alphaU (1, 1) double
        alphaBeta (1, 1) double
        homDOFs (1, :) {mustBeInteger, mustBePositive}
        computeBasisFunctionsAndDerivs (1, 1) function_handle
    end
%[text] %[text:anchor:H_C22335FA] ### Read input
    numNodes = numel(msh.nodes(:, 1));
    numDOFs = 2*numNodes;
%[text] %[text:anchor:H_7BD6F95D] ### Initialization of the Penalty matrix
    KPenalty = zeros(numDOFs);
%[text] %[text:anchor:H_AE26E32D] ### Loop over all the DOFs on the Dirichlet boundary to compute the Penalty contribution
    for ii = 1:numel(homDOFs)
%[text] %[text:anchor:H_FAAD8BF9] #### Element containing the corresponding DOF
        idDOF = homDOFs(ii);
        idNode = ceil(idDOF/2);
        if mod(idDOF, 2)
            alpha = alphaU;
        else
            alpha = alphaBeta;
        end
        [idElem1, idElem2] = find(msh.elements==idNode);
        if isempty(idElem1) || isempty(idElem2)
            error("Node %d associated with DOF %d does not appear to belong to any element", idNode, idDOF)
        end
        if numel(idElem1) > 1 || numel(idElem2) >1
            if numel(idElem1) == numel(idElem2)
                error("Node %d associated with DOF %d is found in %d elements. Thus it cannot be a boundary DOF", idNode, idDOF, numel(idElem1))
            else
                error("This is unexpected")
            end
        end
        idElem = idElem1(1);
%[text] %[text:anchor:H_D6541941] #### Parametric coordinate of the node on the boundary associated with the constrained DOF
        if idElem2 == 1
            xi = -1;
        elseif idElem2 == 2
            xi = 1;
        else
            error("The boundary nodes must be either of the vertex nodes of the element")
        end
%[text] %[text:anchor:H_1484C9D2] #### Global numbering of the nodes at the associated element
        id = msh.elements(idElem, :);
        numNodesEl = numel(id);
%[text] %[text:anchor:H_B4F52171] #### Element freedom table for the pair of primal unknown displacement and rotation fields at the element on the Dirichlet boundary
        EFT = idDOF;
%[text] %[text:anchor:H_A2BFA7EE] #### Basis functions at the parametric coordinate associated with the constrained DOF
        dN = computeBasisFunctionsAndDerivs(xi);
        if numel(dN(:, 1)) ~= numNodesEl
            error("The chosen function handle for the computation of the basis functions %s returns %d basis " + ...
                "functions but the element has %d nodes", func2str(computeBasisFunctionsAndDerivs), numel(dN(:, 1)), ...
                numNodesEl);
        end
%[text] %[text:anchor:H_289CBFD3] #### Basis function matrix associated with the pair of primal unknown rotation and displacement fields
        NmtxU = dN(idElem2, 1);
%[text] %[text:anchor:H_AE936E77] #### Assembly of the element Penalty matrix contribution
        KPenalty(EFT, EFT) = KPenalty(EFT, EFT) + transpose(NmtxU)*alpha*NmtxU;
%[text] 
    end
%[text] %[text:anchor:H_2AA488B8] ### 
end

%[appendix]{"version":"1.0"}
%---
