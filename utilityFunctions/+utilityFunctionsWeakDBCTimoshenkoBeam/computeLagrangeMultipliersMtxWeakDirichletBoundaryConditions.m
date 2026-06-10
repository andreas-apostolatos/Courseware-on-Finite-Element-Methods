%[text] %[text:anchor:H_EF553DC7] # Lagrange Multipliers contribution matrix for the imposition of the weak Dirichlet boundary conditions
%[text:tableOfContents]{"heading":"**Table of Contents**"}
%[text] Considering the Penalty method for the application of weak Dirichlet boundary conditions for the Timoshenko beam problem, the following function returns the discrete Penalty matrix that results from the two Langrange Multipliers terms, see Eqs. (5.1) and (5.2) [here](file:..\main_Chapter4_WeakBoundaryConditions_TimoshenkoBeam.mlx:M_1434D535), for the displacement and the rotation field, respectively.
function Lambda = computeLagrangeMultipliersMtxWeakDirichletBoundaryConditions ...
    (msh, homDOFs, computeBasisFunctionsAndDerivs)
%[text] %[text:anchor:H_CEB4EB3D] ## Function description
%[text] This function returns the Lagrange Multiplier's matrix $\\Lambda$ to the system where the Dirichlet boundary conditions are imposed weakly, namely,
%[text] $\\left\[ \\begin{array}{cc} \\mathbf{K} & \\mathbf{\\Lambda} \\\\ \\mathbf{\\Lambda}^{\\text{T}} & \\mathbf{0} \\end{array} \\right\] \\left\[ \\begin{array}{c} \\mathbf{u} \\\\ \\mathbf{\\beta} \\end{array} \\right\] = \\left\[ \\begin{array}{c} \\mathbf{F} + \\overline{\\mathbf{F}} \\\\ \\mathbf{0} \\end{array} \\right\]$
%[text] This contribution stems from the discretization of the following term,
%[text] $\\int\_0^L \\left\[ \\begin{array}{cc} \\delta w & \\delta \\beta \\end{array} \\right\] \\left\[ \\begin{array}{c} \\lambda \\\\ \\mu \\end{array} \\right\] \\; \\text{d} \\Gamma$,
%[text] where $\\lambda$ and $\\mu$ stand for the Lagrange Multipliers regarding the imposition of the Dirichlet boundary conditions for the vertical deflection and cross-sectional rotation, respectively. These Lagrange Multipliers correspond to boundary tranverse shear forces and boundary cross-sectional rotations that are needed for the weak satisfaction of the Dirichlet boundary conditions.
%[text]  **Input :**
%[text]  `msh` : Struct-variable representing the Finite Element mesh containing the following information:
%[text]                                                                          .nodes : Two-dimensional array containing the Cartesian coordinates of each node in the mesh
%[text]                                                                     .elements : Two-dimensional array containing the global numbering of the nodes that belong to each element
%[text]  `homDOFs` : Global numbering of the DOFs that lie on the Dirichlet boundary
%[text] `computeBasisFunctionsAndDerivs` : Function handle for the computation of the basis functions and their derivatives
%[text] 
%[text]  **Output :**
%[text]  `Lambda` : The Lagrange Multipliers matrix $\\Lambda$
%[text] %[text:anchor:H_E0FBFAD1] ## Function implementation
%[text] %[text:anchor:H_F1A03792] ### Input validation
    arguments
        msh (1, 1) {mustHaveNodesAndElements}
        homDOFs (1, :) {mustBeInteger, mustBePositive}
        computeBasisFunctionsAndDerivs (1, 1) function_handle
    end
%[text] %[text:anchor:H_1BE1C17A] ### Read input
    numNodes = numel(msh.nodes(:, 1));
    numDOFs = 2*numNodes;
%[text] %[text:anchor:H_E35E1577] ### Initialization of the Lagrange Multipliers matrix
    numDOFsLM = numel(homDOFs);
    Lambda = zeros(numDOFs, numDOFsLM);
%[text] %[text:anchor:H_4036325D] ### Loop over all the DOFs on the Dirichlet boundary and assemble the Lagrange Multipliers matrix
    for ii = 1:numel(homDOFs)
%[text] %[text:anchor:H_A3A81DEE] #### Element that contains the corresponding DOF
        idDOF = homDOFs(ii);
        idNode = ceil(idDOF/2);
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
%[text] %[text:anchor:H_BCCF591E] #### Parametric coordinate of the node on the boundary associated with the constrained DOF
        if idElem2 == 1
            xi = -1;
        elseif idElem2 == 2
            xi = 1;
        else
            error("The boundary nodes must be either of the vertex nodes of the element")
        end
%[text] %[text:anchor:H_4E3CA74D] #### Global numbering of the nodes at the associated element
        id = msh.elements(idElem, :);
        numNodesEl = numel(id);
%[text] %[text:anchor:H_D9B2437F] #### Element Freedom Table (EFT) for the pair of primal unknown displacement and rotation fields at the element on the Dirichlet boundary
        EFT = idDOF;
%[text] %[text:anchor:H_EB9883E0] #### Element Freedom Table (EFT) for the pair of primal unknown displacement and rotation fields at the element on the Dirichlet boundary
        EFTLM = ii;
%[text] %[text:anchor:H_467AFC21] #### Basis functions at the parametric coordinate associated with the constrained DOF
        dN = computeBasisFunctionsAndDerivs(xi);
        if numel(dN(:, 1)) ~= numNodesEl
            error("The chosen function handle for the computation of the basis functions %s returns %d basis " + ...
                "functions but the element has %d nodes", func2str(computeBasisFunctionsAndDerivs), numel(dN(:, 1)), ...
                numNodesEl);
        end
%[text] %[text:anchor:H_73CFDC97] #### Basis function matrix associated with the pair of primal unknown displacement and rotation fields
        if idElem2 == 1
            NmtxU  = dN(1, 1);
        elseif idElem2 == 2
            NmtxU = dN(2, 1);
        end
%[text] %[text:anchor:H_5307087B] #### Basis function matrix associated with the Lagrange Multipliers
        NmtxBeta = 1;
%[text] %[text:anchor:H_68334639] #### Assembly of the element Lagrange Multipliers matrix contribution
        Lambda(EFT, EFTLM) = Lambda(EFT, EFTLM) + transpose(NmtxU)*NmtxBeta;
%[text] 
    end
%[text] 
end

%[appendix]{"version":"1.0"}
%---
