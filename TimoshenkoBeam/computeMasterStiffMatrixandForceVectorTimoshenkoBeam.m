%[text] %[text:anchor:T_30109DD6] # Master Stiffness Matrix and Body Load Vector for the Timoshenko beam
%[text] Computation of the master stiffness matrix and master load vector corresponding to the standard Timoshenko beam formulation with two degrees of freedom per node $\\left(w,\\beta \\right)${"editStyle":"visual"}, the vertical deflection and the cross-section rotation. The master stiffness matrix and master load vector for the standard Timoshenko beam are given by the following relations:
%[text] $\\begin{array}{ll}\n\\mathbf{K} = \\int\_{L} \\mathbf{B}\_{\\text{s}}^{\\text{T}} \\mathbf{C}\_{\\text{s}} \\mathbf{B}\_{\\text{s}} \\, \\text{d} L+ \\int\_{L} \\mathbf{B}\_{\\text{b}}^{\\text{T}} \\mathbf{C}\_{\\text{b}} \\mathbf{B}\_{\\text{b}} \\: \\text{d} L , & \\text{(1.1)} \\\\\n\\mathbf{F} = \\int\_{L} \\mathbf{N}^{\\text{T}} \\left\[ \\begin{array}{c} \\bar{p} \\\\ \\bar{m} \\end{array} \\right\] \\: \\text{d} L + \\left\[\\mathbf{N}^{\\text{T}} \\left\[ \\begin{array}{c} \\bar{P} \\\\ \\bar{M} \\end{array} \\right\]  \\right\]\_0^{L}. & \\text{(1.2)}\n\\end{array}$
%[text:tableOfContents]{"heading":"**Table of Contents**"}
function [K, F, invHG] = computeMasterStiffMatrixandForceVectorTimoshenkoBeam ...
    (msh, computeElStiffMtxForceVct, computeBasisFunctionsAndDerivsU, ...
    computeBasisFunctionsAndDerivsBeta, propStr)
%[text] %[text:anchor:H_E15E957B] ## **Function Description**
%[text] Returns the master stiffness matrix corresponding to a Timoshenko beam.
%[text]  **Input** :
%[text]  `msh` : Structure containing the nodes and the elements in a quadrilateral finite element mesh
%[text]                                                                              .`nodes` : Array with the nodal coordinates in the finite element mesh
%[text]                                                                        .`elements` : Array with the nodal numbering per element in the finite element mesh
%[text]  `computeElStiffMtxForceVct` : Function handle to the computation of the element stiffness matrix of the Timoshenko beam element
%[text]  `computeBasisFunctionsAndDerivsU` : Function handle to the computation of the basis functions and their derivatives for the pair of primal unknown fields $\\left( w, \\beta \\right)$
%[text]  `computeBasisFunctionsAndDerivsBeta` : Function handle to the computation of the basis functions and their derivatives for the pair of the secondary unknown fields $\\left( Q, M \\right)$ - Only needed for the Hellinger-Reissner formulation of the Timoshenko beam
%[text]  `propStr` : Structure containing the following fields:
%[text]                                                                           .`A` : Cross sectional area of the beam
%[text]                                                                           .`q` : Distributed load on the beam
%[text]                                                                           .`m` : Distributed moment on the beam
%[text]                                                                           .`E` : Young's modulus
%[text]                                                                         .`nu` : Poisson's ratio
%[text]                                                                           .`I` : Moment of inertia
%[text]                                                                           .`G` : Shear modulus
%[text]                                                                   .`alpha` : Shear correction factor
%[text] 
%[text]  **Output :**
%[text]  `K` : Master stiffness matrix of a Timoshenko beam
%[text]  `F` : Master force vector of a Timoshenko beam
%[text]  `invHG` : The mapping matrix $\\mathbf{M}^{-1} \\mathbf{G}$ regarding the static condensation of the pair of secondary unknown fields $\\left( Q, M\\right)$ when using the Hellinger-Reissner principle, namely, $\\mathbf{\\beta} = \\mathbf{H}^{-1} \\mathbf{G} \\, \\mathbf{u}$ - Needed to compute the DOFs of the secondary unknown fields $\\left( Q, M\\right)$ given the DOFs of the primal unknown fields $\\left( w, \\beta \\right)$
%[text] %[text:anchor:H_4D5E6417] ## Function Implementation
%[text] %[text:anchor:H_9A915165] ### Input validation
    arguments
        msh (1, 1) struct {mustHaveNodesAndElements}
        computeElStiffMtxForceVct (1, 1) function_handle
        computeBasisFunctionsAndDerivsU (1, 1) function_handle
        computeBasisFunctionsAndDerivsBeta
        propStr (1, 1) struct {mustHaveTimoshenkoBeamProperties}
    end
%[text] %[text:anchor:H_E1EFAA4F] ### Read input
    numEl = length(msh.elements(:, 1));
    numNodes_u = length(msh.nodes(:, 1));
    numDOFs_u = 2*numNodes_u;
    isHR = false;
    if isa(computeBasisFunctionsAndDerivsBeta, 'function_handle')
        isHR = true;
    end
%[text] %[text:anchor:H_2C4267C4] ### Initialization of the global master stiffness matrix and master load vector
    K = zeros(numDOFs_u);
    F = zeros(numDOFs_u, 1);
    if isHR
        numNodes_beta_el = height(computeBasisFunctionsAndDerivsBeta(0));
        if numNodes_beta_el > 1
            error("The output variable 'invHG' is only valid for piecewise constant " + ...
                "approximation of the Lagrange Multipliers");
        end
        numNodes_beta = numNodes_beta_el*numEl;
        numDOFs_beta = 2*numNodes_beta;
        invHG = zeros(numDOFs_beta, numDOFs_u);
    else
        invHG = 'undefined';
    end
%[text] %[text:anchor:H_100CBAA1] ### Loop over all the elements in the mesh
    for ii = 1:numEl        
%[text] %[text:anchor:H_E457EA46] #### Nodal indices and coordinates
        id_nodes_el = msh.elements(ii, :);
        X = msh.nodes(id_nodes_el);
        propStr.idEl = ii; % Record the element numbering only for unit testing purposes
%[text] %[text:anchor:H_D22B23C7] #### Element stiffness matrix
        [Ke, Fe, invHeGe] = computeElStiffMtxForceVct ...
            (X, computeBasisFunctionsAndDerivsU, computeBasisFunctionsAndDerivsBeta, propStr);
%[text] %[text:anchor:H_CC777245] #### Element Freedom Tables (EFTs)
        EFTU = [2*id_nodes_el-1; 2*id_nodes_el];
        EFTU = EFTU(:);
        if isHR
            EFTBeta = [2*ii - 1; 2*ii]; % Because a constant approximation of the Lagrange Multipliers per element is chosen
        end
%[text] %[text:anchor:H_9D09FAC7] #### Assembly
        K(EFTU, EFTU) = K(EFTU, EFTU) + Ke;
        F(EFTU) = F(EFTU) + Fe;
        if isHR
            invHG(EFTBeta, EFTU) = invHG(EFTBeta, EFTU) + invHeGe;
        end
%[text] 
    end
%[text] 
end

%[appendix]{"version":"1.0"}
%---
