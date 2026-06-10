%[text] %[text:anchor:T_4AD0C9DA] # Master Stiffness Matrix and Body Load Vector for the Reissner-Mindlin Plate
%[text] Computation of the master stiffness matrix and master load vector corresponding to the standard Reissner-Mindlin plate formulation with three degrees of freedom per node $\\left(w,\\beta\_x ,\\beta\_y \\right)${"editStyle":"visual"}, the vertical deflection and two rotations of the plate's cross section around the $X${"editStyle":"visual"}- and the $Y${"editStyle":"visual"}-axis, respectively. The master stiffness matrix and master load vector for the standard Reissner-Mindlin plate are given by the following relations:
%[text] $\\begin{array}{ll}\n\\mathbf{K} = \\int\_{A} \\mathbf{B}\_{\\text{s}}^{\\text{T}} \\mathbf{C}\_{\\text{s}} \\mathbf{B}\_{\\text{s}} + \\int\_{A} \\mathbf{B}\_{\\text{b}}^{\\text{T}} \\mathbf{C}\_{\\text{b}} \\mathbf{B}\_{\\text{b}} \\: \\text{d} A , & \\text{(1.1)} \\\\\n\\mathbf{F} = \\int\_{A} \\mathbf{N}^{\\text{T}} \\left\[ \\begin{array}{c} \\bar{p} \\\\ \\bar{m}\_x \\\\ \\bar{m}\_y \\end{array} \\right\] \\: \\text{d} A + \\int\_{\\subset\\partial A} \\mathbf{N}^{\\text{T}} \\left\[ \\begin{array}{c} \\bar{P} \\\\ \\bar{M}\_x \\\\ \\bar{M}\_y \\end{array} \\right\] \\: \\text{d} \\partial A . & \\text{(1.2)}\n\\end{array}$
%[text:tableOfContents]{"heading":"**Table of Contents**"}
function [K, F] = computeMasterStiffMatrixandForceVectorReissnerMindlinPlate ...
    (msh, computeBasisFunctionsAndDerivs, computeStiffMatrixandForceVector, propStr)
%[text] %[text:anchor:H_9EECB707] ## **Function Description**
%[text] Returns the master stiffness matrix corresponding to a Reissner-Mindlin plate.
%[text]   **Input :**
%[text]  `msh` : Structure containing the nodes and the elements in a quadrilateral finite element mesh
%[text]                                                              .`nodes` : Array with the nodal coordinates in the finite element mesh
%[text]                                                         .`elements` : Array with the nodal numbering per element in the finite element mesh
%[text]  `computeBasisFunctionsAndDerivs` : Either a function handle to the computation of the basis functions and their derivatives or a cell array with the symbolic expression of the basis functions
%[text] `computeStiffMatrixandForceVector` : Function handle to the computation of the element stiffness matrix
%[text]  `propStr` : Structure containing the following fields,
%[text]                                                                       .`t` : Thickness
%[text]                                                                       .`q` : Distributed load
%[text]                                                                     .`mx` : Distributed moment around x-axis
%[text]                                                                     .`my` : Distributed moment around y-axis
%[text]                                                                       .`E` : Young's modulus
%[text]                                                                     .`nu` : Poisson's ratio
%[text]                                                                       .`G` : Shear modulus
%[text]                                                                       .`D` : Plate's stiffness
%[text]                                                               .`alpha` : Shear correction factor
%[text]  **Output :**
%[text]                                                              K : Master stiffness matrix of a Reissner-Mindlin plate
%[text]                                                              F : Master force vector of a Reissner-Mindlin plate
%[text] %[text:anchor:H_A79ECC99] ## **Function Implementation**
%[text] %[text:anchor:H_9016D179] ### Input validation
arguments
    msh (1, 1) struct {mustHaveNodesAndElements}
    computeBasisFunctionsAndDerivs {mustBeFunctionHandleOrSymfunCellArray}
    computeStiffMatrixandForceVector (1, 1) function_handle
    propStr (1, 1) struct {mustHaveReissnerMindlinPlateProperties}
end
%[text] %[text:anchor:H_45CC665F] ### Read input
    numEl = length(msh.elements(:, 1));
    numNodes = length(msh.nodes(:, 1));
    numDOFs = 3*numNodes;
%[text] %[text:anchor:H_87730FEB] ### Initialization of the global master stiffness matrix and master load vector
    K = zeros(numDOFs);
    F = zeros(numDOFs, 1);
%[text] %[text:anchor:H_FAA22C07] ### Loop over all the elements in the mesh
    for ii = 1:numEl    
%[text] %[text:anchor:H_8F9FF9F0] #### Get nodal indices and nodal coordinates
        id_nodes_el = msh.elements(ii, :);
        X = msh.nodes(id_nodes_el, 1:2);
        propStr.idEl = ii; % Record the element numbering only for unit testing purposes
%[text] %[text:anchor:H_6137285B] #### Element stiffness matrix
        [Ke, Fe] = computeStiffMatrixandForceVector ...
            (X, computeBasisFunctionsAndDerivs, propStr);
%[text] %[text:anchor:H_B06738C8] #### Element Freedom Table (EFT)
        EFT = [3*id_nodes_el-2; 3*id_nodes_el-1; 3*id_nodes_el];
        EFT = EFT(:);
%[text] %[text:anchor:H_60C11ACD] #### Assembly
        K(EFT, EFT) = K(EFT, EFT) + Ke;
        F(EFT) = F(EFT) + Fe;
%[text] 
    end
%[text] 
end
%%
%[text] %[text:anchor:H_17333DFC] ### Custom validation function
function mustBeFunctionHandleOrSymfunCellArray(inp)
    if ~isa(inp, 'function_handle') && ...
            ~isa(inp, 'cell')
        error("Input variable %s should be either a function handle or " + ...
            "a cell array of symbolic functions", inputname(1))
    end
end

%[appendix]{"version":"1.0"}
%---
