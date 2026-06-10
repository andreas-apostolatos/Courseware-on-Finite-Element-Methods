%[text] %[text:anchor:T_5EBA6AD8] # Element Stiffness Matrix and Body Load Vector for the Timoshenko beam
%[text] Returns the element stiffness matrix and body force load vector corresponding to the Timoshenko beam formulation with the standard approach. The master stiffness matrix and master load vector for the standard Timoshenko beam are given by the following relations:
%[text] $\\begin{array}{ll}\n\\mathbf{K}^e = \\int\_{L^e} \\left(\\mathbf{B}\_{\\text{s}}^e\\right)^{\\text{T}} \\mathbf{C}\_{\\text{s}} \\mathbf{B}\_{\\text{s}}^e \\, \\text{d} L+ \\int\_{L^e} \\left(\\mathbf{B}\_{\\text{b}}^e\\right)^{\\text{T}} \\mathbf{C}\_{\\text{b}} \\mathbf{B}\_{\\text{b}}^e \\: \\text{d} L , & \\text{(1.1)} \\\\\n\\mathbf{F}^e = \\int\_{L^e} \\left(\\mathbf{N}^e\\right)^{\\text{T}} \\left\[ \\begin{array}{c} \\bar{p} \\\\ \\bar{m} \\end{array} \\right\] \\: \\text{d} L + \\left\[\\left( \\mathbf{N}^e\\right)^{\\text{T}} \\left\[ \\begin{array}{c} \\bar{P} \\\\ \\bar{M} \\end{array} \\right\]  \\right\]\_0^{L^e}. & \\text{(1.2)}\n\\end{array}$
%[text:tableOfContents]{"heading":"**Table of Contents**"}
function [Ke, Fe, invHeGe] = computeElementStiffMatrixandForceVectorTimoshenkoBeam ...
    (X, computeBasisFunctionsAndDerivsU, ~, propStr)
%[text] %[text:anchor:H_457B96B8] ### Function description
%[text] Computation of the element stiffness matrix for an element based on the 1D Timoshenko theory with two degrees of freedom per node $\\left(w,\\beta \\right)${"editStyle":"visual"} the vertical deflection and the rotations of the beam's cross section around the X-axis, respectively using full integration, that is, 2 Gauß points per element.
%[text]  **Input :**
%[text]  `X` : Array containing the nodal coordinates of the nodes in the element
%[text]  `computeBasisFunctionsAndDerivsU` : Function handle to the computation of the basis functions and their derivatives
%[text]  `propStr` : Structure containing the following fields,
%[text]                                                                           .`A` : Cross sectional area of the beam
%[text]                                                                     .`qBar` : Distributed load on the beam
%[text]                                                                     .`mBar` : Distributed moment on the beam
%[text]                                                                           .`E` : Young's modulus
%[text]                                                                         .`nu` : Poisson's ratio
%[text]                                                                           .`I` : Moment of inertia
%[text]                                                                           .`G` : Shear modulus
%[text]                                                                   .`alpha` : Shear correction factor
%[text] 
%[text]  **Output :**
%[text]  `Ke` : Element stiffness matrix of a Timoshenko beam
%[text]  `Fe` : Element force vector of a Timoshenko beam
%[text]  `invHG` : Dummy output
%[text] %[text:anchor:H_659956C4] ### Function implementation
%[text] %[text:anchor:H_731E177F] ### Input validation
    arguments
        X (:, 1) double
        computeBasisFunctionsAndDerivsU (1, 1) function_handle
        ~
        propStr (1, 1) struct {mustHaveTimoshenkoBeamProperties}
    end
%[text] %[text:anchor:H_9EC643B6] ### Check input
    if strcmp(char(computeBasisFunctionsAndDerivsU), ...
            "computeLinearBasisFunctionsAndFirstDerivatives")
        if height(X) ~= 2
            error("Linear shape functions have been selected but the element " + ...
                "appears not be two-noded");
        end
    elseif strcmp(char(computeBasisFunctionsAndDerivsU), ...
            "computeQuadraticBasisFunctionsAndFirstDerivatives")
        if height(X) ~= 3
            error("Quadratic shape functions have been selected but the element " + ...
                "appears not be three-noded");
        end
    else
        error("Only linear and quadratic one-dimensional elements can be considered");
    end
%[text] %[text:anchor:H_0226BA10] ### Initialize output
    numNodes = height(X);
    numDOFs = 2*numNodes;
    Ke = zeros(numDOFs);
    Fe = zeros(numDOFs, 1);
    invHeGe = 'undefined';
%[text] %[text:anchor:H_BC9CB8C5] ### Material matrix
    C = [propStr.alpha*propStr.G*propStr.A 0
         0                                 propStr.E*propStr.I];
%[text] %[text:anchor:H_9D718938] ### Gauss points and weights
    numGP = height(X);
    [xiGP, GWxi] = getGaussPointsAndWeightsOverUnitDomain(numGP);
%[text] %[text:anchor:H_06C890A6] ### Loop over all the Gauss points
    for ii = 1:height(GWxi)
%[text] %[text:anchor:H_3632F41E] #### Basis functions at the Gauss point
        dN = computeBasisFunctionsAndDerivsU(xiGP(ii));
%[text] %[text:anchor:H_2531F8B9] #### Basis function matrix at the Gauss point 
        Nmtx = zeros(2, numDOFs);
        for jj = 1:numNodes
            Nmtx(1, 2*jj - 1) = dN(jj, 1);
            Nmtx(2, 2*jj) = dN(jj, 1);
        end
%[text] %[text:anchor:H_781352D8] #### Jacobian matrix at the Gauss point 
        Jmtx = dN(:, 2)'*X;
%[text] %[text:anchor:H_CC8CD106] #### Derivatives of the basis functions in the physical space at the Gauss point
        dNdX = (Jmtx\dN(:, 2)')';
%[text] %[text:anchor:H_2601F3CC] #### B-operator matrix at the Gauss point
        Bmtx = zeros(2, numDOFs);
        for jj = 1:numNodes
            Bmtx(1, 2*jj - 1) = dNdX(jj);
            Bmtx(1, 2*jj) = dN(jj, 1);
            Bmtx(2, 2*jj) = dNdX(jj);
        end
%[text] %[text:anchor:H_A8973D87] #### Element stiffness matrix and body force vector at the Gauss point and assembly to the global counterparts
        Ke = Ke + transpose(Bmtx)*C*Bmtx*det(Jmtx)*GWxi(ii);
        Fe = Fe + transpose(Nmtx)*[propStr.qBar; propStr.mBar]*det(Jmtx)*GWxi(ii);
%[text] 
    end
%[text] 
end

%[appendix]{"version":"1.0"}
%---
