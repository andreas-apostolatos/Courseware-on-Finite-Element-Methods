%[text] %[text:anchor:T_9D6E7E0C] # Element Stiffness Matrix and Body Load Vector for the Timoshenko beam using the Hellinger-Reissner Principle
%[text] Returns the element stiffness matrix for an element based on the Hellinger-Reissner formulation of the 1D Timoshenko theory with four degrees of freedom per node $\\left(w,\\beta ,Q,M\\right)${"editStyle":"visual"}, the vertical deflection, the rotation of the beam's cross section around the X-axis, the reaction transverse shear forces and the reaction bending moments, respectively, using full integration. The pair of unknowns for the secondary unknown fields $\\left(Q,M\\right)${"editStyle":"visual"} is condensed out from the equation system. Assuming that the vector of DOFs related to the primal unknown fields $\\left(w,\\beta \\right)${"editStyle":"visual"} is denoted as $\\mathbf{u}$ whereas the vector of DOFs associated with the secondary unknown fields $\\left(Q,M\\right)${"editStyle":"visual"} is denoted by $\\mathbf{\\beta}$, the discrete equation system regarding the Hellinger-Reissner formulation of the Timoshenko beam problem reads as follows:
%[text] $\\left\[ \\begin{array}{cc} -\\mathbf{H} & \\mathbf{G} \\\\ \n\\mathbf{G}^{\\text{T}} & \\mathbf{0} \\end{array}\\right\] \\left\[ \\begin{array}{c} \\mathbf{\\beta} \\\\ \\mathbf{u} \\end{array} \\right\] = \\left\[ \\begin{array}{c} \\mathbf{0} \\\\ \\mathbf{F} \\end{array} \\right\]$
%[text] where,
%[text] $\\mathbf{M}^e = \\int\_{e} \\mathbf{N}\_{\\mathbf{\\beta}}^{\\text{T}} \\, \\mathbf{S} \\, \\mathbf{N}\_{\\mathbf{\\beta}} \\; \\text{d} e$,
%[text] $\\mathbf{G}^e = \\int\_{e}  \\mathbf{N}\_{\\mathbf{\\beta}}^{\\text{T}} \\, \\mathbf{B}\_{\\mathbf{\\mathbf{u}}} \\; \\text{d} e$.
%[text] B-operator matrix $\\mathbf{B}\_{\\mathbf{\\beta}}$ is related to the discretization of the pair $\\left(Q,M\\right)${"editStyle":"visual"} and it is given by the following equation,
%[text] $\\mathbf{N}\_{\\mathbf{\\beta}} = \\left\[ \\begin{array}{ccccc} N\_1^{\\mathbf{\\beta}} & 0 & \\ldots & N\_{n\_{\\mathbf{\\beta}}}^{\\mathbf{\\beta}} & 0 \\\\ \n0 & N\_1^{\\mathbf{\\beta}} & \\ldots & 0 & N\_{n\_{\\mathbf{\\beta}}}^{\\mathbf{\\beta}}\\end{array} \\right\]$,
%[text] since derivatives on the secondary pair of unknown fields do not appear, also piecewise constant functions can be used as basis functions. B-operator matrix $\\mathbf{B}\_{\\mathbf{u}}$ is related to the discretization of the shear deformation $\\gamma = \\frac{\\text{d} w}{\\text{d} x} + \\beta$ and it has been introduced already in the standard variational formulation of the Timoshenko beam. Finally, the compliance matrix $\\mathbf{S}$ is given by the following relationship,
%[text] $\\mathbf{S} = \\left\[ \\begin{array}{cc} \\frac{1}{\\alpha G A} & 0 \\\\ 0 & \\frac{1}{EI} \\end{array} \\right\]$, which is the inverse of the material matrix in this case.
%[text:tableOfContents]{"heading":"**Table of Contents**"}
%[text] %[text:anchor:H_1654FB05] ## Function implementation
function [Ke, Fe, invHeGe] = ... 
    computeElStiffMtxForceVctHellingerReissnerTimoshenkoBeam ...
    (X, computeBasisFunctionsAndDerivsU, computeBasisFunctionsAndDerivsBeta, propStr)
%[text] %[text:anchor:H_B0A7A0B3] ## Function Description
%[text] Returns the element stiffness matrix for an element based on the Hellinger-Reissner formulation of the 1D Timoshenko theory with four degrees of freedom per node $\\left(w,\\beta ,Q,M\\right)${"editStyle":"visual"}, the vertical deflection, the rotation of the beam's cross section around the X-axis, the reaction transverse shear forces and the reaction bending moments, respectively, using full integration.
%[text] %[text:anchor:H_D5E11DB6]  **Input** :
%[text]  `X` : Array containing the nodal coordinates of the nodes in the element
%[text]  `computeBasisFunctionsAndDerivsU` : Function handle to the computation of the basis functions and their derivatives for the primary pair of unknown fields $\\left( w, \\beta \\right)$
%[text] `computeBasisFunctionsAndDerivsBeta` : Function handle to the computation of the basis functions and their derivatives for the secondary pair of unknown fields $\\left( Q, M \\right)$
%[text]  `propStr` : Structure containing the following fields,
%[text]                                                                            .`A` : Cross sectional area of the beam
%[text]                                                                      .`qBar` : Distributed load on the beam
%[text]                                                                      .`mBar` : Distributed moment on the beam
%[text]                                                                            .`E` : Young's modulus
%[text]                                                                          .`nu` : Poisson's ratio
%[text]                                                                            .`I` : Moment of inertia
%[text]                                                                            .`G` : Shear modulus
%[text]                                                                    .`alpha` : Shear correction factor
%[text] 
%[text] %[text:anchor:H_25E71D8D]  **Output** :
%[text]  `Ke` : Element stiffness matrix of a Timoshenko beam
%[text]  `Fe` : Element force vector of a Timoshenko beam
%[text]  `invHeGe` : The mapping matrix $\\mathbf{M}^{-1} \\mathbf{G}$ regarding the static condensation of the pair of secondary unknown fields $\\left( Q, M\\right)$, namely, $\\mathbf{\\beta} = \\mathbf{H}^{-1} \\mathbf{G} \\, \\mathbf{u}$ - Needed to compute the DOFs of the secondary unknown fields $\\left( Q, M\\right)$ given the DOFs of the primal unknown fields $\\left( w, \\beta \\right)$
%[text] %[text:anchor:H_A140D0A3] ## Function Implementation
%[text] %[text:anchor:H_731E177F] ### Input validation
    arguments
        X (:, 1) double
        computeBasisFunctionsAndDerivsU (1, 1) function_handle
        computeBasisFunctionsAndDerivsBeta (1, 1) function_handle
        propStr (1, 1) struct {mustHaveTimoshenkoBeamProperties}
    end
%[text] %[text:anchor:H_70CBBDFC] ### Check input
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
%[text] %[text:anchor:H_C30984B8] ### Initialize output
    numNodesU = height(X);
    numDOFsU = 2*numNodesU;
    numNodesBeta = height(computeBasisFunctionsAndDerivsBeta(0));
    numDOFsBeta = 2*numNodesBeta;
    He = zeros(numDOFsBeta);
    Ge = zeros(numDOFsBeta, numDOFsU);
    Fe = zeros(numDOFsU, 1);
%[text] %[text:anchor:H_CC557936] ### Compliance matrix
    S = [1/(propStr.alpha*propStr.G*propStr.A) 0
         0                                     1/(propStr.E*propStr.I)];
%[text] %[text:anchor:H_5CC03240] ### Gauss points and weights
    numGP = height(X);
    [xiGP, GWxi] = getGaussPointsAndWeightsOverUnitDomain(numGP);
%[text] %[text:anchor:H_D2C850A6] ### Loop over all the Gauss points
    for ii = 1:height(GWxi)
%[text] %[text:anchor:H_9809A1C8] #### Basis functions at the Gauss point        
        dNU = computeBasisFunctionsAndDerivsU(xiGP(ii));
        dNBeta = computeBasisFunctionsAndDerivsBeta(xiGP(ii));
%[text] %[text:anchor:H_838B6673] #### Basis function matrix at the Gauss point
        Nmtx = zeros(2, numDOFsU);
        for jj = 1:numNodesU
            Nmtx(1, 2*jj - 1) = dNU(jj, 1);
            Nmtx(2, 2*jj) = dNU(jj, 1);
        end
%[text] %[text:anchor:H_80769174] #### Jacobian matrix at the Gauss point
        Jmtx = transpose(dNU(:, 2))*X;
%[text] %[text:anchor:H_7CCA12B9] #### Derivatives of the basis functions for the pair of primal unknown fields in the physical space at the Gauss point
%[text] $\\left\[ \\begin{array}{c} \\frac{\\text{d} w}{ \\text{d} X} \\\\ \\frac{\\text{d} \\beta}{ \\text{d} X} \\end{array} \\right\] = \\mathbf{J}^{-T} \\left\[ \\begin{array}{c} \\frac{\\text{d} w}{ \\text{d} \\xi} \\\\ \\frac{\\text{d} \\beta}{ \\text{d} \\xi} \\end{array} \\right\]$
        dNdX = transpose((Jmtx\transpose(dNU(:, 2))));
%[text] %[text:anchor:H_475A2471] #### B-operator matrix for the strain vector associated with the primal unknown fields at the Gauss point
%[text] $\\left\[ \\begin{array}{c} \\gamma \\left( w, \\frac{\\text{d} w}{\\text{d} X}, \\beta, \\frac{\\text{d} \\beta}{ \\text{d} X}\\right) \\\\  \\kappa \\left( w, \\frac{\\text{d} w}{ \\text{d} X}, \\beta, \\frac{\\text{d} \\beta}{ \\text{d} X} \\right) \\end{array} \\right\] = \\mathbf{B}\_u \\, \\mathbf{u}$
        BmtxU = zeros(2, numDOFsU);
        for jj = 1:numNodesU
            BmtxU(1, 2*jj - 1) = dNdX(jj);
            BmtxU(1, 2*jj) = dNU(jj, 1);
            BmtxU(2, 2*jj) = dNdX(jj);
        end
%[text] %[text:anchor:H_22F62C9D] #### B-operator matrix for the pair of secondary unknown fields at the Gauss point
%[text] $\\left\[ \\begin{array}{c} \\frac{\\text{d} Q}{ \\text{d} X} \\\\ \\frac{\\text{d} M}{ \\text{d} X} \\end{array} \\right\] = \\mathbf{N}\_\\beta \\, \\mathbf{\\beta}$
        NmtxBeta = zeros(2, numDOFsBeta);
        for jj = 1:numNodesBeta
            NmtxBeta(1, 2*jj - 1) = dNBeta(jj, 1);
            NmtxBeta(2, 2*jj) = dNBeta(jj, 1);
        end
%[text] %[text:anchor:H_5FB5B609] #### Element matrices corresponding to the Hellinger-Reissner formulation and body force vector at the Gauss point and assembly to the global counterparts
        He = He + transpose(NmtxBeta)*S*NmtxBeta*det(Jmtx)*GWxi(ii);
        Ge = Ge + transpose(NmtxBeta)*BmtxU*det(Jmtx)*GWxi(ii);
        Fe = Fe + transpose(Nmtx)*[propStr.qBar; propStr.mBar]*det(Jmtx)*GWxi(ii);
%[text] 
    end
%[text] %[text:anchor:H_E9EC4493] ### Compute the part of the stiffness matrix associated purely with the vertical deflection - cross sectional rotation DOFs
    invHeGe = He\Ge;
    Ke = transpose(Ge)*invHeGe; % Simplification for the expression transpose(Ge)*inv(He)*Ge 
%[text] 
end

%[appendix]{"version":"1.0"}
%---
