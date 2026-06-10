%[text] %[text:anchor:T_364076B1] # Element Stiffness Matrix and Body Load Vector for the Reissner-Mindlin Plate using selective-reduced integration
%[text] Returns the element stiffness matrix and body force load vector corresponding to the Reissner-Mindlin plate formulation with the standard approach but selective-reduced integration is used for the transverse-shear part of the virtual work. Since the polynomial order of the bending part of the virtual work is one-order less than for the transverse-shear part of the virtual work, one can reduce the number of Gauss points equally for both parts. Herein the separate Gauss integration for the two parts of the virtual work is demonstrated for the sake of generality. The master stiffness matrix and master load vector for the standard Reissner-Mindlin plate are given by the following relations:
%[text] $\\begin{array}{ll}\n\\mathbf{K}^e = \\int\_{A^e} \\left( \\mathbf{B}\_{\\text{s}}^e\\right)^{\\text{T}} \\mathbf{C}\_{\\text{s}} \\mathbf{B}\_{\\text{s}}^e + \\int\_{A^e} \\left( \\mathbf{B}\_{\\text{b}}^e\\right)^{\\text{T}} \\mathbf{C}\_{\\text{b}} \\mathbf{B}\_{\\text{b}}^e \\: \\text{d} A , & \\text{(1.1)} \\\\\n\\mathbf{F} = \\int\_{A^e} \\left( \\mathbf{N}^e\\right)^{\\text{T}} \\left\[ \\begin{array}{c} \\bar{p} \\\\ \\bar{m}\_x \\\\ \\bar{m}\_y \\end{array} \\right\] \\: \\text{d} A + \\int\_{\\subset\\partial A^e} \\left( \\mathbf{N}^e\\right)^{\\text{T}} \\left\[ \\begin{array}{c} \\bar{P} \\\\ \\bar{M}\_x \\\\ \\bar{M}\_y \\end{array} \\right\] \\: \\text{d} \\partial A . & \\text{(1.2)}\n\\end{array}$
%[text] In the frame of the selective-reduced integration the first term Eq. (1.1), namely, the transverse-shear part of the stiffness, is integrated using one less Gauss point than the second term in Eq. (1.1), namely, the bending part of the stiffness.
%[text:tableOfContents]{"heading":"**Table of Contents**"}
function [Ke, Fe] = computeElementStiffMatrixandForceVectorReissnerMindlinPlateRI ...
    (X, computeBasisFunctionsAndDerivs, propStr)
%[text] %[text:anchor:H_B4ABEB49] ## Function Description
%[text] Computation of the element stiffness matrix for an element based on the 2D Reissner-Mindlin plate theory with three degrees of freedom per node $\\left(w,\\beta\_x ,\\beta\_y \\right)${"editStyle":"visual"}, the vertical deflection and two rotations of the plate's cross section around the X- and the Y-axis, respectively. It is used selective-reduced integration for the transverse-shear part of the virtual work.
%[text] %[text:anchor:H_5E6E8C79]  **Input :**
%[text]  `X` : Array containing the nodal coordinates of the nodes in the element
%[text] `computeBasisFunctionsAndDerivs` : Function handle to the computation of the basis functions
%[text]  `propStr` : Structure containing the following fields:
%[text]                                                                   .`t` : Thickness of the shell
%[text]                                                             .`pBar` : Distributed load on the shell
%[text]                                                           .`mxBar` : Distributed moment on the beam around x-axis
%[text]                                                           .`myBar` : Distributed moment on the beam around y-axis
%[text]                                                                   .`E` : Young's modulus
%[text]                                                                 .`nu` : Poisson's ratio
%[text]                                                                   .`G` : Shear modulus
%[text]                                                                   .`D` : Plate's stiffness
%[text]                                                           .`alpha` : Shear correction factor
%[text] 
%[text] %[text:anchor:H_354D7949]  **Output :**
%[text]  `Ke` : Element stiffness matrix of a Reissner-Mindlin shell
%[text]  `Fe` : Element force vector of a Reissner-Mindlin shell
%[text] 
%[text] %[text:anchor:H_728009A7] ## Function Implementation
%[text] %[text:anchor:H_B2611BC6] ### Input validation
    arguments
        X (:, 2) double
        computeBasisFunctionsAndDerivs (1, 1) function_handle
        propStr (1, 1) {mustHaveReissnerMindlinPlateProperties}
    end
%[text] %[text:anchor:H_A583777B] ### Check input
    if strcmp(func2str(computeBasisFunctionsAndDerivs), ...
            "computeBilinearBasisFunctionsAndFirstDerivatives")
        numNodesEl = 4;
    elseif strcmp(func2str(computeBasisFunctionsAndDerivs), ...
            "computeBiquadraticBasisFunctionsAndFirstDerivatives")
        numNodesEl = 9;
    else
        error("Input 'computeBasisFunctionsAndDerivs' is defined as %s but only " + ...
            "values 'computeBilinearBasisFunctionsAndFirstDerivatives' and " + ...
            "'computeBiquadraticBasisFunctionsAndFirstDerivatives' are supported", ...
            computeBasisFunctionsAndDerivs);
    end
    if numNodesEl ~= numel(X(:, 1))
        error("It is provided an element with %d nodes and the function handle %s " + ...
            "for the computation of the basis functions", numel(X(:, 1)), ...
            func2str(computeBasisFunctionsAndDerivs));
    end
%[text] %[text:anchor:H_5CFBEA59] ### Initialize output
    Ke = zeros(12);
    Fe = zeros(12, 1);
%[text] %[text:anchor:H_4CF8352A] ### Material matrix
    C = [propStr.alpha*propStr.G*propStr.t 0                                 0                    0                     0
         0                                 propStr.alpha*propStr.G*propStr.t 0                    0                     0
         0                                 0                                 propStr.D            propStr.nu*propStr.D  0
         0                                 0                                 propStr.nu*propStr.D propStr.D             0
         0                                 0                                 0                    0                     propStr.D/2*(1 - propStr.nu)];
%[text] %[text:anchor:H_CD73EF2A] ### Gauss points and weights for the integration of the bending part of the virtual work
    if strcmp(func2str(computeBasisFunctionsAndDerivs), ...
            "computeBilinearBasisFunctionsAndFirstDerivatives")
        numGPShear = 2;
    elseif strcmp(func2str(computeBasisFunctionsAndDerivs), ...
            "computeBiquadraticBasisFunctionsAndFirstDerivatives")
        numGPShear = 4;
    end
    [xiGP, GWxi] = getGaussPointsAndWeightsOverUnitDomain(numGPShear);
    [etaGP, GWeta] = getGaussPointsAndWeightsOverUnitDomain(numGPShear);
%[text] %[text:anchor:H_588718D3] ### Loop over all the Gauss points to integrand the bending part of the virtual work
    for jj = 1:numel(GWeta)
        for ii = 1:numel(GWxi)
%[text] %[text:anchor:H_8753B626] #### Basis functions at the Gauss point
            [dN, ~] = computeBasisFunctionsAndDerivs(xiGP(ii), etaGP(jj));
%[text] %[text:anchor:H_BF8DA08D] #### Shape function matrix at the Gauss point
%             Nmtx = [dN(1, 1) 0         0        dN(2, 1) 0         0        dN(3, 1) 0         0        dN(4, 1) 0         0
%                     0        dN(1, 1)  0        0        dN(2, 1)  0        0        dN(3, 1)  0        0        dN(4, 1)  0
%                     0        0         dN(1, 1) 0        0         dN(2, 1) 0        0         dN(3, 1) 0        0         dN(4, 1)];
            Nmtx = zeros(3, 3*numNodesEl);
            for kk = 1:numNodesEl
                Nmtx(1, 3*kk - 2) = dN(kk, 1);
                Nmtx(2, 3*kk - 1) = dN(kk, 1);
                Nmtx(3, 3*kk) = dN(kk, 1);
            end
%[text] %[text:anchor:H_C15A1044] #### Jacobian matrix at the Gauss point
            J = zeros(2);
            for kk = 1:numNodesEl
                J(1, 1) = J(1, 1) + dN(kk, 2)*X(kk, 1);
                J(1, 2) = J(1, 2) + dN(kk, 2)*X(kk, 2);
                J(2, 1) = J(2, 1) + dN(kk, 3)*X(kk, 1);
                J(2, 2) = J(2, 2) + dN(kk, 3)*X(kk, 2);
            end
%[text] %[text:anchor:H_5BA72743] #### Derivatives of the basis functions with respect to the physical space at the Gauss point
            dNdX = (J\dN(:, 2:3)')';
%[text] %[text:anchor:H_7868B3D4] #### B-operator matrix for the bending part of the virtual work at the Gauss point
%             Bmtx = [  zeros(2, 12)
%                       0          0           dNdX(1, 1) 0           0           dNdX(2, 1)  0          0           dNdX(3, 1)  0           0           dNdX(4, 1)
%                       0          -dNdX(1, 2) 0          0           -dNdX(2, 2) 0           0          -dNdX(3, 2) 0           0           -dNdX(4, 2) 0
%                       0          -dNdX(1, 1) dNdX(1, 2) 0           -dNdX(2, 1) dNdX(2, 2)  0          -dNdX(3, 1) dNdX(3, 2)  0           -dNdX(4, 1) dNdX(4, 2)];
            Bmtx = zeros(5, 3*numNodesEl);
            for kk = 1:numNodesEl
                Bmtx(3, 3*kk) = dNdX(kk, 1);
                Bmtx(4, 3*kk - 1) = -dNdX(kk,2);
                Bmtx(5, 3*kk - 1) = -dNdX(kk, 1);
                Bmtx(5, 3*kk) = dNdX(kk,2);
            end
            Bmtx(1:2, :) = 0*Bmtx(1:2, :);
%[text] %[text:anchor:H_33F77F57] #### Computation and assembly of the element stiffness matrix for the bending part of the virtual work and body load vector to the global counterparts
            Ke = Ke + transpose(Bmtx)*C*Bmtx*det(J)*GWxi(ii)*GWeta(jj);
            Fe = Fe + transpose(Nmtx)*[propStr.q; propStr.mx; propStr.my]*det(J)*GWxi(ii)*GWeta(jj);
%[text] 
        end
%[text] 
    end
%[text] %[text:anchor:H_F42CEC4E] ### Gauss points and weights for the integration of the transverse-shear part of the virtual work
    numGPShear = numGPBending - 1;
    [xiGP, GWxi] = getGaussPointsAndWeightsOverUnitDomain(numGPShear);
    [etaGP, GWeta] = getGaussPointsAndWeightsOverUnitDomain(numGPShear);
%[text] %[text:anchor:H_18976F4F] ### Loop over all the Gauss points to integrand the bending part of the virtual work using selective-reduced integration
    for jj = 1:numel(GWeta)
        for ii = 1:numel(GWxi)
%[text] %[text:anchor:H_2E435A4C] #### Basis functions at the Gauss point
            [dN, ~] = computeBasisFunctionsAndDerivs(xiGP(ii), etaGP(jj));
%[text] %[text:anchor:H_43AB4948] #### Jacobian matrix at the Gauss point
            J = zeros(2);
            for kk = 1:numNodesEl
                J(1, 1) = J(1, 1) + dN(kk, 2)*X(kk, 1);
                J(1, 2) = J(1, 2) + dN(kk, 2)*X(kk, 2);
                J(2, 1) = J(2, 1) + dN(kk, 3)*X(kk, 1);
                J(2, 2) = J(2, 2) + dN(kk, 3)*X(kk, 2);
            end
%[text] %[text:anchor:H_5DEDA60A] #### Derivatives of the basis functions with respect to the physical space at the Gauss point
            dNdX = (J\dN(:, 2:3)')';
%[text] %[text:anchor:H_8631DBA6] #### B-operator matrix for the transverse-shear part of the virtual work at the Gauss point
%             Bmtx = [  dNdX(1, 1) 0           dN(1, 1)   dNdX(2, 1)  0           dN(2, 1)    dNdX(3, 1) 0           dN(3, 1)    dNdX(4, 1)  0           dN(4, 1)
%                       dNdX(1, 2) -dN(1, 1)   0          dNdX(2, 2)  -dN(2, 1)   0           dNdX(3, 2) -dN(3, 1)   0           dNdX(4, 2)  -dN(4, 1)   0
%                       zeros(3, 12)];
            Bmtx = zeros(5, 3*numNodesEl);
            for kk = 1:numNodesEl
                Bmtx(1, 3*kk - 2) = dNdX(kk, 1);
                Bmtx(1, 3*kk) = dN(kk, 1);
                Bmtx(2, 3*kk - 2) = dNdX(kk,2);
                Bmtx(2, 3*kk - 1) = -dN(kk, 1);
            end
%[text] %[text:anchor:H_BEC5F1E4] #### Computation and assembly of the element stiffness matrix for the transverse-shear part of the virtual work to the global counterpart
            Ke = Ke + transpose(Bmtx)*C*Bmtx*det(J)*GWxi(ii)*GWeta(jj);
%[text] 
        end
%[text] 
    end
%[text] 
end

%[appendix]{"version":"1.0"}
%---
