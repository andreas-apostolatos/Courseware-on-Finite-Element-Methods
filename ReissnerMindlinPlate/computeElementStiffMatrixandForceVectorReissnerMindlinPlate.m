%[text] %[text:anchor:T_5EBA6AD8] # Element Stiffness Matrix and Body Load Vector for the Reissner-Mindlin Plate
%[text] Returns the element stiffness matrix and body force load vector corresponding to the Reissner-Mindlin plate formulation with the standard approach. The master stiffness matrix and master load vector for the standard Reissner-Mindlin plate are given by the following relations:
%[text] $\\begin{array}{ll}\n\\mathbf{K} = \\int\_{A} \\mathbf{B}\_{\\text{s}}^{\\text{T}} \\mathbf{C}\_{\\text{s}} \\mathbf{B}\_{\\text{s}} \\, \\text{d} A+ \\int\_{A} \\mathbf{B}\_{\\text{b}}^{\\text{T}} \\mathbf{C}\_{\\text{b}} \\mathbf{B}\_{\\text{b}} \\: \\text{d} A , & \\text{(1.1)} \\\\\n\\mathbf{F} = \\int\_{A} \\mathbf{N}^{\\text{T}} \\left\[ \\begin{array}{c} \\bar{p} \\\\ \\bar{m}\_x \\\\ \\bar{m}\_y \\end{array} \\right\] \\: \\text{d} A + \\int\_{\\subset\\partial A} \\mathbf{N}^{\\text{T}} \\left\[ \\begin{array}{c} \\bar{P} \\\\ \\bar{M}\_x \\\\ \\bar{M}\_y \\end{array} \\right\] \\: \\text{d} \\partial A . & \\text{(1.2)}\n\\end{array}$
%[text:tableOfContents]{"heading":"**Table of Contents**"}
function [Ke, Fe] = computeElementStiffMatrixandForceVectorReissnerMindlinPlate ...
    (X, computeBasisFunctionsAndDerivs, propStr)
%[text] %[text:anchor:H_B4ABEB49] ## Function Description
%[text] Computation of the element stiffness matrix for an element based on the 2D Reissner-Mindlin plate theory with three degrees of freedom per node $\\left(w,\\beta\_x ,\\beta\_y \\right)${"editStyle":"visual"}, the vertical deflection and two rotations of the plate's cross section around the X- and the Y-axis, respectively.
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
%[text] %[text:anchor:H_D9FE555F] ### Check input
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
%[text] %[text:anchor:H_4BA8D8E8] ### Initialize output
    Ke = zeros(3*numNodesEl);
    Fe = zeros(3*numNodesEl, 1);
%[text] %[text:anchor:H_55F17252] ### Material matrix
    C = [propStr.alpha*propStr.G*propStr.t 0                                 0                    0                     0
         0                                 propStr.alpha*propStr.G*propStr.t 0                    0                     0
         0                                 0                                 propStr.D            propStr.nu*propStr.D  0
         0                                 0                                 propStr.nu*propStr.D propStr.D             0
         0                                 0                                 0                    0                     propStr.D/2*(1 - propStr.nu)];
%[text] %[text:anchor:H_553CE5F6] ### Gauss points and weights
    if strcmp(func2str(computeBasisFunctionsAndDerivs), ...
            "computeBilinearBasisFunctionsAndFirstDerivatives")
        numGP = 2;
    elseif strcmp(func2str(computeBasisFunctionsAndDerivs), ...
            "computeBiquadraticBasisFunctionsAndFirstDerivatives")
        numGP = 4;
    end
    [xiGP, GWxi] = getGaussPointsAndWeightsOverUnitDomain(numGP);
    [etaGP, GWeta] = getGaussPointsAndWeightsOverUnitDomain(numGP);
%[text] %[text:anchor:H_3148F91D] ### Loop over all the Gauss points
    for jj = 1:height(GWeta)
        for ii = 1:height(GWxi)
%[text] %[text:anchor:H_10B028EF] #### Computation of the the basis functions at the Gauss point
            [dN, ~] = computeBasisFunctionsAndDerivs(xiGP(ii), etaGP(jj));
            if numel(X(:, 1)) ~= numel(dN(:, 1))
                error("Number of nodes does not match the number of basis functions");
            end
%[text] %[text:anchor:H_B995C283] #### Shape function matrix at the Gauss point
%             Nmtx = [dN(1, 1) 0         0        dN(2, 1) 0         0        dN(3, 1) 0         0        dN(4, 1) 0         0
%                     0        dN(1, 1)  0        0        dN(2, 1)  0        0        dN(3, 1)  0        0        dN(4, 1)  0
%                     0        0         dN(1, 1) 0        0         dN(2, 1) 0        0         dN(3, 1) 0        0         dN(4, 1)];
            Nmtx = zeros(3, 3*numNodesEl);
            for kk = 1:numNodesEl
                Nmtx(1, 3*kk - 2) = dN(kk, 1);
                Nmtx(2, 3*kk - 1) = dN(kk, 1);
                Nmtx(3, 3*kk) = dN(kk, 1);
            end
%[text] %[text:anchor:H_94D5143A] #### Jacobian matrix at the Gauss point
            J = zeros(2);
            for kk = 1:numNodesEl
                J(1, 1) = J(1, 1) + dN(kk, 2)*X(kk, 1);
                J(1, 2) = J(1, 2) + dN(kk, 2)*X(kk, 2);
                J(2, 1) = J(2, 1) + dN(kk, 3)*X(kk, 1);
                J(2, 2) = J(2, 2) + dN(kk, 3)*X(kk, 2);
            end

%[text] %[text:anchor:H_4A39B992] #### Derivatives of the basis functions with respect to the physical space at the Gauss point
            dNdX = transpose(J\transpose(dN(:, 2:3)));
%[text] %[text:anchor:H_1DE4C7A6] #### B-operator matrix at the Gauss point
%             Bmtx = [  dNdX(1, 1) 0           dN(1, 1)   dNdX(2, 1)  0           dN(2, 1)    dNdX(3, 1) 0           dN(3, 1)    dNdX(4, 1)  0           dN(4, 1)
%                       dNdX(1, 2) -dN(1, 1)   0          dNdX(2, 2)  -dN(2, 1)   0           dNdX(3, 2) -dN(3, 1)   0           dNdX(4, 2)  -dN(4, 1)   0
%                       0          0           dNdX(1, 1) 0           0           dNdX(2, 1)  0          0           dNdX(3, 1)  0           0           dNdX(4, 1)
%                       0          -dNdX(1, 2) 0          0           -dNdX(2, 2) 0           0          -dNdX(3, 2) 0           0           -dNdX(4, 2) 0
%                       0          -dNdX(1, 1) dNdX(1, 2) 0           -dNdX(2, 1) dNdX(2, 2)  0          -dNdX(3, 1) dNdX(3, 2)  0           -dNdX(4, 1) dNdX(4, 2)];
            Bmtx = zeros(5, 3*numNodesEl);
            for kk = 1:numNodesEl
                Bmtx(1, 3*kk - 2) = dNdX(kk, 1);
                Bmtx(1, 3*kk) = dN(kk, 1);
                Bmtx(2, 3*kk - 2) = dNdX(kk, 2);
                Bmtx(2, 3*kk - 1) = -dN(kk, 1);
                Bmtx(3, 3*kk) = dNdX(kk, 1);
                Bmtx(4, 3*kk - 1) = -dNdX(kk, 2);
                Bmtx(5, 3*kk - 1) = -dNdX(kk, 1);
                Bmtx(5, 3*kk) = dNdX(kk,2);
            end
%[text] %[text:anchor:H_7A8E7A6E] #### Computation of the element stiffness matrix and addition of the contribution from the Gauss points and body load vector
            Ke = Ke + transpose(Bmtx)*C*Bmtx*det(J)*GWxi(ii)*GWeta(jj);
            Fe = Fe + transpose(Nmtx)*[propStr.pBar; propStr.mxBar; propStr.myBar]* ...
                det(J)*GWxi(ii)*GWeta(jj);
%[text] 
        end
%[text] 
    end
%[text] 
end

%[appendix]{"version":"1.0"}
%---
