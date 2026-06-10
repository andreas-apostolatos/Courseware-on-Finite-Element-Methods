%[text] %[text:anchor:H_CAFF79D5] # Computation of the L2-norm of the error for a Timoshenko beam
%[text] This function returns the $L^2(0,L)$-norm of the error for a Timoshenko beam problem given a Finite Element mesh of a Timoshenko beam problem, its discrete solution in terms of the displacement/rotation pair $\\left( u\_h,\\beta\_h \\right)$ and an analytical expression for the solution. The $L^2(0,L)$-norm of a tensorial function $\\mathbf{f}: \\left(0,L \\right) \\rightarrow \\mathbb{R}$ is defined as $\\left \\| \\mathbf{f} \\right\\|\_{L^2\\left( 0,L \\right)} = \\left( \\int\_{0}^L \\sum\_i \\sum\_j \\ldots \\sum\_k \\, f\_{ij \\ldots k} \\, f\_{ij \\ldots k} \\, \\text{d} X \\right)^{\\frac{1}{2}}$
%[text:tableOfContents]{"heading":"**Table of Contents**"}
function [errL2U, errL2Beta] = computeL2ErrorTimoshenkoBeam ...
    (msh, uh, computeBasisFunctionsAndDerivs, wEx, betaEx, numGP)
%[text] %[text:anchor:H_C497B93D] ## Function Description
%[text] This function returns the $L^2(0, L)$-norm of the error corresponding to the solution in terms of the vertical deflection $w$ of a beam modelled with the Timoshenko beam theory, namely,
%[text] $\\left\\| w - w\_h \\right\\|\_{L^2(0,L)} = \\left( \\int\_0^L \\left( w - w\_h \\right)^2 \\; \\text{d} x \\right)^{1/2} \\\\\n\\left\\| \\beta - \\beta\_h \\right\\|\_{L^2(0,L)} = \\left( \\int\_0^L \\left( \\beta - \\beta\_h \\right)^2 \\; \\text{d} x \\right)^{1/2}$
%[text]  **Input :**
%[text]  `msh` : Struct-variable representing the Finite Element mesh containing the following information:
%[text]                                                                          .nodes : Two-dimensional array containing the Cartesian coordinates of each node in the mesh
%[text]                                                                     .elements : Two-dimensional array containing the global numbering of the nodes that belong to each element
%[text]  `uh` : Solution vector corresponding to the Finite Element solution
%[text] `computeBasisFunctionsAndDerivs` : Function handle for the computation of the basis functions and their derivatives
%[text]  `wEx` : Symbolic expression of the analytical solution in terms of the vertical deflection $w\\left(X\\right)${"editStyle":"visual"} for the problem at hand
%[text]  `betaEx` : Symbolic expression of the analytical solution in terms of the cross-sectional rotation $\\beta \\left(X\\right)${"editStyle":"visual"} for the problem at hand
%[text]  `numGP` : Number of Gauss points for the evaluation of the integral regarding the $L^2$-norm
%[text]  **Output :**
%[text]  `errL2U` : The error of the vertical deflection field in the $L^2$-norm
%[text]  `errL2Beta` : The error of the cross-sectional rotation field in the $L^2$-norm
%[text] %[text:anchor:H_7DD4E5D1] ## Function Implementation
%[text] %[text:anchor:H_D4ED1162] ### Input validation
    arguments
        msh (1, 1) struct {mustHaveNodesAndElements}
        uh (:, 1) double
        computeBasisFunctionsAndDerivs (1, 1) function_handle
        wEx (1, 1) symfun
        betaEx (1, 1) symfun
        numGP (1, 1) double {mustBeInteger mustBeGreaterThanOrEqual(numGP, 1)}
    end
%[text] %[text:anchor:H_9270E4AF] ### Read input
    numEl = numel(msh.elements(:, 1));
%[text] %[text:anchor:H_9213DAB7] ## Initialize outputs
    errL2U = 0;
    errL2Beta = 0;
%[text] %[text:anchor:H_CDF10458] ## Loop over all the elements in the Finite Element mesh and compute the integral of the L2-norm
    for ii = 1:numEl
%[text] %[text:anchor:H_7BE60882] ### Element nodes by id
        id = msh.elements(ii, :);
%[text] %[text:anchor:H_61AF98FE] ### Element Freedom Table (EFT)
        EFT = vertcat(2*id - 1, 2*id);
        EFT = EFT(:);
%[text] %[text:anchor:H_E21494E5] ### Nodal discrete solution vector
        uEl = uh(EFT);
%[text] %[text:anchor:H_BC3BD820] ### Gauss points for the integration in the element level
        [GP, GW] = getGaussPointsAndWeightsOverUnitDomain(numGP);
%[text] %[text:anchor:H_4002E7C4] ### Loop over all the Gauss points and compute the error in the element level
        for jj = 1:height(GW)
%[text] %[text:anchor:H_B0ED1CE7] #### Basis functions at the Gauss point
            dN = computeBasisFunctionsAndDerivs(GP(jj));
            if numel(dN(:, 1)) ~= numel(id)
                error("The number of nodes %d and basis functions %d given the function handle %s do not match", ...
                    numel(id), numel(dN(:, 1)), func2str(computeBasisFunctionsAndDerivs));
            end
%[text] %[text:anchor:H_EA6D7FFB] #### Jacobian of the geometric transformation
            Jmtx = transpose(dN(:, 2))*msh.nodes(id);
%[text] %[text:anchor:H_1D0E61B9] #### Physical image of the Gauss point
            xCoord = transpose(dN(:, 1))*msh.nodes(id);
%[text] %[text:anchor:H_E11E1EED] #### Analytical solutions at the Gauss point
            wExVal = double(wEx(xCoord));
            betaExVal = double(betaEx(xCoord));
%[text] %[text:anchor:H_C45A6E79] #### Discrete numerical solutions at the Gauss point
            wh = transpose(dN(:, 1))*uEl(1:2:end);
            betah = transpose(dN(:, 1))*uEl(2:2:end);
%[text] %[text:anchor:H_7A66C9E0] #### Gauss point contribution to the L2-norm
            errL2U = errL2U + (wExVal - wh)^2*det(Jmtx)*GW(jj);
            errL2Beta = errL2Beta + (betaExVal - betah)^2*det(Jmtx)*GW(jj);
%[text] 
        end
%[text] 
    end
%[text] %[text:anchor:H_1049D8A9] ## Compute the square roots necessary for the L2-norm expression
    errL2U = sqrt(errL2U);
    errL2Beta = sqrt(errL2Beta);
%[text] 
end

%[appendix]{"version":"1.0"}
%---
