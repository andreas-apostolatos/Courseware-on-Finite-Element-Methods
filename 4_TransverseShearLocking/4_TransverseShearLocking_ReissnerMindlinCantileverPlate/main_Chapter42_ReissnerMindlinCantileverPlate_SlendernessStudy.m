%[text] %[text:anchor:H_AE8EB025] # Slenderness study on transverse shear-locking
%[text:tableOfContents]{"heading":"Table of Contents"} %[text:anchor:M_6CC50A24]
%[text] In this section it is investigated the behavior of the convergence as the slenderness (length over thickness) goes to infinity $L/t \\rightarrow \\infty$, namely, as the thickness goes to zero $t \\rightarrow 0$. As the thickness goes to zero, the spurious forces due to transverse shear locking increase, leading to an ever increasing artificially stiffening behavior. To keep the solution constant the same scaling factor is applied both to the applied force and the thickness, so that the ratio $p/t$, and thus the resulting deformation, remains constant.
%[text] 
%[text] Go to [THIS](file:.\main_Chapter4_Locking_ReissnerMindlinCantileverPlate.m) Section in the main driver script *Transverse Shear Locking in Shear-Deformable Plates* *`-`* *Cantilever Plate Subject to Vertical Distributed Load*
%[text] Go to [PREVIOUS](file:.\main_Chapter41_ReissnerMindlinCantileverPlate_ShearForces.m) Section *Transverse-shear locking due incompatible discretization spaces*
%[text] Go to [NEXT](file:.\main_Chapter43_ReissnerMindlinCantileverPlate_ConvergenceStudy.m) Section *Convergence study on transverse shear-locking*
%[text] 
%[text] **Author: Andreas Apostolatos, PhD (aapostol@mathworks.com)**
%[text] **Date: 10.06.2026**
%%
%[text] %[text:anchor:H_F203EBCC] ## Preamble
vars = ["propStr" "XLx" "YLy" "numEl" "mshQ1" "mshQ2" ...
    "numNodesElQ1" "freeDOFsQ1" "freeDOFsQ2" "wEx"];
str = "This m-file should be ran from parent m-file " + ...
    "main_Chapter4_Locking_ReissnerMindlinCantileverPlate.m";
for ii = 1:length(vars)
    if ~exist(vars(ii), 'var') 
        error(strcat("Variable ", vars(ii), " is not defined. ", str));
    end
end
%%
%[text] %[text:anchor:H_D7EDC448] ## Initialization of the resulting vertical deflection with respect to the slenderness $L/t$
propStrSlend = propStr; %[text:anchor:H_F1F3C605]
numPtsSlend = mPowU - mPowL + 1;
slenderness_ratio = 10.^linspace(mPowL, mPowU, numPtsSlend);
numSlendLevels = numel(slenderness_ratio);
tipDeflectionSlend = zeros(numSlendLevels, 2);
tipDeflectionAnalyticalSlend = zeros(numSlendLevels, 1);
%[text] %[text:anchor:H_F396834F] ## Initialize auxiliary array
syms wExSlend(X)
%[text] %[text:anchor:H_73C5CB82] ## Point for the displacement evaluation
xyPostProc = [XLx (0 + YLy)/2];
if ~isequal(size(xyPostProc), [1 2])
    error("The provided points should be a row-vector with two components")
end
%[text] %[text:anchor:H_E38F9D89] ## Element in mesh containing the evaluation point
elPstProc = nan;
for ii = 1:numEl
    [in, ~] = inpolygon(xyPostProc(1), xyPostProc(2), ...
        mshQ1.nodes(mshQ1.elements(ii, :), 1), ...
        mshQ1.nodes(mshQ1.elements(ii, :), 2));
    if in
        elPstProc = ii;
        break;
    end
end
if isnan(elPstProc)
    error("No element found containing the point where " + ...
        "the postprocessing results are evaluated");
end
propNR.eps = 1e-12;
propNR.maxIt = 20;
[xiEta, isConverged] = ...
    computePointCoordinatesOnCanonicalBilinearQuadrilateral ...
    (xyPostProc', mshQ1.nodes(mshQ1.elements(elPstProc, :), 1:2)', ...
    propNR);
if ~isConverged
    error("The parametric coordinates of point (%d, %d) can be found " + ...
        "in no element %d", xyPostProc(1), xyPostProc(2), elPstProc);
end
%[text] %[text:anchor:H_5CE3A2DD] ## Element Freedom Tables (EFTs) for the bilinear and biquadratic quadrilateral meshes for the element containing the evaluation point
EFTQ1PostprocSlend = [3*mshQ1.elements(elPstProc, :) - 2
                      3*mshQ1.elements(elPstProc, :) - 1
                      3*mshQ1.elements(elPstProc, :)];
EFTQ1PostprocSlend = EFTQ1PostprocSlend(:);
EFTQ2PostprocSlend = [3*mshQ2.elements(elPstProc, :) - 2
                      3*mshQ2.elements(elPstProc, :) - 1
                      3*mshQ2.elements(elPstProc, :)];
EFTQ2PostprocSlend = EFTQ2PostprocSlend(:);
%[text] %[text:anchor:H_C7F3552B] ## Basis functions and set-up the basis function matrices at the parametric location of the evaluation point
%[text] %[text:anchor:H_1662A777] ### Basis functions matrix for the bilinear quadrilateral mesh
ShapeFunctionValuesQ1 = ...
    computeBilinearBasisFunctionsAndFirstDerivatives ...
    (xiEta(1, 1), xiEta(2, 1));
NmtxQ1 = zeros(3, 3*numNodesElQ1);
for ii = 1:numel(ShapeFunctionValuesQ1(:, 1))
    NmtxQ1(1, 3*ii - 2) = ShapeFunctionValuesQ1(ii, 1);
    NmtxQ1(2, 3*ii - 1) = ShapeFunctionValuesQ1(ii, 1);
    NmtxQ1(3, 3*ii) = ShapeFunctionValuesQ1(ii, 1);
end
%[text] %[text:anchor:H_094BDDFB] ### Basis functions matrix for the biquadratic quadrilateral mesh
ShapeFunctionValuesQ2 = ...
    computeBiquadraticBasisFunctionsAndFirstDerivatives ...
    (xiEta(1, 1), xiEta(2, 1));
NmtxQ2 = zeros(3, 3*numNodesElQ2);
for ii = 1:numel(ShapeFunctionValuesQ2(:, 1))
    NmtxQ2(1, 3*ii - 2) = ShapeFunctionValuesQ2(ii, 1);
    NmtxQ2(2, 3*ii - 1) = ShapeFunctionValuesQ2(ii, 1);
    NmtxQ2(3, 3*ii) = ShapeFunctionValuesQ2(ii, 1);
end
%%
%[text] %[text:anchor:H_76F3F67A] ## Loop over all the slenderness levels
for ii = 1:numSlendLevels
%[text] %[text:anchor:H_A31D31E7] ### **Scaling of the externally applied load and the thickness of the plate based on the slederness to obtain dimensionless quantities**
    propStrSlend.t = propStr.t/slenderness_ratio(ii);
    propStrSlend.pBar = propStr.pBar*propStrSlend.t^3*XLx^3;
    propStrSlend.D = ...
        propStrSlend.E*propStrSlend.t^3/12/(1 - propStrSlend.nu^2);
%[text] %[text:anchor:H_EC1E2D66] ### Master stiffness matrix  corresponding to the bilinear quadrilateral mesh
    computeStiffMatrixandForceVectorQ1 = 'undefined';
    computeBasisFunctionsAndDerivsQ1 = ...
        @computeBilinearBasisFunctionsAndFirstDerivatives;
    [KQ1Sled, FQ1Sled] = ...
        computeMasterStiffMatrixandForceVctReissnerMindlinPlatePageWise ...
        (mshQ1, computeBasisFunctionsAndDerivsQ1, ...
        computeStiffMatrixandForceVectorQ1, propStrSlend);
%[text] ![](text:image:9e6a) **Try**
%[text] Compute the master stiffness matrix using reduced integration
%[text] ```matlabCodeExample
%[text] computeStiffMatrixandForceVectorRI = ...
%[text]     @computeElementStiffMatrixandForceVectorReissnerMindlinPlateRI;
%[text] [KRI, FRI] = ...
%[text]     computeMasterStiffMatrixandForceVectorReissnerMindlinPlate ...
%[text]     (msh, computeStiffMatrixandForceVectorRI, propStr);
%[text] ```
%[text] %[text:anchor:H_7FF41760] ### Master stiffness matrix  corresponding to the biquadratic quadrilateral mesh
    computeStiffMatrixandForceVectorQ2 = 'undefined';
    computeBasisFunctionsAndDerivsQ2 = ...
        @computeBiquadraticBasisFunctionsAndFirstDerivatives;
    [KQ2Sled, FQ2Sled] = ...
        computeMasterStiffMatrixandForceVctReissnerMindlinPlatePageWise ...
        (mshQ2, computeBasisFunctionsAndDerivsQ2, ...
        computeStiffMatrixandForceVectorQ2, propStrSlend);
%[text] %[text:anchor:H_C8F31CD1] ### Solution of the linear equation system corresponding to the bilinear quadrilateral mesh
    uQ1Sled = zeros(numDOFsQ1, 1);
    uQ1Sled(freeDOFsQ1, 1) = ... 
        KQ1Sled(freeDOFsQ1, freeDOFsQ1)\FQ1Sled(freeDOFsQ1);
%[text] ![](text:image:2360) **Try**
%[text] Solve the linear equation system corresponding to the stiffness matrix with reduced integration
%[text] ```matlabCodeExample
%[text] try
%[text]   u(freeDOFs, 2) = KRI(freeDOFs, freeDOFs)\FRI(freeDOFs);
%[text] catch ME
%[text]     rethrow(ME)
%[text] end
%[text] ```
%[text] %[text:anchor:H_1CA2EDD2] ### Solution of the linear equation system corresponding to the biquadratic quadrilateral mesh
    uQ2Sled = zeros(numDOFsQ2, 1);
    uQ2Sled(freeDOFsQ2, 1) = ...
        KQ2Sled(freeDOFsQ2, freeDOFsQ2)\FQ2Sled(freeDOFsQ2);
%[text] %[text:anchor:H_5806E697] ### Vertical deflection at the selected point corresponding to the bilinear quadrilateral mesh
    wPostprocQ1 = NmtxQ1*uQ1Sled(EFTQ1PostprocSlend);
    tipDeflectionSlend(ii, 1) = wPostprocQ1(1);
%[text] %[text:anchor:H_1A937D5F] ### Vertical deflection at the selected point corresponding to the biquadratic quadrilateral mesh
    wPostprocQ2 = NmtxQ2*uQ2Sled(EFTQ2PostprocSlend);
    tipDeflectionSlend(ii, 2) = wPostprocQ2(1);
%[text] %[text:anchor:H_82F5AAB4] ### Analytical solution for the vertical deflection field
    wExSlend(X) = subs(wEx,{Lsym, psym, Gsym, Esym, Isym, Aqsym}, ...
        {XLx, propStrSlend.pBar, propStrSlend.G, propStrSlend.E, ...
        YLy*(propStrSlend.t^3)/12, propStrSlend.alpha*YLy*propStrSlend.t});
    tipDeflectionAnalyticalSlend(ii, 1) = double(wExSlend(XLx));
%[text] 
end
%%
%[text] %[text:anchor:H_EC3F1B32] ## Plot the results of the slenderness study
loglog(slenderness_ratio, abs(tipDeflectionSlend(:, 1)), "-bo", LineWidth=2) %[text:anchor:M_A0FC8936]
hold on;
loglog(slenderness_ratio, abs(tipDeflectionSlend(:, 2)), "-ro", LineWidth=2)
loglog(slenderness_ratio, abs(tipDeflectionAnalyticalSlend), "-k", LineWidth=2)
hold off;
grid on;
legend("Q1", "Q2", "analytical", Location="best")
xlabel("Slenderness $L/t$", interpreter='latex')        
ylabel("Tip displacement")
title("Slenderness study")
%[text] ![](text:image:91f0) **Try**
%[text] Revise the vector containing the values of the analytical solution throughout the slenderness levels saved in variable `tip_deflection_analytical_slend`
%[text] ```matlabCodeExample
%[text] loglog(slenderness_ratio, abs(tipDeflectionSlend(:, 1)), "-bo", ...
%[text]     slendernessRatio, abs(tipDeflectionSlend(:, 2)), "-ro", ...
%[text]     slendernessRatio, abs(tipDeflectionAnalyticalSlend(:, 1)), "-ko");
%[text] legend("Q1", "Q2", "Analytical")
%[text] xlabel("Slenderness L/t")
%[text] ylabel("Tip displacement")
%[text] title("Slenderness level")
%[text] ```
%[text] ![](text:image:1408) **Reflect**
%[text] How does the analytical solution change with respect to the different slenderness levels and why?

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"inline","rightPanelPercent":40}
%---
%[text:image:9e6a]
%   data: {"align":"baseline","height":26,"src":"data:image\/png;base64,iVBORw0KGgoAAAANSUhEUgAAAB4AAAApCAYAAAAiT5m3AAADOElEQVR42r1XSWgUQRSNGndR3Ii7V\/EmCCp4DaIe1Oh4EMHBwMjsq47Eg+OOJ0VIICAeFAIeXEADHgyoqChE1ASNhISAxH0jg5gxq+\/FX1A009PTM93T8OiuX7\/q1f+\/+tevqqoCTyKRWBmNRo8Cj4GPwF++I5HII7zT8Xh8eZWTTzgcnh6Lxc5j8kFgHMgCT4G7wBPgl8j\/YBGnfT7f1LJJg8HgQrGQE3diAXVciK5DIsi3o7+devh+4Pf755dMmslkpmGih0LaiHa1hf5k6J0T\/TYrfdMHKz8pFjTr1kFWD9ltvJ8Bt\/Dt1d0rYeG4Y7ZJMbAG+A300HLKuHkw2SuxKAe8lQ3Gdjv6llLP4\/FMQfs1MMBQ2SUOyaoPihur0X4BDGEDHVGLYbwhawCG6QG6W6yu43jo+uwSt3IyWLlA2vt09+F7CwlDodAaaZ+Q\/j3inZn0GENil\/gd0Ke1rwEjgUBgDqzYIe4lBiBbkkqlZov7r2hjXtLldom\/0LVa+x7QK9+NGjHduUvkfcAdbUwb0G+XuNtg8WV6QeLXbCDeKzo9QJM2phPosEt8ExijG9nG5NuKJK6VTTeXrofudbvEh2TS49r\/eakQMeQXNV2v9Nfbzc9c8U\/GWv2LKkmYEat+r9c7Q+L9nZvOdhLBhOF86c+MWJ5JkLVIX6jUdM1Jbsgk3NU1RRBPHCzQ2crxZR0U3CBC8saKGO3NCNNGx85kpk4QPLfaXJr8jGPkKs4FfqesJv9kcP1OvFeXa70Z8Ygmz1KWTCYX4fuDqlAQhnUVIWafQfesXjjgF5zlCrGWSBQuyOm1nuEARqGTqRgxKxejrm1iYL+URKPFEktJrGSDpRJ3aYlGYRhFwjK8T+lyle+1yrUsYjMM5ZG9l0tBzk3iYvGfWE6XVWYniyvEci1RgjEpZwP6Itwi7jfppHy31Mz3HSeGNZvw8aOAUs5h0olUqqqPFVIhjlcI34wXsMMuWWhEb76ibwPw1U1ibmCzwm+tdil3HEyfhcrdBhetbrW6pHe5RNxiVeTXukTcVEyevupCjA9YEqfT6XlyGXOKuINnQ1GnEy9y6uJeBuFn3qXxxyzmnP8Amz0wy7TM4YgAAAAASUVORK5CYII=","width":19}
%---
%[text:image:2360]
%   data: {"align":"baseline","height":26,"src":"data:image\/png;base64,iVBORw0KGgoAAAANSUhEUgAAAB4AAAApCAYAAAAiT5m3AAADOElEQVR42r1XSWgUQRSNGndR3Ii7V\/EmCCp4DaIe1Oh4EMHBwMjsq47Eg+OOJ0VIICAeFAIeXEADHgyoqChE1ASNhISAxH0jg5gxq+\/FX1A009PTM93T8OiuX7\/q1f+\/+tevqqoCTyKRWBmNRo8Cj4GPwF++I5HII7zT8Xh8eZWTTzgcnh6Lxc5j8kFgHMgCT4G7wBPgl8j\/YBGnfT7f1LJJg8HgQrGQE3diAXVciK5DIsi3o7+devh+4Pf755dMmslkpmGih0LaiHa1hf5k6J0T\/TYrfdMHKz8pFjTr1kFWD9ltvJ8Bt\/Dt1d0rYeG4Y7ZJMbAG+A300HLKuHkw2SuxKAe8lQ3Gdjv6llLP4\/FMQfs1MMBQ2SUOyaoPihur0X4BDGEDHVGLYbwhawCG6QG6W6yu43jo+uwSt3IyWLlA2vt09+F7CwlDodAaaZ+Q\/j3inZn0GENil\/gd0Ke1rwEjgUBgDqzYIe4lBiBbkkqlZov7r2hjXtLldom\/0LVa+x7QK9+NGjHduUvkfcAdbUwb0G+XuNtg8WV6QeLXbCDeKzo9QJM2phPosEt8ExijG9nG5NuKJK6VTTeXrofudbvEh2TS49r\/eakQMeQXNV2v9Nfbzc9c8U\/GWv2LKkmYEat+r9c7Q+L9nZvOdhLBhOF86c+MWJ5JkLVIX6jUdM1Jbsgk3NU1RRBPHCzQ2crxZR0U3CBC8saKGO3NCNNGx85kpk4QPLfaXJr8jGPkKs4FfqesJv9kcP1OvFeXa70Z8Ygmz1KWTCYX4fuDqlAQhnUVIWafQfesXjjgF5zlCrGWSBQuyOm1nuEARqGTqRgxKxejrm1iYL+URKPFEktJrGSDpRJ3aYlGYRhFwjK8T+lyle+1yrUsYjMM5ZG9l0tBzk3iYvGfWE6XVWYniyvEci1RgjEpZwP6Itwi7jfppHy31Mz3HSeGNZvw8aOAUs5h0olUqqqPFVIhjlcI34wXsMMuWWhEb76ibwPw1U1ibmCzwm+tdil3HEyfhcrdBhetbrW6pHe5RNxiVeTXukTcVEyevupCjA9YEqfT6XlyGXOKuINnQ1GnEy9y6uJeBuFn3qXxxyzmnP8Amz0wy7TM4YgAAAAASUVORK5CYII=","width":19}
%---
%[text:image:91f0]
%   data: {"align":"baseline","height":26,"src":"data:image\/png;base64,iVBORw0KGgoAAAANSUhEUgAAAB4AAAApCAYAAAAiT5m3AAADOElEQVR42r1XSWgUQRSNGndR3Ii7V\/EmCCp4DaIe1Oh4EMHBwMjsq47Eg+OOJ0VIICAeFAIeXEADHgyoqChE1ASNhISAxH0jg5gxq+\/FX1A009PTM93T8OiuX7\/q1f+\/+tevqqoCTyKRWBmNRo8Cj4GPwF++I5HII7zT8Xh8eZWTTzgcnh6Lxc5j8kFgHMgCT4G7wBPgl8j\/YBGnfT7f1LJJg8HgQrGQE3diAXVciK5DIsi3o7+devh+4Pf755dMmslkpmGih0LaiHa1hf5k6J0T\/TYrfdMHKz8pFjTr1kFWD9ltvJ8Bt\/Dt1d0rYeG4Y7ZJMbAG+A300HLKuHkw2SuxKAe8lQ3Gdjv6llLP4\/FMQfs1MMBQ2SUOyaoPihur0X4BDGEDHVGLYbwhawCG6QG6W6yu43jo+uwSt3IyWLlA2vt09+F7CwlDodAaaZ+Q\/j3inZn0GENil\/gd0Ke1rwEjgUBgDqzYIe4lBiBbkkqlZov7r2hjXtLldom\/0LVa+x7QK9+NGjHduUvkfcAdbUwb0G+XuNtg8WV6QeLXbCDeKzo9QJM2phPosEt8ExijG9nG5NuKJK6VTTeXrofudbvEh2TS49r\/eakQMeQXNV2v9Nfbzc9c8U\/GWv2LKkmYEat+r9c7Q+L9nZvOdhLBhOF86c+MWJ5JkLVIX6jUdM1Jbsgk3NU1RRBPHCzQ2crxZR0U3CBC8saKGO3NCNNGx85kpk4QPLfaXJr8jGPkKs4FfqesJv9kcP1OvFeXa70Z8Ygmz1KWTCYX4fuDqlAQhnUVIWafQfesXjjgF5zlCrGWSBQuyOm1nuEARqGTqRgxKxejrm1iYL+URKPFEktJrGSDpRJ3aYlGYRhFwjK8T+lyle+1yrUsYjMM5ZG9l0tBzk3iYvGfWE6XVWYniyvEci1RgjEpZwP6Itwi7jfppHy31Mz3HSeGNZvw8aOAUs5h0olUqqqPFVIhjlcI34wXsMMuWWhEb76ibwPw1U1ibmCzwm+tdil3HEyfhcrdBhetbrW6pHe5RNxiVeTXukTcVEyevupCjA9YEqfT6XlyGXOKuINnQ1GnEy9y6uJeBuFn3qXxxyzmnP8Amz0wy7TM4YgAAAAASUVORK5CYII=","width":19}
%---
%[text:image:1408]
%   data: {"align":"baseline","height":24,"src":"data:image\/png;base64,iVBORw0KGgoAAAANSUhEUgAAACgAAAAsCAYAAAAXb\/p7AAADlklEQVR42u1YS0hUYRQejR4UvR+SvaCyWiS2DIoCKWjVYzG2SsYWQ+a8HRwXFkOgqdGLSoSeqyBaJUGLCKGkNpGtgiyScnpYBpaU5au+E+eH03Wu93rvP2OEFw46955z\/u8\/7\/\/3eDQ+4XD4Eui65199AK4d1DEJcBLgfw0wGo2ujkQi1yoqKlY5AFgDqnWw5lpaMxaLrbBkDgaDm7HIEKjTCcjxPqFQqABrdYEGsHaRXaESCAyCXpNFMwwuReDw\/77xusurQGJna0ysPR2K94CnGdQKesnUyu92J5PJaelk4c51Chzcu9dp0O9nkO3yvd\/vnwqlEbzvAf2yoB5sIkwyBt0doJ+0Qbdu2AXyC8UrOVuNQH5Q3IJegfrTfH8ikwC\/A9hjsdaYCQQCG6D4vVj0DShB7+HKXMGaQ\/EFADEGrPhJdn1GAhoZvRDKX4jFTvt8vhlWchSD4D0OGiE5AH9eXl4+PxNF+IoCB8sccyBfreQB8qJWcAC0CYqHeYE2cqEDNTmQvcc6hlC+CnXWq\/PCtTtdeGGr0HNWFz7aeTcr\/eT1eqe41PVWJYxd62zkovwXqQ6Cv8vErm9r8MZNEctLVf9Ph4GwecSOjPSQu8V28e6chmRrFPq28btHJhhSBKDIBH0B726HEDilAWCtsGCx6MmjMNgaHqBkiygPl3UmHE1OrjOEiqoqsqDHGizYxrqGAXCOriL9jJUOoq3lO9UDQItpQGD3PtXZRapFHB51oSchwiWsDWA8Hl8Cpb2svA+WWO4AXJ7Q0Z1IJObq7sUBN\/UQLr0hrHfANSCaQihrobiBf+dC+QOxSNk4MrdEyN1VvRy6T+DdVbOp2wpcCyutN8yDahjtVZ3A4sS2ALwflYxhaK1n0C22QfL8dkuNRYZBlKxxRLj6gg3XNgj+w8beLOriHcv50gocJ8wsfP\/ASgfMDlRs8XzwfGPeThMr2QcpEqFprJmPSoSwSvMYfCcFX6nFlNPEfAGrUuC1GkhxQpspTnTfKysrFxl5qqqqZouy0mU81ZmApB6cp2vKlrFVk+Z7TGRuPOsXRZSNfGYmEO\/ofoXmOiKKS3Ga6wPYeRNxk9Vo49BONAIL1mUdIBa9bxPgn+zMOkAqIXDdQbp1sKAy6uWeyWd0Cyuk5LDh3pTtez\/NSVJqNwbh5kNZB8htkTrP1zHA9dM1nZ07nExmc50ZQGA7M+GxyLden9MA\/OLm\/KKe34xPoKAD7jPkAAAAAElFTkSuQmCC","width":22}
%---
