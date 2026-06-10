%[text] %[text:anchor:H_733CD629] # Convergence study on transverse shear-locking
%[text:tableOfContents]{"heading":"Table of Contents"}
%[text] Another way to identify transverse shear locking is to perform a convergence study for different slenderness levels of the structure. A reduced convergence rate as the slenderness level increases indicates the presence of transverse shear-locking. In extreme transverse shear-locking cases the convergence rate drops to nearly zero
%[text] 
%[text] Go to [THIS](file:.\main_Chapter4_Locking_ReissnerMindlinCantileverPlate.m) Section in the main driver script *Transverse Shear Locking in Shear-Deformable Plates* *`-`* *Cantilever Plate Subject to Vertical Distributed Load*
%[text] Go to [PREVIOUS](file:.\main_Chapter42_ReissnerMindlinCantileverPlate_SlendernessStudy.m) Section *Slenderness study on transverse shear-locking*
%[text] 
%[text] **Author: Andreas Apostolatos, PhD (aapostol@mathworks.com)**
%[text] **Date: 10.06.2026**
%%
%[text] %[text:anchor:H_F203EBCC] ## Preamble
vars = ["xyPostProc" "wEx" "propStr" "X0" "Y0" "XLx" "YLy" ...
    "numNodesElQ1" "numNodesElQ2" "numEly"];
str = "This m-file should be ran from parent m-file " + ...
    "main_Chapter4_Locking_ReissnerMindlinCantileverPlate.m";
for ii = 1:length(vars)
    if ~exist(vars(ii), 'var') 
        error(strcat("Variable ", vars(ii), " is not defined. ", str));
    end
end
%%
%[text] %[text:anchor:H_E9F2F5E6] ## Analytical solution for the vertical deflection field
syms wExCon(X);
wExCon(X) = subs(wEx,{Lsym, psym, Gsym, Esym, Isym, Aqsym}, ...
    {XLx, propStr.pBar, propStr.G, propStr.E, YLy*(propStr.t^3)/12, ...
    propStr.alpha*YLy*propStr.t});
%[text] %[text:anchor:H_17D0DB6A] ## Initialize auxiliary arrays
tipDeflection = zeros(numRef, 2);
tipDeflectionAnalytical = zeros(numRef, 1);
numElmnts = zeros(numRef, 1);
%%
%[text] %[text:anchor:H_707313A0] ## Loop over all the refinement levels
for ii = 1:numRef
%[text] %[text:anchor:H_BE9F5861] ### **Quadrilateral mesh**
    numElxRef = ii;
    numElyRef = numEly;
    numElmnts(ii, 1) = numElxRef*numElyRef;
    [nodesXRef, nodesYRef] = ...
        generateQuadrilateralMesh(X0, XLx, Y0, YLy, numElxRef, numElyRef);
%[text] %[text:anchor:H_D8741169] ### **Bilinear quadrilateral mesh**
    mshQ1Ref = generateBilinearQuadrilateralMesh ...
        (numElxRef, numElyRef, numElmnts(ii, 1), ...
        nodesXRef, nodesYRef, numNodesElQ1);
%[text] %[text:anchor:H_24D456AB] ### **Biquadratic quadrilateral mesh**
    [mshQ2Ref, ~] = generateBiquadraticQuadrilateralMesh ...
        (numElxRef, numElyRef, nodesXRef, nodesYRef, mshQ1Ref, ...
        numNodesElQ1, numNodesElQ2);
%[text] %[text:anchor:H_B00CC111] ### Element in mesh containing the evaluation point
    elPstProc = nan;
    for jj = 1:numElmnts(ii, 1)
        [in, ~] = inpolygon(xyPostProc(1), xyPostProc(2), ...
            mshQ1Ref.nodes(mshQ1Ref.elements(jj, :), 1), ...
            mshQ1Ref.nodes(mshQ1Ref.elements(jj, :), 2));
        if in
            elPstProc = jj;
            break;
        end
    end
    if isnan(elPstProc)
        error("No element found containing the point where the " + ...
            "postprocessing results are evaluated");
    end
    propNewtonRaphson.eps = 1e-12;
    propNewtonRaphson.maxIt = 20;
    [xiEta, isConverged] = ...
        computePointCoordinatesOnCanonicalBilinearQuadrilateral ...
        (xyPostProc', mshQ1Ref.nodes(mshQ1Ref.elements(elPstProc, :), 1:2)', ...
        propNewtonRaphson);
    if ~isConverged
        error("The parametric coordinates of point (%d, %d) can be found " + ...
            "in no element %d", xyPostProc(1), xyPostProc(2), elPstProc);
    end
%[text] %[text:anchor:H_CC39AB40] ### Element Freedom Tables (EFTs) for both meshes at evaluation point
    EFTQ1Postproc = [3*mshQ1Ref.elements(elPstProc, :) - 2
                     3*mshQ1Ref.elements(elPstProc, :) - 1
                     3*mshQ1Ref.elements(elPstProc, :)];
    EFTQ1Postproc = EFTQ1Postproc(:);
    EFTQ2Postproc = [3*mshQ2Ref.elements(elPstProc, :) - 2
                     3*mshQ2Ref.elements(elPstProc, :)-  1
                     3*mshQ2Ref.elements(elPstProc, :)];
    EFTQ2Postproc = EFTQ2Postproc(:);
%[text] %[text:anchor:H_5979C943] ### Basis functions matrices at evaluation point
%[text] ... for the bilinear mesh
    ShapeFunctionValuesQ1 = ...
        computeBilinearBasisFunctionsAndFirstDerivatives ...
        (xiEta(1, 1), xiEta(2, 1));
    NmtxQ1 = zeros(3, 3*numNodesElQ1);
    for jj = 1:numel(ShapeFunctionValuesQ1(:, 1))
        NmtxQ1(1, 3*jj-2) = ShapeFunctionValuesQ1(jj, 1);
        NmtxQ1(2, 3*jj-1) = ShapeFunctionValuesQ1(jj, 1);
        NmtxQ1(3, 3*jj) = ShapeFunctionValuesQ1(jj, 1);
    end
%[text] %[text:anchor:H_2DE0F805] ... for the bilinear mesh
    ShapeFunctionValuesQ2 = ...
        computeBiquadraticBasisFunctionsAndFirstDerivatives ...
        (xiEta(1, 1), xiEta(2, 1));
    NmtxQ2 = zeros(3, 3*numNodesElQ2);
    for jj = 1:numel(ShapeFunctionValuesQ2(:, 1))
        NmtxQ2(1, 3*jj-2) = ShapeFunctionValuesQ2(jj, 1);
        NmtxQ2(2, 3*jj-1) = ShapeFunctionValuesQ2(jj, 1);
        NmtxQ2(3, 3*jj) = ShapeFunctionValuesQ2(jj, 1);
    end
%[text] %[text:anchor:H_8529AAF4] ### Boundary conditions
%[text] ... for the bilinear quadrilateral mesh
    numNodesQ1Ref = height(mshQ1Ref.nodes);    
    numDOFsQ1Ref = 3*numNodesQ1Ref;
    freeDOFsQ1Ref = ...
        utilityFunctionsCantileverPlate.getFreeDOFsBilinearQuadrilateralMesh ...
        (numNodesQ1Ref, numEly);
%[text] %[text:anchor:H_789012AE] ... for the bilinear quadrilateral mesh
    numNodesQ2Ref = height(mshQ2Ref.nodes);
    numDOFsQ2Ref = 3*numNodesQ2Ref;
    freeDOFsQ2Ref = ...
        utilityFunctionsCantileverPlate.getFreeDOFsBiquadraticQuadrilateralMesh ...
        (numNodesQ1Ref, numNodesQ2Ref, numEly);
%[text] %[text:anchor:H_2B81D4BD] ### Master stiffness matrices
%[text] ... corresponding to the bilinear quadrilateral mesh
    computeStiffMatrixandForceVectorQ1 = 'undefined';
    computeBasisFunctionsAndDerivsQ1 = ... 
        @computeBilinearBasisFunctionsAndFirstDerivatives;
    [KQ1Ref, FQ1Ref] = ...
        computeMasterStiffMatrixandForceVctReissnerMindlinPlatePageWise ...
        (mshQ1Ref, computeBasisFunctionsAndDerivsQ1, ...
        computeStiffMatrixandForceVectorQ1, propStr);
%[text] ... corresponding to the biquadratic quadrilateral mesh
    computeStiffMatrixandForceVectorQ2 = 'undefined';
    computeBasisFunctionsAndDerivsQ2 = ...
        @computeBiquadraticBasisFunctionsAndFirstDerivatives;
    [KQ2Ref, FQ2Ref] = ...
        computeMasterStiffMatrixandForceVctReissnerMindlinPlatePageWise ...
        (mshQ2Ref, computeBasisFunctionsAndDerivsQ2, ...
        computeStiffMatrixandForceVectorQ2, propStr);
%[text] ![](text:image:5f63) **Try**
%[text] Compute the master stiffness matrix for the bilinear elements using reduced integration
%[text] ```matlabCodeExample
%[text] computeStiffMatrixandForceVectorRI = ...
%[text]     @computeElementStiffMatrixandForceVectorReissnerMindlinPlateRI;
%[text] [KQ1RI, FQ1RI] = ...
%[text]     computeMasterStiffMatrixandForceVectorReissnerMindlinPlate ...
%[text]     (msh, computeStiffMatrixandForceVectorRI, propStr);
%[text] ```
%[text] %[text:anchor:H_5C523FE4] ### Solution of the linear equation systems
%[text] ... corresponding to the bilinear quadrilateral mesh
    uQ1Ref = zeros(numDOFsQ1Ref, 1);
    uQ1Ref(freeDOFsQ1Ref, 1) = ...
        KQ1Ref(freeDOFsQ1Ref, freeDOFsQ1Ref)\FQ1Ref(freeDOFsQ1Ref);
%[text] ... corresponding to the biquadratic quadrilateral mesh
    uQ2Ref = zeros(numDOFsQ2Ref, 1);
    uQ2Ref(freeDOFsQ2Ref, 1) = ...
        KQ2Ref(freeDOFsQ2Ref, freeDOFsQ2Ref)\FQ2Ref(freeDOFsQ2Ref);
%[text] ![](text:image:0004) **Try**
%[text] Solve the linear equation system corresponding to the bilinear mesh using selective reduced-integration
%[text] ```matlabCodeExample
%[text] try
%[text]     uQ1RIRef = KQ1RI(freeDOFsQ1Ref, freeDOFsQ1Ref)\FQ1RI(freeDOFsQ1Ref);
%[text] catch ME
%[text]     rethrow(ME)
%[text] end
%[text] ```
%[text] %[text:anchor:H_323B4196] ### Vertical deflections at the evaluation point
%[text] ... corresponding to the bilinear quadrilateral mesh
    wPostprocQ1 = NmtxQ1*uQ1Ref(EFTQ1Postproc);
    tipDeflection(ii, 1) = wPostprocQ1(1);
%[text] %[text:anchor:H_E352EDAD] ... corresponding to the bilinear quadrilateral mesh
    wPostprocQ2 = NmtxQ2*uQ2Ref(EFTQ2Postproc);
    tipDeflection(ii, 2) = wPostprocQ2(1);
%[text] ![](text:image:392b) **Try**
%[text] Compute the vertical deflection at the selected point using the bilinear quadrilateral elements and reduced integration
%[text] ```matlabCodeExample
%[text] try
%[text]     wPostprocQ1Red = NmtxQ1*uQ1RIRef(EFTQ1Postproc);
%[text] catch ME
%[text]     rethrow(ME)
%[text] end
%[text] ```
%[text] ![](text:image:0be1) **Reflect**
%[text] Can the solution be computed for this kind of elements and if not, why?
%[text] ... corresponding to the analytical solution
tipDeflectionAnalytical(ii, 1) = double(wExCon(XLx));
%[text] 
end
%%
%[text] %[text:anchor:H_F41A9C84] ## Convergence graph regarding the refinement study
plot(numElmnts, abs(tipDeflection(:, 1)), "-bo", LineWidth=2) %[text:anchor:M_B8524D00]
hold on;
plot(numElmnts, abs(tipDeflection(:, 2)), "-ro", LineWidth=2)
plot(numElmnts, abs(tipDeflectionAnalytical), "-k", LineWidth=2)
hold off;
grid on;
legend("Q1", "Q2", "analytical", Location='best')
xlabel("Number of elements")
ylabel("Tip displacement")
title("Convergence study")
%[text] ![](text:image:0f79) **Try**
%[text] Change the value of the thickness by means of the provided slider in code section **Problem parameters** of [this](file:.\main_Chapter4_Locking_ReissnerMindlinCantileverPlate.m) script. Change its value to the extrems and the script is going to be automatically executed.
%[text] ![](text:image:40e4) **Reflect**
%[text] When changing the value of the thickness to the minimum possible through the slider:
%[text] - What do you observe regarding the distribution of the transverse-shear force field in the corresponding figure?
%[text] - Which solution is off judging from the [convergence graph](internal:M_B8524D00) and why? \
%[text] When changing the value of the thickness to the maximum possible through the slider:
%[text] - Which solution is off in regard to the distribution of the transverse-shear force field in the corresponding figure?
%[text] - What can you tell about the convergence rates in the [convergence graph](internal:M_B8524D00) of the different elements? Is this to be expected? \

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"inline","rightPanelPercent":40}
%---
%[text:image:5f63]
%   data: {"align":"baseline","height":26,"src":"data:image\/png;base64,iVBORw0KGgoAAAANSUhEUgAAAB4AAAApCAYAAAAiT5m3AAADOElEQVR42r1XSWgUQRSNGndR3Ii7V\/EmCCp4DaIe1Oh4EMHBwMjsq47Eg+OOJ0VIICAeFAIeXEADHgyoqChE1ASNhISAxH0jg5gxq+\/FX1A009PTM93T8OiuX7\/q1f+\/+tevqqoCTyKRWBmNRo8Cj4GPwF++I5HII7zT8Xh8eZWTTzgcnh6Lxc5j8kFgHMgCT4G7wBPgl8j\/YBGnfT7f1LJJg8HgQrGQE3diAXVciK5DIsi3o7+devh+4Pf755dMmslkpmGih0LaiHa1hf5k6J0T\/TYrfdMHKz8pFjTr1kFWD9ltvJ8Bt\/Dt1d0rYeG4Y7ZJMbAG+A300HLKuHkw2SuxKAe8lQ3Gdjv6llLP4\/FMQfs1MMBQ2SUOyaoPihur0X4BDGEDHVGLYbwhawCG6QG6W6yu43jo+uwSt3IyWLlA2vt09+F7CwlDodAaaZ+Q\/j3inZn0GENil\/gd0Ke1rwEjgUBgDqzYIe4lBiBbkkqlZov7r2hjXtLldom\/0LVa+x7QK9+NGjHduUvkfcAdbUwb0G+XuNtg8WV6QeLXbCDeKzo9QJM2phPosEt8ExijG9nG5NuKJK6VTTeXrofudbvEh2TS49r\/eakQMeQXNV2v9Nfbzc9c8U\/GWv2LKkmYEat+r9c7Q+L9nZvOdhLBhOF86c+MWJ5JkLVIX6jUdM1Jbsgk3NU1RRBPHCzQ2crxZR0U3CBC8saKGO3NCNNGx85kpk4QPLfaXJr8jGPkKs4FfqesJv9kcP1OvFeXa70Z8Ygmz1KWTCYX4fuDqlAQhnUVIWafQfesXjjgF5zlCrGWSBQuyOm1nuEARqGTqRgxKxejrm1iYL+URKPFEktJrGSDpRJ3aYlGYRhFwjK8T+lyle+1yrUsYjMM5ZG9l0tBzk3iYvGfWE6XVWYniyvEci1RgjEpZwP6Itwi7jfppHy31Mz3HSeGNZvw8aOAUs5h0olUqqqPFVIhjlcI34wXsMMuWWhEb76ibwPw1U1ibmCzwm+tdil3HEyfhcrdBhetbrW6pHe5RNxiVeTXukTcVEyevupCjA9YEqfT6XlyGXOKuINnQ1GnEy9y6uJeBuFn3qXxxyzmnP8Amz0wy7TM4YgAAAAASUVORK5CYII=","width":19}
%---
%[text:image:0004]
%   data: {"align":"baseline","height":26,"src":"data:image\/png;base64,iVBORw0KGgoAAAANSUhEUgAAAB4AAAApCAYAAAAiT5m3AAADOElEQVR42r1XSWgUQRSNGndR3Ii7V\/EmCCp4DaIe1Oh4EMHBwMjsq47Eg+OOJ0VIICAeFAIeXEADHgyoqChE1ASNhISAxH0jg5gxq+\/FX1A009PTM93T8OiuX7\/q1f+\/+tevqqoCTyKRWBmNRo8Cj4GPwF++I5HII7zT8Xh8eZWTTzgcnh6Lxc5j8kFgHMgCT4G7wBPgl8j\/YBGnfT7f1LJJg8HgQrGQE3diAXVciK5DIsi3o7+devh+4Pf755dMmslkpmGih0LaiHa1hf5k6J0T\/TYrfdMHKz8pFjTr1kFWD9ltvJ8Bt\/Dt1d0rYeG4Y7ZJMbAG+A300HLKuHkw2SuxKAe8lQ3Gdjv6llLP4\/FMQfs1MMBQ2SUOyaoPihur0X4BDGEDHVGLYbwhawCG6QG6W6yu43jo+uwSt3IyWLlA2vt09+F7CwlDodAaaZ+Q\/j3inZn0GENil\/gd0Ke1rwEjgUBgDqzYIe4lBiBbkkqlZov7r2hjXtLldom\/0LVa+x7QK9+NGjHduUvkfcAdbUwb0G+XuNtg8WV6QeLXbCDeKzo9QJM2phPosEt8ExijG9nG5NuKJK6VTTeXrofudbvEh2TS49r\/eakQMeQXNV2v9Nfbzc9c8U\/GWv2LKkmYEat+r9c7Q+L9nZvOdhLBhOF86c+MWJ5JkLVIX6jUdM1Jbsgk3NU1RRBPHCzQ2crxZR0U3CBC8saKGO3NCNNGx85kpk4QPLfaXJr8jGPkKs4FfqesJv9kcP1OvFeXa70Z8Ygmz1KWTCYX4fuDqlAQhnUVIWafQfesXjjgF5zlCrGWSBQuyOm1nuEARqGTqRgxKxejrm1iYL+URKPFEktJrGSDpRJ3aYlGYRhFwjK8T+lyle+1yrUsYjMM5ZG9l0tBzk3iYvGfWE6XVWYniyvEci1RgjEpZwP6Itwi7jfppHy31Mz3HSeGNZvw8aOAUs5h0olUqqqPFVIhjlcI34wXsMMuWWhEb76ibwPw1U1ibmCzwm+tdil3HEyfhcrdBhetbrW6pHe5RNxiVeTXukTcVEyevupCjA9YEqfT6XlyGXOKuINnQ1GnEy9y6uJeBuFn3qXxxyzmnP8Amz0wy7TM4YgAAAAASUVORK5CYII=","width":19}
%---
%[text:image:392b]
%   data: {"align":"baseline","height":26,"src":"data:image\/png;base64,iVBORw0KGgoAAAANSUhEUgAAAB4AAAApCAYAAAAiT5m3AAADOElEQVR42r1XSWgUQRSNGndR3Ii7V\/EmCCp4DaIe1Oh4EMHBwMjsq47Eg+OOJ0VIICAeFAIeXEADHgyoqChE1ASNhISAxH0jg5gxq+\/FX1A009PTM93T8OiuX7\/q1f+\/+tevqqoCTyKRWBmNRo8Cj4GPwF++I5HII7zT8Xh8eZWTTzgcnh6Lxc5j8kFgHMgCT4G7wBPgl8j\/YBGnfT7f1LJJg8HgQrGQE3diAXVciK5DIsi3o7+devh+4Pf755dMmslkpmGih0LaiHa1hf5k6J0T\/TYrfdMHKz8pFjTr1kFWD9ltvJ8Bt\/Dt1d0rYeG4Y7ZJMbAG+A300HLKuHkw2SuxKAe8lQ3Gdjv6llLP4\/FMQfs1MMBQ2SUOyaoPihur0X4BDGEDHVGLYbwhawCG6QG6W6yu43jo+uwSt3IyWLlA2vt09+F7CwlDodAaaZ+Q\/j3inZn0GENil\/gd0Ke1rwEjgUBgDqzYIe4lBiBbkkqlZov7r2hjXtLldom\/0LVa+x7QK9+NGjHduUvkfcAdbUwb0G+XuNtg8WV6QeLXbCDeKzo9QJM2phPosEt8ExijG9nG5NuKJK6VTTeXrofudbvEh2TS49r\/eakQMeQXNV2v9Nfbzc9c8U\/GWv2LKkmYEat+r9c7Q+L9nZvOdhLBhOF86c+MWJ5JkLVIX6jUdM1Jbsgk3NU1RRBPHCzQ2crxZR0U3CBC8saKGO3NCNNGx85kpk4QPLfaXJr8jGPkKs4FfqesJv9kcP1OvFeXa70Z8Ygmz1KWTCYX4fuDqlAQhnUVIWafQfesXjjgF5zlCrGWSBQuyOm1nuEARqGTqRgxKxejrm1iYL+URKPFEktJrGSDpRJ3aYlGYRhFwjK8T+lyle+1yrUsYjMM5ZG9l0tBzk3iYvGfWE6XVWYniyvEci1RgjEpZwP6Itwi7jfppHy31Mz3HSeGNZvw8aOAUs5h0olUqqqPFVIhjlcI34wXsMMuWWhEb76ibwPw1U1ibmCzwm+tdil3HEyfhcrdBhetbrW6pHe5RNxiVeTXukTcVEyevupCjA9YEqfT6XlyGXOKuINnQ1GnEy9y6uJeBuFn3qXxxyzmnP8Amz0wy7TM4YgAAAAASUVORK5CYII=","width":19}
%---
%[text:image:0be1]
%   data: {"align":"baseline","height":24,"src":"data:image\/png;base64,iVBORw0KGgoAAAANSUhEUgAAACgAAAAsCAYAAAAXb\/p7AAADlklEQVR42u1YS0hUYRQejR4UvR+SvaCyWiS2DIoCKWjVYzG2SsYWQ+a8HRwXFkOgqdGLSoSeqyBaJUGLCKGkNpGtgiyScnpYBpaU5au+E+eH03Wu93rvP2OEFw46955z\/u8\/7\/\/3eDQ+4XD4Eui65199AK4d1DEJcBLgfw0wGo2ujkQi1yoqKlY5AFgDqnWw5lpaMxaLrbBkDgaDm7HIEKjTCcjxPqFQqABrdYEGsHaRXaESCAyCXpNFMwwuReDw\/77xusurQGJna0ysPR2K94CnGdQKesnUyu92J5PJaelk4c51Chzcu9dp0O9nkO3yvd\/vnwqlEbzvAf2yoB5sIkwyBt0doJ+0Qbdu2AXyC8UrOVuNQH5Q3IJegfrTfH8ikwC\/A9hjsdaYCQQCG6D4vVj0DShB7+HKXMGaQ\/EFADEGrPhJdn1GAhoZvRDKX4jFTvt8vhlWchSD4D0OGiE5AH9eXl4+PxNF+IoCB8sccyBfreQB8qJWcAC0CYqHeYE2cqEDNTmQvcc6hlC+CnXWq\/PCtTtdeGGr0HNWFz7aeTcr\/eT1eqe41PVWJYxd62zkovwXqQ6Cv8vErm9r8MZNEctLVf9Ph4GwecSOjPSQu8V28e6chmRrFPq28btHJhhSBKDIBH0B726HEDilAWCtsGCx6MmjMNgaHqBkiygPl3UmHE1OrjOEiqoqsqDHGizYxrqGAXCOriL9jJUOoq3lO9UDQItpQGD3PtXZRapFHB51oSchwiWsDWA8Hl8Cpb2svA+WWO4AXJ7Q0Z1IJObq7sUBN\/UQLr0hrHfANSCaQihrobiBf+dC+QOxSNk4MrdEyN1VvRy6T+DdVbOp2wpcCyutN8yDahjtVZ3A4sS2ALwflYxhaK1n0C22QfL8dkuNRYZBlKxxRLj6gg3XNgj+w8beLOriHcv50gocJ8wsfP\/ASgfMDlRs8XzwfGPeThMr2QcpEqFprJmPSoSwSvMYfCcFX6nFlNPEfAGrUuC1GkhxQpspTnTfKysrFxl5qqqqZouy0mU81ZmApB6cp2vKlrFVk+Z7TGRuPOsXRZSNfGYmEO\/ofoXmOiKKS3Ga6wPYeRNxk9Vo49BONAIL1mUdIBa9bxPgn+zMOkAqIXDdQbp1sKAy6uWeyWd0Cyuk5LDh3pTtez\/NSVJqNwbh5kNZB8htkTrP1zHA9dM1nZ07nExmc50ZQGA7M+GxyLden9MA\/OLm\/KKe34xPoKAD7jPkAAAAAElFTkSuQmCC","width":22}
%---
%[text:image:0f79]
%   data: {"align":"baseline","height":26,"src":"data:image\/png;base64,iVBORw0KGgoAAAANSUhEUgAAAB4AAAApCAYAAAAiT5m3AAADOElEQVR42r1XSWgUQRSNGndR3Ii7V\/EmCCp4DaIe1Oh4EMHBwMjsq47Eg+OOJ0VIICAeFAIeXEADHgyoqChE1ASNhISAxH0jg5gxq+\/FX1A009PTM93T8OiuX7\/q1f+\/+tevqqoCTyKRWBmNRo8Cj4GPwF++I5HII7zT8Xh8eZWTTzgcnh6Lxc5j8kFgHMgCT4G7wBPgl8j\/YBGnfT7f1LJJg8HgQrGQE3diAXVciK5DIsi3o7+devh+4Pf755dMmslkpmGih0LaiHa1hf5k6J0T\/TYrfdMHKz8pFjTr1kFWD9ltvJ8Bt\/Dt1d0rYeG4Y7ZJMbAG+A300HLKuHkw2SuxKAe8lQ3Gdjv6llLP4\/FMQfs1MMBQ2SUOyaoPihur0X4BDGEDHVGLYbwhawCG6QG6W6yu43jo+uwSt3IyWLlA2vt09+F7CwlDodAaaZ+Q\/j3inZn0GENil\/gd0Ke1rwEjgUBgDqzYIe4lBiBbkkqlZov7r2hjXtLldom\/0LVa+x7QK9+NGjHduUvkfcAdbUwb0G+XuNtg8WV6QeLXbCDeKzo9QJM2phPosEt8ExijG9nG5NuKJK6VTTeXrofudbvEh2TS49r\/eakQMeQXNV2v9Nfbzc9c8U\/GWv2LKkmYEat+r9c7Q+L9nZvOdhLBhOF86c+MWJ5JkLVIX6jUdM1Jbsgk3NU1RRBPHCzQ2crxZR0U3CBC8saKGO3NCNNGx85kpk4QPLfaXJr8jGPkKs4FfqesJv9kcP1OvFeXa70Z8Ygmz1KWTCYX4fuDqlAQhnUVIWafQfesXjjgF5zlCrGWSBQuyOm1nuEARqGTqRgxKxejrm1iYL+URKPFEktJrGSDpRJ3aYlGYRhFwjK8T+lyle+1yrUsYjMM5ZG9l0tBzk3iYvGfWE6XVWYniyvEci1RgjEpZwP6Itwi7jfppHy31Mz3HSeGNZvw8aOAUs5h0olUqqqPFVIhjlcI34wXsMMuWWhEb76ibwPw1U1ibmCzwm+tdil3HEyfhcrdBhetbrW6pHe5RNxiVeTXukTcVEyevupCjA9YEqfT6XlyGXOKuINnQ1GnEy9y6uJeBuFn3qXxxyzmnP8Amz0wy7TM4YgAAAAASUVORK5CYII=","width":19}
%---
%[text:image:40e4]
%   data: {"align":"baseline","height":24,"src":"data:image\/png;base64,iVBORw0KGgoAAAANSUhEUgAAACgAAAAsCAYAAAAXb\/p7AAADlklEQVR42u1YS0hUYRQejR4UvR+SvaCyWiS2DIoCKWjVYzG2SsYWQ+a8HRwXFkOgqdGLSoSeqyBaJUGLCKGkNpGtgiyScnpYBpaU5au+E+eH03Wu93rvP2OEFw46955z\/u8\/7\/\/3eDQ+4XD4Eui65199AK4d1DEJcBLgfw0wGo2ujkQi1yoqKlY5AFgDqnWw5lpaMxaLrbBkDgaDm7HIEKjTCcjxPqFQqABrdYEGsHaRXaESCAyCXpNFMwwuReDw\/77xusurQGJna0ysPR2K94CnGdQKesnUyu92J5PJaelk4c51Chzcu9dp0O9nkO3yvd\/vnwqlEbzvAf2yoB5sIkwyBt0doJ+0Qbdu2AXyC8UrOVuNQH5Q3IJegfrTfH8ikwC\/A9hjsdaYCQQCG6D4vVj0DShB7+HKXMGaQ\/EFADEGrPhJdn1GAhoZvRDKX4jFTvt8vhlWchSD4D0OGiE5AH9eXl4+PxNF+IoCB8sccyBfreQB8qJWcAC0CYqHeYE2cqEDNTmQvcc6hlC+CnXWq\/PCtTtdeGGr0HNWFz7aeTcr\/eT1eqe41PVWJYxd62zkovwXqQ6Cv8vErm9r8MZNEctLVf9Ph4GwecSOjPSQu8V28e6chmRrFPq28btHJhhSBKDIBH0B726HEDilAWCtsGCx6MmjMNgaHqBkiygPl3UmHE1OrjOEiqoqsqDHGizYxrqGAXCOriL9jJUOoq3lO9UDQItpQGD3PtXZRapFHB51oSchwiWsDWA8Hl8Cpb2svA+WWO4AXJ7Q0Z1IJObq7sUBN\/UQLr0hrHfANSCaQihrobiBf+dC+QOxSNk4MrdEyN1VvRy6T+DdVbOp2wpcCyutN8yDahjtVZ3A4sS2ALwflYxhaK1n0C22QfL8dkuNRYZBlKxxRLj6gg3XNgj+w8beLOriHcv50gocJ8wsfP\/ASgfMDlRs8XzwfGPeThMr2QcpEqFprJmPSoSwSvMYfCcFX6nFlNPEfAGrUuC1GkhxQpspTnTfKysrFxl5qqqqZouy0mU81ZmApB6cp2vKlrFVk+Z7TGRuPOsXRZSNfGYmEO\/ofoXmOiKKS3Ga6wPYeRNxk9Vo49BONAIL1mUdIBa9bxPgn+zMOkAqIXDdQbp1sKAy6uWeyWd0Cyuk5LDh3pTtez\/NSVJqNwbh5kNZB8htkTrP1zHA9dM1nZ07nExmc50ZQGA7M+GxyLden9MA\/OLm\/KKe34xPoKAD7jPkAAAAAElFTkSuQmCC","width":22}
%---
