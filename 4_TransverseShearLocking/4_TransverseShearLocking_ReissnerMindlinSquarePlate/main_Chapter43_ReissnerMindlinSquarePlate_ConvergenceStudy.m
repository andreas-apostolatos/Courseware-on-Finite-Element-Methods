%[text] %[text:anchor:H_4E8CB83F] # Convergence study on transverse shear-locking
%[text:tableOfContents]{"heading":"Table of Contents"}
%[text] Another way to identify transverse shear locking is to perform a convergence study for different slenderness levels of the structure. A reduced convergence rate as the slenderness level increases indicates the presence of transverse shear-locking. In extreme transverse shear-locking cases the convergence rate drops to nearly zero
%[text] 
%[text] Go to [THIS](file:.\main_Chapter4_Locking_ReissnerMindlinSquarePlate.m) Section in the main driver script *Transverse-Shear Locking in Shear-Deformable Plates* *`-`* *Two Sided Clampled Plate Subject to Vertical Distributed Load*
%[text] Go to [PREVIOUS](file:.\main_Chapter42_ReissnerMindlinSquarePlate_SlendernessStudy.m) Section *Slenderness study on transverse shear-locking*
%[text] 
%[text] **Author: Andreas Apostolatos, PhD (aapostol@mathworks.com)**
%[text] **Date: 10.06.2026**
%%
%[text] %[text:anchor:H_F203EBCC] ## Preamble
vars = ["numRef" "wKirchhoff" "X0" "XLx" "Y0" "YLy" ...
    "numNodesElQ1" "numNodesElQ2" "xyPostProc" "propStr"];
str = "This m-file should be run from parent m-file " + ...
    "main_Chapter4_Locking_ReissnerMindlinSquarePlate.m";
for ii = 1:length(vars)
    if ~exist(vars(ii), 'var') 
        error(strcat("Variable ", vars(ii), " is not defined. ", str));
    end
end
%%
%[text] %[text:anchor:H_4F9F49BD] ## Initialize auxiliary arrays
tipDeflection = zeros(numRef, 3);
numElmnts = zeros(numRef, 1);
tipDeflectionKirchhoff = repmat(wKirchhoff, numRef, 1);
%%
%[text] %[text:anchor:H_17C21D6D] ## Loop over all the refinement levels
for ii = 1:numRef
%[text] %[text:anchor:H_F5D504C2] ### **Generation of a quadrilateral mesh**
    numElxRef = ii;
    numElyRef = ii;
    numElRef = numElxRef*numElyRef;
    numElmnts(ii, 1) = numElRef;
    [nodesXRef, nodesYRef] = ...
        generateQuadrilateralMesh(X0, XLx, Y0, YLy, numElxRef, numElyRef);
%[text] %[text:anchor:H_EEFA9339] ### **Bilinear quadrilateral mesh**
    mshQ1Ref = generateBilinearQuadrilateralMesh ...
        (numElxRef, numElyRef, numElRef, ...
        nodesXRef, nodesYRef, numNodesElQ1);
%[text] %[text:anchor:H_B8168A1A] ### **Biquadratic quadrilateral mesh**
    [mshQ2Ref, biasXRef] = generateBiquadraticQuadrilateralMesh ...
        (numElxRef, numElyRef, nodesXRef, nodesYRef, ...
        mshQ1Ref, numNodesElQ1, numNodesElQ2);
%[text] %[text:anchor:H_8FAAB408] ### Element containing the evaluation point
    elPstProc = nan;
    for jj = 1:numElRef
        [in, ~] = inpolygon(xyPostProc(1), xyPostProc(2), ...
            mshQ1Ref.nodes(mshQ1Ref.elements(jj, :), 1), ...
            mshQ1Ref.nodes(mshQ1Ref.elements(jj, :), 2));
        if in
            elPstProc = jj;
            break;
        end
    end
    if isnan(elPstProc)
        error("No element found containing the point where " + ...
            "the postprocessing results are evaluated");
    end
    propNewtonRaphson.eps = 1e-12;
    propNewtonRaphson.maxIt = 20;
    [xiEta, isConverged] = ...
        computePointCoordinatesOnCanonicalBilinearQuadrilateral ...
        (xyPostProc', mshQ1Ref.nodes(mshQ1Ref.elements(elPstProc, :), 1:2)', ...
        propNewtonRaphson);
    if ~isConverged
        error("The parametric coordinates of point (%d, %d) can not be " + ...
            "found in element %d", xyPostProc(1), xyPostProc(2), elPstProc);
    end
%[text] %[text:anchor:H_80B3914D] ### Element Freedom Tables (EFTs) for both meshes at the evaluation points
    EFTQ1Postproc = [3*mshQ1Ref.elements(elPstProc, :) - 2
                     3*mshQ1Ref.elements(elPstProc, :) - 1
                     3*mshQ1Ref.elements(elPstProc, :)];
    EFTQ1Postproc = EFTQ1Postproc(:);
    EFTQ2Postproc = [3*mshQ2Ref.elements(elPstProc, :) - 2
                     3*mshQ2Ref.elements(elPstProc, :) - 1
                     3*mshQ2Ref.elements(elPstProc, :)];
    EFTQ2Postproc = EFTQ2Postproc(:);
%[text] %[text:anchor:H_34CE5031] ### Basis functions matrices at the evaluation point
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
%[text] ... for the bilinear mesh
    ShapeFunctionValuesQ2 = ...
        computeBiquadraticBasisFunctionsAndFirstDerivatives ...
        (xiEta(1, 1), xiEta(2, 1));
    NmtxQ2 = zeros(3, 3*numNodesElQ2);
    for jj = 1:numel(ShapeFunctionValuesQ2(:, 1))
        NmtxQ2(1, 3*jj-2) = ShapeFunctionValuesQ2(jj, 1);
        NmtxQ2(2, 3*jj-1) = ShapeFunctionValuesQ2(jj, 1);
        NmtxQ2(3, 3*jj) = ShapeFunctionValuesQ2(jj, 1);
    end
%[text] %[text:anchor:H_95BC169F] ### Dirichlet boundary conditions
%[text] %[text:anchor:H_6064062F] #### Dirichlet boundary conditions for the bilinear quadrilateral mesh
    numNodesQ1Ref = numel(mshQ1Ref.nodes(:, 1));
    [freeDOFsQ1Ref, ~, idx_x_ref, homDOFs_Q1_X_ref] = ...
        utilityFunctionsTwoSidedClampedPlate.getFreeDOFsBilinearQuadrilateralMesh ...
        (numNodesQ1Ref, numElxRef, numElyRef);
%[text] %[text:anchor:H_789012AE] #### Dirichlet boundary conditions for the bilinear quadrilateral mesh
    numNodesQ2Ref = numel(mshQ2Ref.nodes(:, 1));
    [freeDOFsQ2Ref, ~] = ...
        utilityFunctionsTwoSidedClampedPlate.getFreeDOFsBiquadraticQuadrilateralMesh ...
        (numNodesQ1Ref, numNodesQ2Ref, numElyRef, homDOFs_Q1_X_ref, idx_x_ref, ...
        biasXRef);
%[text] %[text:anchor:H_768F3614] ### Master stiffness matrices
%[text] %[text:anchor:H_3CC2A15D] #### Master stiffness matrix for the bilinear quadrilateral mesh
    computeStiffMatrixandForceVectorQ1 = 'undefined';
    computeBasisFunctionsAndDerivsQ1 = ...
        @computeBilinearBasisFunctionsAndFirstDerivatives;
    [KQ1Ref, FQ1Ref] = ...
        computeMasterStiffMatrixandForceVctReissnerMindlinPlatePageWise ...
        (mshQ1Ref, computeBasisFunctionsAndDerivsQ1, ...
        computeStiffMatrixandForceVectorQ1, propStr);
%[text] %[text:anchor:H_CD8AB6DD] #### Master stiffness matrix for the biquadratic quadrilateral mesh
    computeStiffMatrixandForceVectorQ2 = 'undefined';
    computeBasisFunctionsAndDerivsQ2 = ...
        @computeBiquadraticBasisFunctionsAndFirstDerivatives;
    [KQ2Ref, FQ2Ref] = ...
        computeMasterStiffMatrixandForceVctReissnerMindlinPlatePageWise ...
        (mshQ2Ref, computeBasisFunctionsAndDerivsQ2, ...
        computeStiffMatrixandForceVectorQ2, propStr);
%[text] %[text:anchor:H_1CFB0D47] #### Master stiffness matrix for the bilinear quadrilateral mesh with Assumed Natural Strain (ANS)
    computeStiffMatrixandForceVectorQ1ANS = ...
        @computeElementStiffMatrixandForceVectorReissnerMindlinPlateANS;
    computeBasisFunctionsAndDerivsQ1ANS = ...
        @computeBilinearBasisFunctionsAndFirstDerivatives;
    [KQ1ANSRef, FQ1ANSRef] = ...
        computeMasterStiffMatrixandForceVectorReissnerMindlinPlate ...
        (mshQ1Ref, computeBasisFunctionsAndDerivsQ1ANS, ...
        computeStiffMatrixandForceVectorQ1ANS, propStr);
%[text] %[text:anchor:H_D00EB33D] #### Solution of the linear equation systems
%[text] %[text:anchor:H_DFD35E49] #### Solution of the linear equation system for the bilinear quadrilateral mesh
    uQ1Ref = zeros(3*numNodesQ1Ref, 1);
    uQ1Ref(freeDOFsQ1Ref, 1) = ...
        KQ1Ref(freeDOFsQ1Ref, freeDOFsQ1Ref)\FQ1Ref(freeDOFsQ1Ref);
%[text] %[text:anchor:H_BB0AC0A7] #### Solution of the linear equation system for the biquadratic quadrilateral mesh
    uQ2Ref = zeros(3*numNodesQ2Ref, 1);
    uQ2Ref(freeDOFsQ2Ref, 1) = ...
        KQ2Ref(freeDOFsQ2Ref, freeDOFsQ2Ref)\FQ2Ref(freeDOFsQ2Ref);
%[text] %[text:anchor:H_53B3305A] #### Solution of the linear equation system for the bilinear quadrilateral mesh with Assumed Natural Strain (ANS)
    uQ1ANSRef = zeros(3*numNodesQ1Ref, 1);
    uQ1ANSRef(freeDOFsQ1Ref, 1) = ...
        KQ1ANSRef(freeDOFsQ1Ref, freeDOFsQ1Ref)\FQ1ANSRef(freeDOFsQ1Ref);
%[text] %[text:anchor:H_16046F57] ### Computation of the vertical deflections at the selected point
%[text] %[text:anchor:H_E3726A1C] #### Computation of the vertical deflection for the bilinear quadrilateral mesh
    wPostprocQ1 = NmtxQ1*uQ1Ref(EFTQ1Postproc);
    tipDeflection(ii, 1) = wPostprocQ1(1);
%[text] %[text:anchor:H_76070AB2] #### Computation of the vertical deflection for the bilinear quadrilateral mesh
    wPostprocQ2 = NmtxQ2*uQ2Ref(EFTQ2Postproc);
    tipDeflection(ii, 2) = wPostprocQ2(1);
%[text] %[text:anchor:H_313379B5] #### Computation of the vertical deflection for the bilinear quadrilateral mesh with Assumed Natural Strain (ANS)
    wPostprocQ1ANS = NmtxQ1*uQ1ANSRef(EFTQ1Postproc);
    tipDeflection(ii, 3) = wPostprocQ1ANS(1);
%[text] 
end
%%
%[text] %[text:anchor:H_032FEBD8] ## Convergence graph regarding the refinement study
loglog(numElmnts, abs(tipDeflection(:, 1)), "-bo", LineWidth=2) %[text:anchor:M_662C0768]
hold on;
loglog(numElmnts, abs(tipDeflection(:, 2)), "-ro", LineWidth=2)
loglog(numElmnts, abs(tipDeflection(:, 3)), "-go", LineWidth=2)
loglog(numElmnts, abs(tipDeflectionKirchhoff), "-k", LineWidth=2)
hold off;
grid on;
legend("Q1", "Q2", "Q1 (ANS)", "Kirchhoff theory", Location='best')
xlabel("Number of elements")
ylabel("Tip displacement")
title("Refinment study")
%[text] ![](text:image:7efd) **Try**
%[text] Change the thickness using the provided slider in code section **Problem setup** of [this](file:.\main_Chapter4_Locking_ReissnerMindlinSquarePlate.m) script and observe the convergence curves in the figure about [convergence study](internal:M_662C0768).
%[text] ![](text:image:8a24) **Reflect**
%[text] - What can you say about the convergence of the Q1-elements when the thickness is scaled down/up?
%[text] - What can you say about the convergence of the Q2-elements when the thickness is scaled up?
%[text] - What can you say about the convergence of the Q1 (ANS)-elements when the thickness is scaled up? \

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"inline","rightPanelPercent":40}
%---
%[text:image:7efd]
%   data: {"align":"baseline","height":26,"src":"data:image\/png;base64,iVBORw0KGgoAAAANSUhEUgAAAB4AAAApCAYAAAAiT5m3AAADOElEQVR42r1XSWgUQRSNGndR3Ii7V\/EmCCp4DaIe1Oh4EMHBwMjsq47Eg+OOJ0VIICAeFAIeXEADHgyoqChE1ASNhISAxH0jg5gxq+\/FX1A009PTM93T8OiuX7\/q1f+\/+tevqqoCTyKRWBmNRo8Cj4GPwF++I5HII7zT8Xh8eZWTTzgcnh6Lxc5j8kFgHMgCT4G7wBPgl8j\/YBGnfT7f1LJJg8HgQrGQE3diAXVciK5DIsi3o7+devh+4Pf755dMmslkpmGih0LaiHa1hf5k6J0T\/TYrfdMHKz8pFjTr1kFWD9ltvJ8Bt\/Dt1d0rYeG4Y7ZJMbAG+A300HLKuHkw2SuxKAe8lQ3Gdjv6llLP4\/FMQfs1MMBQ2SUOyaoPihur0X4BDGEDHVGLYbwhawCG6QG6W6yu43jo+uwSt3IyWLlA2vt09+F7CwlDodAaaZ+Q\/j3inZn0GENil\/gd0Ke1rwEjgUBgDqzYIe4lBiBbkkqlZov7r2hjXtLldom\/0LVa+x7QK9+NGjHduUvkfcAdbUwb0G+XuNtg8WV6QeLXbCDeKzo9QJM2phPosEt8ExijG9nG5NuKJK6VTTeXrofudbvEh2TS49r\/eakQMeQXNV2v9Nfbzc9c8U\/GWv2LKkmYEat+r9c7Q+L9nZvOdhLBhOF86c+MWJ5JkLVIX6jUdM1Jbsgk3NU1RRBPHCzQ2crxZR0U3CBC8saKGO3NCNNGx85kpk4QPLfaXJr8jGPkKs4FfqesJv9kcP1OvFeXa70Z8Ygmz1KWTCYX4fuDqlAQhnUVIWafQfesXjjgF5zlCrGWSBQuyOm1nuEARqGTqRgxKxejrm1iYL+URKPFEktJrGSDpRJ3aYlGYRhFwjK8T+lyle+1yrUsYjMM5ZG9l0tBzk3iYvGfWE6XVWYniyvEci1RgjEpZwP6Itwi7jfppHy31Mz3HSeGNZvw8aOAUs5h0olUqqqPFVIhjlcI34wXsMMuWWhEb76ibwPw1U1ibmCzwm+tdil3HEyfhcrdBhetbrW6pHe5RNxiVeTXukTcVEyevupCjA9YEqfT6XlyGXOKuINnQ1GnEy9y6uJeBuFn3qXxxyzmnP8Amz0wy7TM4YgAAAAASUVORK5CYII=","width":19}
%---
%[text:image:8a24]
%   data: {"align":"baseline","height":24,"src":"data:image\/png;base64,iVBORw0KGgoAAAANSUhEUgAAACgAAAAsCAYAAAAXb\/p7AAADlklEQVR42u1YS0hUYRQejR4UvR+SvaCyWiS2DIoCKWjVYzG2SsYWQ+a8HRwXFkOgqdGLSoSeqyBaJUGLCKGkNpGtgiyScnpYBpaU5au+E+eH03Wu93rvP2OEFw46955z\/u8\/7\/\/3eDQ+4XD4Eui65199AK4d1DEJcBLgfw0wGo2ujkQi1yoqKlY5AFgDqnWw5lpaMxaLrbBkDgaDm7HIEKjTCcjxPqFQqABrdYEGsHaRXaESCAyCXpNFMwwuReDw\/77xusurQGJna0ysPR2K94CnGdQKesnUyu92J5PJaelk4c51Chzcu9dp0O9nkO3yvd\/vnwqlEbzvAf2yoB5sIkwyBt0doJ+0Qbdu2AXyC8UrOVuNQH5Q3IJegfrTfH8ikwC\/A9hjsdaYCQQCG6D4vVj0DShB7+HKXMGaQ\/EFADEGrPhJdn1GAhoZvRDKX4jFTvt8vhlWchSD4D0OGiE5AH9eXl4+PxNF+IoCB8sccyBfreQB8qJWcAC0CYqHeYE2cqEDNTmQvcc6hlC+CnXWq\/PCtTtdeGGr0HNWFz7aeTcr\/eT1eqe41PVWJYxd62zkovwXqQ6Cv8vErm9r8MZNEctLVf9Ph4GwecSOjPSQu8V28e6chmRrFPq28btHJhhSBKDIBH0B726HEDilAWCtsGCx6MmjMNgaHqBkiygPl3UmHE1OrjOEiqoqsqDHGizYxrqGAXCOriL9jJUOoq3lO9UDQItpQGD3PtXZRapFHB51oSchwiWsDWA8Hl8Cpb2svA+WWO4AXJ7Q0Z1IJObq7sUBN\/UQLr0hrHfANSCaQihrobiBf+dC+QOxSNk4MrdEyN1VvRy6T+DdVbOp2wpcCyutN8yDahjtVZ3A4sS2ALwflYxhaK1n0C22QfL8dkuNRYZBlKxxRLj6gg3XNgj+w8beLOriHcv50gocJ8wsfP\/ASgfMDlRs8XzwfGPeThMr2QcpEqFprJmPSoSwSvMYfCcFX6nFlNPEfAGrUuC1GkhxQpspTnTfKysrFxl5qqqqZouy0mU81ZmApB6cp2vKlrFVk+Z7TGRuPOsXRZSNfGYmEO\/ofoXmOiKKS3Ga6wPYeRNxk9Vo49BONAIL1mUdIBa9bxPgn+zMOkAqIXDdQbp1sKAy6uWeyWd0Cyuk5LDh3pTtez\/NSVJqNwbh5kNZB8htkTrP1zHA9dM1nZ07nExmc50ZQGA7M+GxyLden9MA\/OLm\/KKe34xPoKAD7jPkAAAAAElFTkSuQmCC","width":22}
%---
