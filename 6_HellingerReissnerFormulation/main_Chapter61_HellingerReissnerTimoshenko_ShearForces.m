%[text] %[text:anchor:H_4E46353A] # Transverse-shear locking due incompatible discretization spaces
%[text:tableOfContents]{"heading":"Table of Contents"}
%[text] The matter of transverse-shear locking and its presence in Finite Element formulation for shear-deformable structural elements is discussed in detail in Live Script [`main_Chapter4_TransverseShearLocking_TimoshenkoBeam.m`](file:..\4_TransverseShearLocking\4_TransverseShearLocking_TimoshenkoBeam\main_Chapter4_TransverseShearLocking_TimoshenkoBeam.m). As mentioned therein the issue is the incompatibility of the discretization spaces for the displacement and the cross-sectional rotation fields because of the use of the same basis functions for the discretization of both fields. If it is used linear shape functions then the transverse-shear deformation cannot represent zero-states by simply adding piecewise constant with piecewise linear functions, namely:
%[text] $\\gamma = \\frac{\\text{d}w}{\\text{d}x} + \\beta = \\sum\_{i=1}^n \\overbrace{\\frac{\\text{d}N\_i^w(\\xi)}{\\text{d}x}}^{\\text{constant}} w\_i + \\sum\_{i = 1}^n \\overbrace{N\_i^{\\beta}(\\xi)}^{\\text{linear}} \\beta\_i \\qquad (2)$
%[text] This leads to the fact that residual spurious transverse-shear forces remain that eventually cause the transverse-shear locking numerical phenomenon. Another way to identify/explain transverse-shear locking is that the problem is overconstrained by means of the Kirchhoff constrain $\\gamma \\overbrace{=}^{!} 0$ in case of thin shear-deformable elements. A way to circumvent this issue to avoid in the first place having to incorporate this contrain in the variational form of the problem that describes the Partial Differential Equation (PDE) weakly. In the case of the Timoshenko beam one does not need to introduce the transverse-shear deformations into the variational form of the problem and leverage the inverse of the material law. The underlying derivation leads in this case to the following variational formulation leveraging the Hellinger-Reissner principle:
%[text] %[text:anchor:H_32BAF93F] $\\int\_0^L \\delta \\left( \\frac{\\text{d} w}{\\text{d} x} + \\beta \\right) \\, Q  \\; \\text{d} x - \\delta Q \\, \\frac{1}{\\alpha G A} \\, Q + \\delta Q \\, \\left( \\frac{\\text{d} w}{\\text{d} x} + \\beta \\right) + \\delta \\frac{\\text{d} \\beta}{\\text{d} x} \\, M - \\delta M \\, \\frac{1}{EI} \\, M + \\delta M \\, \\frac{\\text{d} \\beta}{\\text{d} x} = \\overbrace{\\int\_0^L \\overline{q} \\, \\delta w + \\overline{m} \\, \\delta \\beta \\; \\text{d}x }^{\\text{source loads}} + \\overbrace{\\left\[ \\delta w \\, Q \\right\]\_0^L + \\left\[ \\delta \\beta \\, M \\right\]\_0^L}^{\\text{boundary loads}} \\;. \\quad (3)$
%[text] The latter variational formulation is proven to be transverse-shear locking-free, however additional challenges are posed:
%[text] - The system involves additional unknowns (in this case the Lagrange Multipliers) which can be condensed out in particular situations
%[text] - The system is no longer necessarily positive definite meaning that the discretization spaces need to fullfil the [Ladyzenskaja-Babuska-Brezzi (LBB) (or alternative known as inf-sup) condition](https://en.wikipedia.org/wiki/Ladyzhenskaya%E2%80%93Babu%C5%A1ka%E2%80%93Brezzi_condition) so that a unique solution to the problem can be guaranteed \
%[text] 
%[text] Go to [THIS](file:.\main_Chapter6_HellingerReissnerFormulation_TimoshenkoBeam_Study.m) Section in the main driver script *Transverse-shear locking due incompatible discretization spaces*
%[text] Go to [NEXT](file:.\main_Chapter62_HellingerReissnerTimoshenko_ConvergenceStudy.m) Section *Convergence study on transverse shear-locking*
%[text] 
%[text] **Author: Andreas Apostolatos, PhD (aapostol@mathworks.com)**
%[text] **Date: 20.02.2023**
%%
%[text] %[text:anchor:H_7236A159] ## Preamble
vars = ["numPts" "mshL1" "mshL2" "numElL1" "NL1" "NL2" ...
    "propStr" "choiceShapeFunctionsU" "choiceShapeFunctionsLM"];
str = "This mlx-file should be run from parent mlx-file " + ...
    "main_Chapter5_HellingerReissnerFormulation_TimoshenkoBeam_Study.mlx";
for ii = 1:length(vars)
    if ~exist(vars(ii), 'var') 
        error(strcat("Variable ", vars(ii), " is not defined. ", str));
    end
end
%%
%[text] %[text:anchor:H_376C7AF5] ## Initialization of the evaluation points along the beam where the shear forces are computed
%[text] Definition of the vertical deflection $w(\\xi)$, rotation $\\beta$, position $X(\\xi)$, jacobian $J (\\xi) = \\frac{\\text{d} X}{\\text{d} \\xi}$, the derivative of the vertical deflection with respect to the physical space $\\frac{\\partial w (\\xi)}{\\partial X}$ and shear force $Q(\\xi)$ function fields
syms wL1(xi) wL1RI(xi) wL2(xi) betaL1(xi) betaL1RI(xi) betaL2(xi) ...
    posVctL1(xi) posVctL2(xi) jacobianMtxL1(xi) jacobianMtxL2(xi) ...
    dwdxL1(xi) dwdxL1RI(xi) dwdxL2(xi) QL1(xi) QL1RI(xi) QL2(xi) ...
    phi1(xi) phi2(xi) phi3(xi) phi4(xi) Phi1(X) Phi2(X) Phi3(X) ...
    Phi4(X) Phi(X) wEB(X) QEB(X)
%[text] Compute the bounding box of the Finite Element mesh and generate a sequence of sampling points onto which to compute the shear forces:
boundingBox = [min(mshL1.nodes(:, 1)) max(mshL1.nodes(:, 1))];
xVct = linspace(boundingBox(1), boundingBox(2), numPts);
QPts = zeros(numPts, 5);
%[text] %[text:anchor:H_7F73DFF9] ## Loop over all the sampling points for the computation of the transverse-shear forces
for ii = 1:numPts
%[text] %[text:anchor:H_A61AB742] ### Element where the evaluation point lies
    el = nan;
    for jj = 1:numElL1
%[text] %[text:anchor:H_F91C4FB2] ### Set-up the position vector corresponding
%[text] The position vector is given by the linear combination $x(\\xi) = \\sum\_{i = 1}^2 N\_{i}(\\xi) x\_i^{\\text{e}}$ for both the linear and the bilinear mesh
        posVctL1(xi) = NL1*mshL1.nodes(mshL1.elements(jj, :));
        posVctL2(xi) = NL2*mshL2.nodes(mshL2.elements(jj, :));
%[text] %[text:anchor:H_6022508A] ### Inverse position vector corresponding to the current element
%[text] The inverse mapping can be constructed in the one-dimensional case analytically as $\\xi(x) = x(\\xi)^{-1}$ **Note** This is not possible for higher dimensions, where a matrix inversion with nonlinear entries might be needed
        posVctInv = finverse(posVctL1);
%[text] %[text:anchor:H_8EDCC784] ### Parametric coordinate of the evaluation point in the current element using the element's inverse parametric representation
        xiParam = posVctInv(xVct(ii));
%[text] %[text:anchor:H_D1B377FD] ### Verify whether the computed parametric coordinate lies in the parametric space
%[text] Verify whether the computed parametric $\\xi$ lies within the closed space $\[-1,1\]$. If it is, then the evaluation point belongs to the current element, otherwise not
        if xiParam >= -1 && xiParam <= 1
            el = jj;
            break;
        end
%[text] 
    end
%[text] %[text:anchor:H_544B08E0] ### Verification of the element's validity containing the current point
    if isnan(el)
        error("The point with coordinate %d was found in no element " + ...
            "in the mesh", x(ii));
    end
%[text] %[text:anchor:H_30FD46CC] ### Element Freedom Table (EFT) for the one-dimensional linear element
    EFTL1 = [2*mshL1.elements(el, :) - 1
             2*mshL1.elements(el, :)];
    EFTL1 = EFTL1(:);
%[text] %[text:anchor:H_952478D3] ### Element Freedom Table (EFT) for the two-dimensional quadratic element
    EFTL2 = [2*mshL2.elements(el, :) - 1
             2*mshL2.elements(el, :)];
    EFTL2 = EFTL2(:);
%[text] %[text:anchor:H_159A0535] ### Vertical deflection and rotation fields for the one-dimensional linear element
    wL1(xi) = NL1*uL1(EFTL1(1:2:end), 1); % w(x)
    betaL1(xi) = NL1*uL1(EFTL1(2:2:end), 1); % β(x)
%[text] %[text:anchor:H_54DA1CA6] ### Vertical deflection and rotation fields for the one-dimensional linear element with the selective-reduced integration
    wL1RI(xi) = NL1*uL1RI(EFTL1(1:2:end), 1); % w(x)
    betaL1RI(xi) = NL1*uL1RI(EFTL1(2:2:end), 1); % β(x)
%[text] %[text:anchor:H_6B18C1D6] ### Vertical deflection and rotation fields for the one-dimensional quadratic element
    wL2(xi) = NL2*uL2(EFTL2(1:2:end), 1); % w(x)
    betaL2(xi) = NL2*uL2(EFTL2(2:2:end), 1); % β(x)
%[text] %[text:anchor:H_CD409334] ### Derivatives of the vertical deflection and rotation fields with respect to the physical space
%[text] To compute the derivatives of the basis functions with respect to the physical space given that they are defined in the parametric space, the chain rule in one-dimensiona can be used, namely:
%[text] $J^{-1}(\\xi) = \\frac{\\text{d} \\xi}{\\text{d} x}$
%[text] $\\frac{\\text{d} w}{\\text{d} x} = \\frac{\\text{d} w}{\\text{d} \\xi} \\overbrace{\\frac{\\text{d} \\xi}{\\text{d} x}}^{J^{-1}}$
%[text] i.e. it is simply a scalar quantity. The Jacobian matrix in this case can be computed using function [`jacobian`](https://www.mathworks.com/help/symbolic/sym.jacobian.html) in MATLAB:
    jacobianMtxL1 = jacobian(posVctL1, xi);
    jacobianMtxL2 = jacobian(posVctL2, xi);
%[text] ... setting-up the derivatives of the vertical deflection field corresponding to the one-dimensional linear elements based on full integration with respect to the physical space
    dwdxL1(xi) = transpose(inv(jacobianMtxL1(xi)))*diff(wL1, xi);
%[text] ... setting-up the derivatives of the vertical deflection field corresponding to the one-dimensional linear elements based on selective-reduced integration with respect to the physical space
    dwdxL1RI(xi) = transpose(inv(jacobianMtxL1(xi)))*diff(wL1RI, xi);
%[text] ... setting-up the derivatives of the vertical deflection field corresponding to the one-dimensional quadratic elements based on full integration with respect to the physical space
    dwdxL2(xi) = transpose(inv(jacobianMtxL2(xi)))*diff(wL2, xi);
%[text] %[text:anchor:H_3605C384] ### Resultant transverse-shear force for the standard Timoshenko beam Finite Element formulation at the evaluation point
%[text] $Q (\\xi) = \\alpha G A \\gamma (\\xi) = \\alpha G A \\overbrace{ \\left( \\frac{\\text{d} w (\\xi)}{\\text{d} x} + \\beta(\\xi)\\right)}^{\\gamma(\\xi)}$
%[text] ... corresponding to the one-dimensional linear elements based on full integration with respect to the physical space
    QL1(xi) = propStr.alpha*propStr.G*propStr.A*(dwdxL1(xi) + ...
        betaL1(xi));
    QPts(ii, 1) = double(QL1(xiParam));
%[text] ... corresponding to the one-dimensional linear elements based on selective-reduced integration with respect to the physical space
    QL1RI(xi) = propStr.alpha*propStr.G*propStr.A*(dwdxL1RI(xi) + ...
        betaL1RI(xi));
    QPts(ii, 2) = double(QL1RI(xiParam));
%[text] ... corresponding to the one-dimensional quadratic elements based on full integration with respect to the physical space
    QL2(xi) = propStr.alpha*propStr.G*propStr.A*(dwdxL2(xi) + ...
        betaL2(xi));
    QPts(ii, 3) = double(QL2(xiParam));
%[text] ... corresponding to the Hellinger-Reissner formulation
    QPts(ii, 4) = QPtsL01HR(el, 1);
%[text] %[text:anchor:H_57C315ED] ### Displacement field at the element level with respect to the physical space using the cubic Hermitian basis functions
%[text] $w(X) = N\_1(X)w\_i + (X\_j - X\_i)N\_2(X)\\phi\_i + N\_3(X)w\_j + (X\_j - X\_i)N\_4(X) \\phi\_j$
    phi1(xi) = 2*xi^3 - 3*xi^2 + 1;
    phi2(xi) = diff(mshL1.nodes(mshL1.elements(el, :)))*...
        (xi.^3 - 2*xi.^2 + xi);
    phi3(xi) = -2*xi.^3 + 3*xi.^2;
    phi4(xi) = diff(mshL1.nodes(mshL1.elements(el, :)))*(xi.^3 - xi.^2);
    Phi1(X) = subs(phi1, xi, (X - mshL1.nodes(mshL1.elements(el, 1)))/ ...
        diff(mshL1.nodes(mshL1.elements(el, :))));
    Phi2(X) = subs(phi2, xi, (X - mshL1.nodes(mshL1.elements(el, 1)))/ ...
        diff(mshL1.nodes(mshL1.elements(el, :))));
    Phi3(X) = subs(phi3, xi, (X - mshL1.nodes(mshL1.elements(el, 1)))/ ...
        diff(mshL1.nodes(mshL1.elements(el, :))));
    Phi4(X) = subs(phi4, xi, (X - mshL1.nodes(mshL1.elements(el, 1)))/ ...
        diff(mshL1.nodes(mshL1.elements(el, :))));
    Phi(X) = [Phi1(X) Phi2(X) Phi3(X) Phi4(X)];
    wEB(X) = Phi(X)*uEB(EFTL1, 1);
%[text] %[text:anchor:H_511D846A] ### Resultant transverse-shear force for the Euler-Bernoulli beam Finite Element formulation at the evaluation point
%[text] $Q = -\\frac{\\text{d} }{\\text{d} x} \\left( E I \\frac{\\text{d}^2 w}{\\text{d} x^2} \\right)$
    QEB(X) = - propStr.E*propStr.I * diff(wEB(X), X, 3);
    QPts(ii, 5) = double(QEB(xiParam));
%[text] 
end
%%
%[text] %[text:anchor:H_7CF6B51F] ## Shear force distribution along the beam's length
plotMesh1d(mshL1);
hold on;
plot(xVct, QPts(:, 1), '-bo', LineWidth=2);
plot(xVct, QPts(:, 2), '-bv', LineWidth=2);
plot(xVct, QPts(:, 3), '-ro', LineWidth=2);
plot(xVct, QPts(:, 4), '-m^', LineWidth=2);
plot(xVct, QPts(:, 5), '-go', LineWidth=2);
fplot(QEx, [X0, XL], '-k', LineWidth=2)
legend("", "L1", "L1 (RI)", "L2", "L" + num2str(choiceShapeFunctionsU) + ...
    '/' + num2str(choiceShapeFunctionsLM) + " H/R", "E/B", "analytical", ...
    Location="best");
hold off;
grid on;
xlabel("x")
ylabel("y")
zlabel("q_x")
title("Shear forces Q along the beam")
%[text] ![](text:image:822a) **Try**
%[text] Change the thickness of the beam (and eventually its slenderness) using the numeric slider in code section **Problem parameters** of [this](file:.\main_Chapter6_HellingerReissnerFormulation_TimoshenkoBeam_Study.m) script. Try the extreme values provided in the slider
%[text] ![](text:image:8fcf) **Reflect**
%[text] - How do the numerical solutions corresponding to the Timosheko beam theory compare to the numerical solution of the Euler-Bernoulli beam theory for thick beams?
%[text] - How do the resultant transverse shear forces behave for the linear Finite Element discretization for thin beams?
%[text] - How do the resultant transverse shear forces behave for the linear Finite Element discretization for thin beams when using selective-reduced integration?
%[text] - What can you comment about the distribution of the resultant transverse-shear forces when using the Hellinger-Reissner principle as compared to the equivalent distribution when using the Euler-Bernoulli beam theory for thin beams? \

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"inline","rightPanelPercent":40}
%---
%[text:image:822a]
%   data: {"align":"baseline","height":26,"src":"data:image\/png;base64,iVBORw0KGgoAAAANSUhEUgAAAB4AAAApCAYAAAAiT5m3AAADOElEQVR42r1XSWgUQRSNGndR3Ii7V\/EmCCp4DaIe1Oh4EMHBwMjsq47Eg+OOJ0VIICAeFAIeXEADHgyoqChE1ASNhISAxH0jg5gxq+\/FX1A009PTM93T8OiuX7\/q1f+\/+tevqqoCTyKRWBmNRo8Cj4GPwF++I5HII7zT8Xh8eZWTTzgcnh6Lxc5j8kFgHMgCT4G7wBPgl8j\/YBGnfT7f1LJJg8HgQrGQE3diAXVciK5DIsi3o7+devh+4Pf755dMmslkpmGih0LaiHa1hf5k6J0T\/TYrfdMHKz8pFjTr1kFWD9ltvJ8Bt\/Dt1d0rYeG4Y7ZJMbAG+A300HLKuHkw2SuxKAe8lQ3Gdjv6llLP4\/FMQfs1MMBQ2SUOyaoPihur0X4BDGEDHVGLYbwhawCG6QG6W6yu43jo+uwSt3IyWLlA2vt09+F7CwlDodAaaZ+Q\/j3inZn0GENil\/gd0Ke1rwEjgUBgDqzYIe4lBiBbkkqlZov7r2hjXtLldom\/0LVa+x7QK9+NGjHduUvkfcAdbUwb0G+XuNtg8WV6QeLXbCDeKzo9QJM2phPosEt8ExijG9nG5NuKJK6VTTeXrofudbvEh2TS49r\/eakQMeQXNV2v9Nfbzc9c8U\/GWv2LKkmYEat+r9c7Q+L9nZvOdhLBhOF86c+MWJ5JkLVIX6jUdM1Jbsgk3NU1RRBPHCzQ2crxZR0U3CBC8saKGO3NCNNGx85kpk4QPLfaXJr8jGPkKs4FfqesJv9kcP1OvFeXa70Z8Ygmz1KWTCYX4fuDqlAQhnUVIWafQfesXjjgF5zlCrGWSBQuyOm1nuEARqGTqRgxKxejrm1iYL+URKPFEktJrGSDpRJ3aYlGYRhFwjK8T+lyle+1yrUsYjMM5ZG9l0tBzk3iYvGfWE6XVWYniyvEci1RgjEpZwP6Itwi7jfppHy31Mz3HSeGNZvw8aOAUs5h0olUqqqPFVIhjlcI34wXsMMuWWhEb76ibwPw1U1ibmCzwm+tdil3HEyfhcrdBhetbrW6pHe5RNxiVeTXukTcVEyevupCjA9YEqfT6XlyGXOKuINnQ1GnEy9y6uJeBuFn3qXxxyzmnP8Amz0wy7TM4YgAAAAASUVORK5CYII=","width":19}
%---
%[text:image:8fcf]
%   data: {"align":"baseline","height":24,"src":"data:image\/png;base64,iVBORw0KGgoAAAANSUhEUgAAACgAAAAsCAYAAAAXb\/p7AAADlklEQVR42u1YS0hUYRQejR4UvR+SvaCyWiS2DIoCKWjVYzG2SsYWQ+a8HRwXFkOgqdGLSoSeqyBaJUGLCKGkNpGtgiyScnpYBpaU5au+E+eH03Wu93rvP2OEFw46955z\/u8\/7\/\/3eDQ+4XD4Eui65199AK4d1DEJcBLgfw0wGo2ujkQi1yoqKlY5AFgDqnWw5lpaMxaLrbBkDgaDm7HIEKjTCcjxPqFQqABrdYEGsHaRXaESCAyCXpNFMwwuReDw\/77xusurQGJna0ysPR2K94CnGdQKesnUyu92J5PJaelk4c51Chzcu9dp0O9nkO3yvd\/vnwqlEbzvAf2yoB5sIkwyBt0doJ+0Qbdu2AXyC8UrOVuNQH5Q3IJegfrTfH8ikwC\/A9hjsdaYCQQCG6D4vVj0DShB7+HKXMGaQ\/EFADEGrPhJdn1GAhoZvRDKX4jFTvt8vhlWchSD4D0OGiE5AH9eXl4+PxNF+IoCB8sccyBfreQB8qJWcAC0CYqHeYE2cqEDNTmQvcc6hlC+CnXWq\/PCtTtdeGGr0HNWFz7aeTcr\/eT1eqe41PVWJYxd62zkovwXqQ6Cv8vErm9r8MZNEctLVf9Ph4GwecSOjPSQu8V28e6chmRrFPq28btHJhhSBKDIBH0B726HEDilAWCtsGCx6MmjMNgaHqBkiygPl3UmHE1OrjOEiqoqsqDHGizYxrqGAXCOriL9jJUOoq3lO9UDQItpQGD3PtXZRapFHB51oSchwiWsDWA8Hl8Cpb2svA+WWO4AXJ7Q0Z1IJObq7sUBN\/UQLr0hrHfANSCaQihrobiBf+dC+QOxSNk4MrdEyN1VvRy6T+DdVbOp2wpcCyutN8yDahjtVZ3A4sS2ALwflYxhaK1n0C22QfL8dkuNRYZBlKxxRLj6gg3XNgj+w8beLOriHcv50gocJ8wsfP\/ASgfMDlRs8XzwfGPeThMr2QcpEqFprJmPSoSwSvMYfCcFX6nFlNPEfAGrUuC1GkhxQpspTnTfKysrFxl5qqqqZouy0mU81ZmApB6cp2vKlrFVk+Z7TGRuPOsXRZSNfGYmEO\/ofoXmOiKKS3Ga6wPYeRNxk9Vo49BONAIL1mUdIBa9bxPgn+zMOkAqIXDdQbp1sKAy6uWeyWd0Cyuk5LDh3pTtez\/NSVJqNwbh5kNZB8htkTrP1zHA9dM1nZ07nExmc50ZQGA7M+GxyLden9MA\/OLm\/KKe34xPoKAD7jPkAAAAAElFTkSuQmCC","width":22}
%---
