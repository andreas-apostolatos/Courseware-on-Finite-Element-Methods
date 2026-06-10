%[text] %[text:anchor:T_E97F0E8D] # Test Quadrature Over Triangle
%[text:tableOfContents]{"heading":"**Table of Contents**"}
%[text] %[text:anchor:H_30C195BF] ## Brief summary of this function
%[text] Tests the Gauss quadrature over a triangle which is constructed using the degenerated quadrilateral for the function $f(\\zeta\_1, \\zeta\_2) = \\ln (\\zeta\_1 + 1)\\sin \\left( \\frac{\\pi}{2} \\zeta\_2 \\right)$
function testQuadratureOverTriangle(testCase)
%[text] %[text:anchor:H_F4BEBCC7] ## Definition of the expected solution
expSol = 0.050909793169365;
%[text] %[text:anchor:H_BD929E27] ## Definition of the integrand
f = @(x,y) log(x + 1)*sin(pi*y/2);
%[text] %[text:anchor:H_DE1EB7BC] ## Definition of tolerances
tol = 1e-1;
%[text] %[text:anchor:H_5B2D15F4] ## Loop over all polynomial orders for the Gauss integration over a triangle using the symmetric rule
for iPolOrder = 1:8
%[text] %[text:anchor:H_48CC629B] ### Initialize the linear combination of the Gauss point coordinates with their weights
    sol = 0;
%[text] %[text:anchor:H_D48DD02F] ### Issue the Gauss point coordinates and weights
    [GP, GW] = getGaussPointsAndWeightsOnCanonicalTriangle(iPolOrder);
%[text] %[text:anchor:H_CDE5FAA2] ### Loop over all Gauss points and compute the integral
    for iGP = 1:height(GP)
        sol = sol + f(GP(iGP, 1), GP(iGP, 2))*GW(iGP, 1);
    end
%[text] %[text:anchor:H_1B4CF589] ### Verify the solution
    tolStricktened = tol*10^(-iPolOrder + 2);
    testCase.verifyEqual(sol, expSol, "AbsTol", tolStricktened);
%[text] 
end
%[text] 
end

%[appendix]{"version":"1.0"}
%---
