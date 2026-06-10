%[text] %[text:anchor:T_FC84E213] # Test Quadrature Over Unit Domain
%[text:tableOfContents]{"heading":"**Table of Contents**"}
%[text] %[text:anchor:H_7E81E991] ## Brief summary of this function
%[text] Tests the Gauss quadrature over a triangle which is constructed using the degenerated quadrilateral for the function $f(\\zeta\_1, \\zeta\_2) = \\ln (\\zeta\_1 + 1)\\sin \\left( \\frac{\\pi}{2} \\zeta\_2 \\right)$
function testQuadratureOverUnitDomain(testCase)
%[text] %[text:anchor:H_EB4BCEAF] ## Definition of the integrand
f = @(x) sin(exp(x));
%[text] %[text:anchor:H_03232FFC] ## Definition of the expected solution
expSol = 1.455915572116364;
%[text] %[text:anchor:H_B74EC61F] ## Definition of tolerances
tol = 1e-1;
tolUlt = 1e-15;
%[text] %[text:anchor:H_9C45AF27] ## Loop over all possible number of Gauss points for the Gauss integration over a unit domain
for iNumGP = 1:50
%[text] %[text:anchor:H_ABAB6BC3] ### Initialize the linear combination of the Gauss point coordinates with their weights
    sol = 0;
%[text] %[text:anchor:H_F65A27CD] ### Issue the Gauss point coordinates and weights
    [GP, GW] = getGaussPointsAndWeightsOverUnitDomain(iNumGP);
%[text] %[text:anchor:H_FDF30409] ### Loop over all Gauss points and compute the integral
    for iGP = 1:height(GP)
        sol = sol + f(GP(iGP))*GW(iGP);
    end
%[text] %[text:anchor:H_A77B03C6] ### Verify the solution
    tolStricktened = max(horzcat(tol*10^(-iNumGP + 2), tolUlt));
    testCase.verifyEqual(sol, expSol, "AbsTol", tolStricktened);
%[text] 
end
%[text] 
end

%[appendix]{"version":"1.0"}
%---
