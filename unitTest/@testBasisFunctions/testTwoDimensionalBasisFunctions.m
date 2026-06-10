%[text] %[text:anchor:T_FC84E213] # Test Two-Dimensional Basis Functions
%[text:tableOfContents]{"heading":"**Table of Contents**"}
%[text] %[text:anchor:H_7E81E991] ## Brief summary of this function
%[text] Tests the bilinear and the biquadratic two-dimensional basis functions which are defined in the bi-unit space $\\textbf{\\xi} = \\left(\\xi .\\eta \\right) \\in \\left\[ -1, 1 \\right\]^2$
function testTwoDimensionalBasisFunctions(testCase)
%[text] %[text:anchor:H_51D4211A] ## Definition of tolerances
    tol = 1e-14;
    tol1 = 1e-14*1e1;
%[text] %[text:anchor:H_A4ABB2E3] ## Test the bilinear two-dimensional basis functions
%[text] %[text:anchor:H_600C14C4] ### Define the expected solution
    expdN_Q1 = [   0.193626775000000  -0.335575000000000  -0.144250000000000
                   0.094873225000000  -0.164425000000000   0.144250000000000
                   0.233976775000000   0.164425000000000   0.355750000000000
                   0.477523225000000   0.335575000000000  -0.355750000000000];
%[text] %[text:anchor:H_5F12C8B8] ### Compute the bilinear two-dimensional basis functions
    dN_Q1 = computeBilinearBasisFunctionsAndFirstDerivatives(0.423,-0.3423);
%[text] %[text:anchor:H_EE4833AD] ### Verify the solution
    testCase.verifyEqual(dN_Q1, expdN_Q1, "RelTol", tol);
%[text] %[text:anchor:H_FF87B8E5] ## Test the biquadratic two-dimensional basis functions
%[text] %[text:anchor:H_E7D4AFA2] ### Define the expected solution
    expdN_Q2 =    [   -0.000075436737750   0.000371128500000   0.075361225500000 %[text:anchor:H_5F91C196]
                       0.150798038762250  -0.741885871500000   0.226385725500000
                      -0.091837582737750   0.256614628500000  -0.137871274500000
                       0.000045941762250  -0.000128371500000  -0.045895774500000
                       0.000301897975500  -0.001485257000000  -0.301746951000000
                      -0.000183859024500   0.000513743000000   0.183767049000000
                      -0.000470005024500  -0.000242757000000   0.469534549000000
                       0.939540043975500   0.485271243000000   1.410485549000000
                       0.001880961049000   0.000971514000000  -1.880020098000000];
%[text] %[text:anchor:H_D3B70F7B] ### Compute the biquadratic two-dimensional basis functions
    dN_Q2 = computeBiquadraticBasisFunctionsAndFirstDerivatives(-0.243, 0.999);
%[text] %[text:anchor:H_45A178E9] ### Verify the solution
    testCase.verifyEqual(dN_Q2, expdN_Q2, "RelTol", tol1);
%[text] 
end

%[appendix]{"version":"1.0"}
%---
