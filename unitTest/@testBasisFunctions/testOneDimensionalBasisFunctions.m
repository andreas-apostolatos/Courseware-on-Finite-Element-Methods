%[text] %[text:anchor:T_FC84E213] # Test One-Dimensional Basis Functions
%[text:tableOfContents]{"heading":"**Table of Contents**"}
%[text] %[text:anchor:H_7E81E991] ## Brief summary of this function
%[text] Tests the constant, the linear and the quadratic one-dimensional basis functions which are defined in the unit space $\\xi \\in \\left\[ -1, 1 \\right\]$
function testOneDimensionalBasisFunctions(testCase)
%[text] %[text:anchor:H_51D4211A] ## Definition of tolerance
    tol = 1e-14;
%[text] %[text:anchor:H_8AAC13A6] ## Test the constant one-dimensional basis functions
%[text] %[text:anchor:H_B2C4AB0A] ### Define the expected solution
    expdN_L0 = [1 0];
%[text] %[text:anchor:H_3DC3AFC6] ### Compute the constant one-dimensional basis functions
    dN_L0 = computeConstantBasisFunctionAndFirstDerivatives(0.3123);
%[text] %[text:anchor:H_A77B03C6] ### Verify the solution
    testCase.verifyEqual(dN_L0, expdN_L0, "RelTol", tol);
%[text] %[text:anchor:H_A4ABB2E3] ## Test the linear one-dimensional basis functions
%[text] %[text:anchor:H_600C14C4] ### Define the expected solution
    expdN_L1 = [0.782250000000000  -0.500000000000000
                0.217750000000000   0.500000000000000];
%[text] %[text:anchor:H_5F12C8B8] ### Compute the linear one-dimensional basis functions
    dN_L1 = computeLinearBasisFunctionsAndFirstDerivatives(-0.5645);
%[text] %[text:anchor:H_EE4833AD] ### Verify the solution
    testCase.verifyEqual(dN_L1, expdN_L1, "RelTol", tol);
%[text] %[text:anchor:H_FF87B8E5] ## Test the quadratic one-dimensional basis functions
%[text] %[text:anchor:H_E7D4AFA2] ### Define the expected solution
    expdN_L2 =  [0.000615756450000  -0.501230000000000 %[text:anchor:H_5F91C196]
                -0.000614243550000   0.498770000000000
                 0.999998487100000   0.002460000000000];
%[text] %[text:anchor:H_D3B70F7B] ### Compute the quadratic one-dimensional basis functions
    dN_L2 = computeQuadraticBasisFunctionsAndFirstDerivatives(-0.00123);
%[text] %[text:anchor:H_45A178E9] ### Verify the solution
    testCase.verifyEqual(dN_L2, expdN_L2, "RelTol", tol);
%[text] 
end

%[appendix]{"version":"1.0"}
%---
