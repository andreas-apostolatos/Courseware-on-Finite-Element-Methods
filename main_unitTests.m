%[text] %[text:anchor:T_2F2203F8] # Unit Tests
%[text:tableOfContents]{"heading":"Table of Contents"}
%[text] Set of [unit tests](https://www.mathworks.com/help/matlab/matlab_prog/ways-to-write-unit-tests.html) to safeguard functionalities provided in this courseware (numerical quadrature, basis functions, stiffness matrices, load vectors, assembly, etc.)
%[text] 
%[text] **Author: Andreas Apostolatos, PhD (aapostol@mathworks.com)**
%[text] **Date: 20.02.2023**
%%
%[text] %[text:anchor:H_B23D1887] ## Importing the [Unit Test Suite of MATLAB](https://www.mathworks.com/help/matlab/ref/matlab.unittest.testsuite-class.html) and accompanying modules
import matlab.unittest.TestSuite;
import matlab.unittest.selectors.HasName;
import matlab.unittest.constraints.EndsWithSubstring;
%%
%[text] %[text:anchor:H_4DA5D08F] ## Run the unit tests
%%
%[text] %[text:anchor:H_E27F0E3B] ### Run the unit tests for the [Gauss Quadrature](https://en.wikipedia.org/wiki/Gaussian_quadrature) numerical integration
suiteQuadrature = TestSuite.fromClass(?testQuadrature);
resultQuadrature = run(suiteQuadrature); %[output:39c621dd]
%%
%[text] %[text:anchor:H_71140377] ### Run the unit tests for the [basis functions](https://en.wikipedia.org/wiki/Finite_element_method#Choosing_a_basis) used for the Finite Element formulations
suiteBasisFunctions = TestSuite.fromClass(?testBasisFunctions);
resultBasisFunctions = run(suiteBasisFunctions); %[output:7026909e]
%%
%[text] %[text:anchor:H_9E1ACF5A] ### Run the unit tests for the Finite Element formulation of the [Timoshenko beam](https://en.wikipedia.org/wiki/Timoshenko%E2%80%93Ehrenfest_beam_theory) problem
suiteTimoshenkoBeam = TestSuite.fromClass(?testTimoshenkoBeam);
resultTimoshenkoBeam = run(suiteTimoshenkoBeam); %[output:096a6d35]
%%
%[text] %[text:anchor:H_2B1F3E9F] ### Run the unit tests for the Finite Element formulation of the [Reissner-Mindlin plate](https://en.wikipedia.org/wiki/Mindlin%E2%80%93Reissner_plate_theory) problem
suiteReissnerMindlinPlate = TestSuite.fromClass(?testReissnerMindlinPlate);
resultReissnerMindlinPlate = run(suiteReissnerMindlinPlate); %[output:5ed9be8f]

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"inline","rightPanelPercent":26.7}
%---
%[output:39c621dd]
%   data: {"dataType":"text","outputData":{"text":"Running testQuadrature\n..\nDone testQuadrature\n__________\n\n","truncated":false}}
%---
%[output:7026909e]
%   data: {"dataType":"text","outputData":{"text":"Running testBasisFunctions\n..\nDone testBasisFunctions\n__________\n\n","truncated":false}}
%---
%[output:096a6d35]
%   data: {"dataType":"text","outputData":{"text":"Running testTimoshenkoBeam\n...\nDone testTimoshenkoBeam\n__________\n\n","truncated":false}}
%---
%[output:5ed9be8f]
%   data: {"dataType":"text","outputData":{"text":"Running testReissnerMindlinPlate\n...\nDone testReissnerMindlinPlate\n__________\n\n","truncated":false}}
%---
