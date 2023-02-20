classdef testReissnerMindlinPlate < matlab.unittest.TestCase
%% Class definition
%
% Test suites for the Finite Element formulation of the Reissner-Mindlin
% plate problem
% - Tests the computation of the element stiffness matrices
% - Tests the assembly of the element stiffness matrices to the master 
%   stiffness matrix

properties
end

%% Method definitions
methods (Test)
    testComputationElStiffMtxLoadVctReissnerMindlin(testCase)
    testAssemblyMasterStiffMtxLoadVctReissnerMindlin(testCase)
    testMasterStiffMtxLoadVctReissnerMindlinPageWise(testCase)
end

end