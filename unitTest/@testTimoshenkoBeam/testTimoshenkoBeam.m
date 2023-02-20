classdef testTimoshenkoBeam < matlab.unittest.TestCase
%% Class definition
%
% Test suites for the Finite Element formulation of the Timoshenko beam
% problem
% - Tests the computation of the element stiffness matrices
% - Tests the assembly of the element stiffness matrices to the master
%   stiffness matrix

properties
end

%% Method definitions
methods (Test)
    testComputationElStiffMtxLoadVctTimoshenko(testCase)
    testAssemblyMasterStiffMtxLoadVctTimoshenko(testCase)
    testMasterStiffMtxLoadVctTimoshenkoPageWise(testCase)
end

end