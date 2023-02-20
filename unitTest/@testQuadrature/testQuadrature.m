classdef testQuadrature < matlab.unittest.TestCase
%% Class definition
%
% Test suites for the Gauss quadrature rule
% - Tests the quadrature rule over a unit domain
% - Tests the quadrature rule over a unit triangle

properties
end

%% Method definitions
methods (Test)
    testQuadratureOverUnitDomain(testCase)
    testQuadratureOverTriangle(testCase)
end

end