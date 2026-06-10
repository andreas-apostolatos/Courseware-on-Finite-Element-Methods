%[text] %[text:anchor:T_BCA49105] # Test the Page-Wise Implementation of the Finite Element Formulation for the Timoshenko Beam Problem
%[text:tableOfContents]{"heading":"**Table of Contents**"}
%[text] %[text:anchor:H_23F418DD] ## Brief summary of this function
%[text] Test the computation of the master stiffness matrix and consistent nodal force vector for the Finite Element formulation of the Timoshenko beam problem using the [page-wise](https://blogs.mathworks.com/loren/2021/01/14/paged-matrix-functions/) computation of the element stiffness matrices and the [`sparse`](https://www.mathworks.com/help/matlab/math/computational-advantages-of-sparse-matrices.html)-matrix construction of the master stiffness matrix.
function testMasterStiffMtxLoadVctTimoshenkoPageWise(testCase)
%[text] %[text:anchor:H_2C2CD2D8] ## **Reading-in the expected solutions**
%[text] %[text:anchor:H_C59195BB] ### Reading-in the expected solution for the master stiffness matrix when using the linear one-dimensional Finite Element mesh and the page-wise implementation
    path_to_file_K_L1 = fullfile("unitTest" ,"@testTimoshenkoBeam", "data", "expKL1PageWise.txt");
    fl = exist(path_to_file_K_L1, "file");
    if fl == 2
        expKL1PageWise = readmatrix(path_to_file_K_L1);
    else
        error("%s does not appear to be an existent file", path_to_file_K_L1)
    end
    path_to_file_F_L1 = fullfile("unitTest" ,"@testTimoshenkoBeam", "data", "expFBodyL1PageWise.txt");
    fl = exist(path_to_file_F_L1, "file");
    if fl == 2
        expFBodyL1PageWise = readmatrix(path_to_file_F_L1);
    else
        error("%s does not appear to be an existent file", path_to_file_F_L1)
    end
%[text] %[text:anchor:H_2B2C9A6D] ### Reading-in the expected solution for the master stiffness matrix when using the quadratic one-dimensional Finite Element mesh and the page-wise implementation
    path_to_file_K_L2 = fullfile("unitTest" ,"@testTimoshenkoBeam", "data", "expKL2PageWise.txt");
    fl = exist(path_to_file_K_L2, "file");
    if fl == 2
        expKL2PageWise = readmatrix(path_to_file_K_L2);
    else
        error("%s does not appear to be an existent file", path_to_file_K_L2)
    end
    path_to_file_F_L2 = fullfile("unitTest" ,"@testTimoshenkoBeam", "data", "expFBodyL2PageWise.txt");
    fl = exist(path_to_file_F_L2, "file");
    if fl == 2
        expFBodyL2PageWise = readmatrix(path_to_file_F_L2);
    else
        error("%s does not appear to be an existent file", path_to_file_F_L2)
    end
%[text] %[text:anchor:H_4BD0AE6F] ## Definition of the tolerances
    scaleTol = 1e1;
    tolKL1PageWise = eps(norm(expKL1PageWise))*scaleTol;
    tolFBodyL1PageWise = eps(norm(expFBodyL1PageWise))*scaleTol;
    tolKL2PageWise = eps(norm(expKL2PageWise))*scaleTol;
    tolFBodyL2PageWise = eps(norm(expFBodyL2PageWise))*scaleTol;
%[text] %[text:anchor:H_3DF377AE] ## **Problem setup**
%[text] %[text:anchor:H_754843BD] ### Definition of the geometric parameters
    propStr.b = .1; % width [m]
    propStr.h = .1; % height [m]
    propStr.A = propStr.b*propStr.h; % height [m]
    X0 = 0; % X-coordinate at the left end of the beam [m]
    XL = 4; % Y-coordinate at the right end of the beam [m]
%[text] %[text:anchor:H_73F29CAE] ### Definition of the loading and material paramerters
    propStr.qBar = -1e2; % distributed load along the beam's length [N/m]
    propStr.mBar = 0; % distributed moment along the beam's length [N]
    propStr.E = 1e7; % Young's modulus [N/m^2]
    propStr.nu = .3; % Poisson ratio [dimensionless]
    propStr.I = propStr.b*(propStr.h^3)/12; % moment of inertia [m^4]
    propStr.G = propStr.E/(2*(1 + propStr.nu)); % shear modulus (connected to epsilon_12 = E/(1+nu)) [N/m^2]
    propStr.alpha = 5/6; % shear correction factor [dimensionless]
%[text] %[text:anchor:H_59A41312] ## Discretization
%[text] %[text:anchor:H_39E22DFB] ### Generation of an one-dimensional mesh using linear two-noded elements
    numEl = 5;
    mshL1 = generateLinearMeshOnLine(X0, XL, numEl);
    varvec_L1 = [-0.478227936430243, -0.254155158477091, -0.226527122592647, -0.094791550357741];
    counterL1 = 1;
    for ii = 2:numel(mshL1.nodes(:, 1)) - 1
        mshL1.nodes(ii, 1) = mshL1.nodes(ii, 1) + varvec_L1(counterL1);
        counterL1 = counterL1 + 1;
    end
%[text] %[text:anchor:H_12DC84DF] ### Generation of an one-dimensional mesh using quadratic three-noded elements
    mshL2 = generateQuadraticMeshOnLine(mshL1);
    varvec_L2 = [0.004911427913199, -0.259052626818722, -0.148448603313542, 0.188134139231087, -0.294237326523659];
    counter_L2 = 1;
    for ii = numel(mshL1.nodes(:, 1)) + 1:numel(mshL2.nodes(:, 1))
        mshL2.nodes(ii, 1) = mshL2.nodes(ii, 1) + varvec_L2(counter_L2);
        counter_L2 = counter_L2 + 1;
    end
%[text] %[text:anchor:H_831A355A] ## Computation of the master stiffness matrices using the page-wise implementation
%[text] %[text:anchor:H_96F96D1B] ### Computation of the master stiffness matrix using one-dimensional linear elements, full-integration and page-wise implementation
    computeBasisFunctionsAndDerivs = @computeLinearBasisFunctionsAndFirstDerivatives;
    computeStiffMatrixandForceVectorTB = 'undefined';
    [KL1PageWise, FBodyL1PageWise] = computeMasterStiffMatrixandForceVectorTimoshenkoBeamPageWise ...
        (mshL1, computeStiffMatrixandForceVectorTB, computeBasisFunctionsAndDerivs, 'undefined', propStr);
%[text] %[text:anchor:H_AAA87A42] ### Computation of the master stiffness matrix using one-dimensional quadratic elements, full-integration and page-wise implementation
    computeBasisFunctionsAndDerivs = @computeQuadraticBasisFunctionsAndFirstDerivatives;
    computeStiffMatrixandForceVectorTB = 'undefined';
    [KL2PageWise, FBodyL2PageWise] = computeMasterStiffMatrixandForceVectorTimoshenkoBeamPageWise ...
        (mshL2, computeStiffMatrixandForceVectorTB, computeBasisFunctionsAndDerivs, 'undefined', propStr);
%[text] %[text:anchor:H_45A178E9] ### Verification of the solutions
%[text] %[text:anchor:H_DB1DECDF] ### Verification of the stiffness matrix and load vector of the one-dimensional linear element using full integration and the page-wise implementation
    testCase.verifyEqual(KL1PageWise, sparse(expKL1PageWise), "AbsTol", tolKL1PageWise);
    testCase.verifyEqual(FBodyL1PageWise, sparse(expFBodyL1PageWise), "AbsTol", tolFBodyL1PageWise);
%[text] %[text:anchor:H_C4B35170] ### Verification of the stiffness matrix and load vector of the one-dimensional quadratic element using full integration and the page-wise implementation
    testCase.verifyEqual(KL2PageWise, sparse(expKL2PageWise), "AbsTol", tolKL2PageWise);
    testCase.verifyEqual(FBodyL2PageWise, sparse(expFBodyL2PageWise), "AbsTol", tolFBodyL2PageWise);
%[text] 
end

%[appendix]{"version":"1.0"}
%---
