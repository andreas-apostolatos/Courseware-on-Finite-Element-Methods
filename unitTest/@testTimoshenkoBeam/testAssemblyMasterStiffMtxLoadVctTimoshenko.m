%[text] %[text:anchor:T_FC84E213] # Test the Assembly to the Master Stiffness Matrix and Load Vector for the Finite Element Formulation of the Timoshenko Beam Problem
%[text:tableOfContents]{"heading":"**Table of Contents**"}
%[text] %[text:anchor:H_7E81E991] ## Brief summary of this function
%[text] Test the assembly to the master stiffness matrix and consistent nodal force vector for the Finite Element formulation of the Timoshenko beam problem. It is used both the linear and the quadratic one-dimensional elements with full-integration for this purpose.
function testAssemblyMasterStiffMtxLoadVctTimoshenko(testCase)
%[text] %[text:anchor:H_2C2CD2D8] ## **Reading-in the expected solutions**
%[text] %[text:anchor:H_C59195BB] ### Reading-in the expected solution for the master stiffness matrix when using the linear one-dimensional Finite Element mesh
    path_to_file_K_L1 = fullfile("unitTest" ,"@testTimoshenkoBeam", "data", "expK_L1.txt");
    fl = exist(path_to_file_K_L1, "file");
    if fl == 2
        expK_L1 = readmatrix(path_to_file_K_L1);
    else
        error("%s does not appear to be an existent file", path_to_file_K_L1)
    end
    path_to_file_F_L1 = fullfile("unitTest" ,"@testTimoshenkoBeam", "data", "expF_body_L1.txt");
    fl = exist(path_to_file_F_L1, "file");
    if fl == 2
        expF_body_L1 = readmatrix(path_to_file_F_L1);
    else
        error("%s does not appear to be an existent file", path_to_file_F_L1)
    end
%[text] %[text:anchor:H_2B2C9A6D] ### Reading-in the expected solution for the master stiffness matrix when using the quadratic one-dimensional Finite Element mesh
    path_to_file_K_L2 = fullfile("unitTest" ,"@testTimoshenkoBeam", "data", "expK_L2.txt");
    fl = exist(path_to_file_K_L2, "file");
    if fl == 2
        expK_L2 = readmatrix(path_to_file_K_L2);
    else
        error("%s does not appear to be an existent file", path_to_file_K_L2)
    end
    path_to_file_F_L2 = fullfile("unitTest" ,"@testTimoshenkoBeam", "data", "expF_body_L2.txt");
    fl = exist(path_to_file_F_L2, "file");
    if fl == 2
        expF_body_L2 = readmatrix(path_to_file_F_L2);
    else
        error("%s does not appear to be an existent file", path_to_file_F_L2)
    end
%[text] %[text:anchor:H_4BD0AE6F] ## Definition of the tolerances
    scaleTol = 1e1;
    tolK_L1 = eps(norm(expK_L1))*scaleTol;
    tolF_L1 = eps(norm(expF_body_L1));
    tolK_L2 = eps(norm(expK_L2))*scaleTol;
    tolF_L2 = eps(norm(expF_body_L2));
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
    msh_L1 = generateLinearMeshOnLine(X0, XL, numEl);
    varvec_L1 = [-0.478227936430243, -0.254155158477091, -0.226527122592647, -0.094791550357741];
    counterL1 = 1;
    for ii = 2:numel(msh_L1.nodes(:, 1)) - 1
        msh_L1.nodes(ii, 1) = msh_L1.nodes(ii, 1) + varvec_L1(counterL1);
        counterL1 = counterL1 + 1;
    end
%[text] %[text:anchor:H_12DC84DF] ### Generation of an one-dimensional mesh using quadratic three-noded elements
    msh_L2 = generateQuadraticMeshOnLine(msh_L1);
    varvec_L2 = [0.004911427913199, -0.259052626818722, -0.148448603313542, 0.188134139231087, -0.294237326523659];
    counter_L2 = 1;
    for ii = numel(msh_L1.nodes(:, 1)) + 1:numel(msh_L2.nodes(:, 1))
        msh_L2.nodes(ii, 1) = msh_L2.nodes(ii, 1) + varvec_L2(counter_L2);
        counter_L2 = counter_L2 + 1;
    end
%[text] %[text:anchor:H_831A355A] ## Computation of the master stiffness matrices
%[text] %[text:anchor:H_E0BA4EA3] ### Computation of the master stiffness matrix using one-dimensional linear elements and full-integration
    computeStiffMatrixandForceVectorTB = @getPrecomputedElementStiffMatrixandForceVectorTimoshenkoBeam;
    computeBasisFunctionsAndDerivs = @computeLinearBasisFunctionsAndFirstDerivatives;
    [K_L1, F_body_L1] = computeMasterStiffMatrixandForceVectorTimoshenkoBeam ...
        (msh_L1, computeStiffMatrixandForceVectorTB, computeBasisFunctionsAndDerivs, 'undefined', propStr);
%[text] %[text:anchor:H_7F0B8BAE] ### Computation of the master stiffness matrix using one-dimensional quadratic elements and full-integration
    computeStiffMatrixandForceVectorTB = @getPrecomputedElementStiffMatrixandForceVectorTimoshenkoBeam;
    computeBasisFunctionsAndDerivs = @computeQuadraticBasisFunctionsAndFirstDerivatives;
    [K_L2, F_body_L2] = computeMasterStiffMatrixandForceVectorTimoshenkoBeam ...
        (msh_L2, computeStiffMatrixandForceVectorTB, computeBasisFunctionsAndDerivs, 'undefined', propStr);
%[text] %[text:anchor:H_45A178E9] ### Verification of the solutions
%[text] %[text:anchor:H_DB1DECDF] ### Verification of the stiffness matrix and load vector of the one-dimensional linear element using full integration
    testCase.verifyEqual(K_L1, expK_L1, "AbsTol", tolK_L1);
    testCase.verifyEqual(F_body_L1, expF_body_L1, "AbsTol", tolF_L1);
%[text] %[text:anchor:H_C4B35170] ### Verification of the stiffness matrix and load vector of the one-dimensional quadratic element using full integration
    testCase.verifyEqual(K_L2, expK_L2, "AbsTol", tolK_L2);
    testCase.verifyEqual(F_body_L2, expF_body_L2, "AbsTol", tolF_L2);
%[text] 
end
%%
%[text] %[text:anchor:H_2BCFFFB9] ## Read-in precomputed values for the element stiffness matrices regarding the one-dimensional linear mesh
function [Ke, Fe, invHeGe] = getPrecomputedElementStiffMatrixandForceVectorTimoshenkoBeam ...
    (~, computeBasisFunctionsAndDerivs_u, ~, propStr)
%[text] %[text:anchor:H_D37AF8C0] ### Check input
if ~isstruct(propStr)
    error("Input variable 'propStr' should be a MATLAB-struct containing the element index")
else
    if ~isfield(propStr, 'idEl')
        error("Input struct-variable 'propStruct' should contain the element index in field 'idEl'")
    end
end
invHeGe = 'undefined';
%[text] %[text:anchor:H_79CFB6CB] ### Find the element index
idEl = propStr.idEl;
%[text] %[text:anchor:H_F06FF521] ### Verification whether the input is a linear or a quadratic one-dimensional mesh
if ~isa(computeBasisFunctionsAndDerivs_u, "function_handle")
    error("Input variable 'computeBasisFunctionsAndDerivs_u' should be a function handle to the computation of the Finite Element basis functions")
else
    if strcmp(func2str(computeBasisFunctionsAndDerivs_u), "computeLinearBasisFunctionsAndFirstDerivatives")
        txtL = "L1";
    elseif strcmp(func2str(computeBasisFunctionsAndDerivs_u), "computeQuadraticBasisFunctionsAndFirstDerivatives")
        txtL = "L2";
    else
        error("Input variable 'computeBasisFunctionsAndDerivs_u' should be either function handle 'computeLinearBasisFunctionsAndFirstDerivatives' (computation of linear basis functions) or 'computeQuadraticBasisFunctionsAndFirstDerivatives' (computation of quadratic basis functions) ")
    end
end
%[text] %[text:anchor:H_31F989D0] ### Reading-in the precomputed element stiffness matrix
filename = ['Ke' num2str(idEl) '_' num2str(txtL) '.txt'];
path_to_file = fullfile("unitTest" ,"@testTimoshenkoBeam", "data", filename);
fl = exist(path_to_file, "file");
if fl == 2
    Ke = readmatrix(path_to_file);
else
    error("%s does not appear to be an existent file", path_to_file)
end
%[text] %[text:anchor:H_9A307628] ### Reading-in the precomputed element load vector
filename = ['Fe' num2str(idEl) '_' num2str(txtL) '.txt'];
path_to_file = fullfile("unitTest" ,"@testTimoshenkoBeam", "data", filename);
fl = exist(path_to_file, "file");
if fl == 2
    Fe = readmatrix(path_to_file);
else
    error("%s does not appear to be an existent file", path_to_file)
end
%[text] 
end

%[appendix]{"version":"1.0"}
%---
