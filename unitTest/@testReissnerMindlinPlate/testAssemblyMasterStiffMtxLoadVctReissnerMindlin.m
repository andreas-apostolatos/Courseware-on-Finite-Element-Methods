%[text] %[text:anchor:T_FC84E213] # Test the Assembly to the Master Stiffness Matrix and Load Vector for the Finite Element Formulation of the Reissner-Mindin Plate Problem
%[text:tableOfContents]{"heading":"**Table of Contents**"}
%[text] %[text:anchor:H_7E81E991] ## Brief summary of this function
%[text] Test the assembly to the master stiffness matrix and load vector for the Finite Element formulation of the Reissner-Mindin plate problem. It is used both the bilinear and the biquadratic two-dimensional quadrilateral elements for this purpose.
function testAssemblyMasterStiffMtxLoadVctReissnerMindlin(testCase)
%[text] %[text:anchor:H_2C2CD2D8] ## **Read-in the expected solutions**
%[text] %[text:anchor:H_C59195BB] ### Read-in the expected solution for the master stiffness matrix when using the two-dimensional bilinear  four-noded Finite Element mesh
    pathToFileKQ1 = fullfile("unitTest" ,"@testReissnerMindlinPlate", "data", "expK_Q1.txt");
    fl = exist(pathToFileKQ1, "file");
    if fl == 2
        expKQ1 = readmatrix(pathToFileKQ1);
    else
        error("%s does not appear to be an existent file", pathToFileKQ1)
    end
    pathToFileFQ1 = fullfile("unitTest" ,"@testReissnerMindlinPlate", "data", "expF_body_Q1.txt");
    fl = exist(pathToFileFQ1, "file");
    if fl == 2
        expFBodyQ1 = readmatrix(pathToFileFQ1);
    else
        error("%s does not appear to be an existent file", pathToFileFQ1)
    end
%[text] %[text:anchor:H_2B2C9A6D] ### Read-in the expected solution for the master stiffness matrix when using the two-dimensional biquadratic nine-noded Finite Element mesh
    pathToFileKQ2 = fullfile("unitTest" ,"@testReissnerMindlinPlate", "data", "expK_Q2.txt");
    fl = exist(pathToFileKQ2, "file");
    if fl == 2
        expKQ2 = readmatrix(pathToFileKQ2);
    else
        error("%s does not appear to be an existent file", pathToFileKQ2)
    end
    pathToFileFQ2 = fullfile("unitTest" ,"@testReissnerMindlinPlate", "data", "expF_body_Q2.txt");
    fl = exist(pathToFileFQ2, "file");
    if fl == 2
        expFBodyQ2 = readmatrix(pathToFileFQ2);
    else
        error("%s does not appear to be an existent file", pathToFileFQ2)
    end
%[text] %[text:anchor:H_4BD0AE6F] ## Definition of the tolerances
    scaleTol1 = 1e1;
    scaleTol2 = scaleTol1*1e1;
    tolK_Q1 = eps(norm(expKQ1))*scaleTol2;
    tolF_Q1 = eps(norm(expFBodyQ1))*scaleTol1;
    tolK_Q2 = eps(norm(expKQ2))*scaleTol2;
    tolF_Q2 = eps(norm(expFBodyQ2))*scaleTol1;
%[text] %[text:anchor:H_3DF377AE] ## **Problem setup**
%[text] %[text:anchor:H_754843BD] ### Definition of the geometric parameters
    propStr.t = .01; % Thickness
    X0 = 0;
    XLx = 4;
    Y0 = 0;
    YLy = 1;
%[text] %[text:anchor:H_73F29CAE] ### Definition of the loading and material paramerters
    propStr.pBar = -1e-1; % distributed load [N/m^2]
    propStr.mxBar = 0; % distributed moment around x-axis [N]
    propStr.myBar = 0; % distributed moment around y-axis [N]
    propStr.E = 1e7; % Young's modulus [N/m^2]
    propStr.nu = 0; % 0.3 Poisson ration [dimensionless]
    propStr.G = propStr.E/2/(1 + propStr.nu); % shear modulus (connected to epsilon_12 = E/(1+nu)) [N/m^2]
    propStr.D = propStr.E*propStr.t^3/12/(1 - propStr.nu^2); % Plate stiffness [Nm]
    propStr.alpha = 5/6; % Shear correction factor [dimensionless]
%[text] %[text:anchor:H_73B279EC] ## Discretization
%[text] %[text:anchor:H_F2FA9FA6] ### Generate a two-dimensional quadrilateral mesh
%[text] %[text:anchor:H_BCC3B8B6] #### Choose number of elements in the mesh
    numElx = 12; % number of elements along the X-direction
    numEly = 3; % number of elements along the Y-direction
    numEl = numElx*numEly; % total number of elements
%[text] %[text:anchor:H_5A9A09EB] #### Generation of a quadrilateral mesh
    [nodesX, nodesY] = generateQuadrilateralMesh(X0, XLx, Y0, YLy, numElx, numEly); 
%[text] %[text:anchor:H_ED2965ED] ### Generate a two-dimensional quadrilateral mesh using bilinear three-noded elements
%[text] %[text:anchor:H_B970E577] #### Mesh generation
    numNodesElQ1 = 4;
    mshQ1 = generateBilinearQuadrilateralMesh(numElx, numEly, numEl, nodesX, nodesY, numNodesElQ1);
    numNodesQ1 = numel(mshQ1.nodes(:, 1));
%[text] %[text:anchor:H_29DC2F44] #### Rendering the mesh non-uniform
    pathToFileVarvec = fullfile("unitTest" ,"@testReissnerMindlinPlate", "data", "varvec.txt");
    fl = exist(pathToFileVarvec, "file");
    if fl == 2
        varvec = readmatrix(pathToFileVarvec);
    else
        error("%s does not appear to be an existent file", pathToFileVarvec)
    end
    for ii = 1:numNodesQ1
        mshQ1.nodes(ii, 1:2) = mshQ1.nodes(ii, 1:2) + varvec(ii, 1:2);
    end
    nodesX = reshape(mshQ1.nodes(:, 1), height(nodesX), []);
    nodesY = mshQ1.nodes(:, 2);
%[text] %[text:anchor:H_7EC8B704] ### Generate a two-dimensional quadrilateral mesh using biquadratic nine-noded elements
%[text] %[text:anchor:H_219C8C46] #### **Mesh generation**
    numNodesElQ2 = 9;
    [mshQ2, ~] = generateBiquadraticQuadrilateralMesh ...
        (numElx, numEly, nodesX, nodesY, mshQ1, numNodesElQ1, numNodesElQ2);
%[text] %[text:anchor:H_831A355A] ## Computation of the master stiffness matrices
%[text] %[text:anchor:H_E0BA4EA3] ### Computation of the master stiffness matrix using two-dimensional bilinear four-noded elements
    computeStiffMatrixandForceVectorQ1 = @getPrecomputedElStiffMatrixandForceVectorReissnerMindlinPlate;
    computeBasisFunctionsAndDerivsQ1 = @computeBilinearBasisFunctionsAndFirstDerivatives;
    [KQ1, FBodyQ1] = computeMasterStiffMatrixandForceVectorReissnerMindlinPlate ...
        (mshQ1, computeBasisFunctionsAndDerivsQ1, computeStiffMatrixandForceVectorQ1, propStr);
%[text] %[text:anchor:H_7F0B8BAE] ### Computation of the master stiffness matrix using two-dimensional bilinear four-noded elements
    computeStiffMatrixandForceVectorQ2 = @getPrecomputedElStiffMatrixandForceVectorReissnerMindlinPlate;
    computeBasisFunctionsAndDerivsQ2 = @computeBiquadraticBasisFunctionsAndFirstDerivatives;
    [KQ2, FBodyQ2] = computeMasterStiffMatrixandForceVectorReissnerMindlinPlate ...
        (mshQ2, computeBasisFunctionsAndDerivsQ2, computeStiffMatrixandForceVectorQ2, propStr);
%[text] %[text:anchor:H_45A178E9] ### Verification of the master stiffness matrices and load vectors
%[text] %[text:anchor:H_DB1DECDF] ### Verification of the master stiffness matrix and load vector of the two-dimensional bilinear four-noded element
    testCase.verifyEqual(KQ1, expKQ1, "AbsTol", tolK_Q1);
    testCase.verifyEqual(FBodyQ1, expFBodyQ1, "AbsTol", tolF_Q1);
%[text] %[text:anchor:H_C4B35170] ### Verification of the master stiffness matrix and load vector of the one-dimensional quadratic element using full integration
    testCase.verifyEqual(KQ2, expKQ2, "AbsTol", tolK_Q2);
    testCase.verifyEqual(FBodyQ2, expFBodyQ2, "AbsTol", tolF_Q2);
%[text] 
end
%%
%[text] %[text:anchor:H_2BCFFFB9] ## Read-in precomputed values for the element stiffness matrices and load vectors corresponding to the Finite Element formulation of the Reissner-Mindlin plate
function [Ke, Fe, invHeGe] = getPrecomputedElStiffMatrixandForceVectorReissnerMindlinPlate ...
    (~, computeBasisFunctionsAndDerivs, propStr)
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
if ~isa(computeBasisFunctionsAndDerivs, "function_handle")
    error("Input variable 'computeBasisFunctionsAndDerivs' should be a function handle to the computation of the Finite Element basis functions")
else
    if strcmp(func2str(computeBasisFunctionsAndDerivs), "computeBilinearBasisFunctionsAndFirstDerivatives")
        txtQ = "Q1";
    elseif strcmp(func2str(computeBasisFunctionsAndDerivs), "computeBiquadraticBasisFunctionsAndFirstDerivatives")
        txtQ = "Q2";
    else
        error("Input variable 'computeBasisFunctionsAndDerivs' should be either function handle 'computeBilinearBasisFunctionsAndFirstDerivatives' (computation of bilinear basis functions) or 'computeBiquadraticBasisFunctionsAndFirstDerivatives' (computation of biquadratic basis functions) ")
    end
end
%[text] %[text:anchor:H_31F989D0] ### Reading-in the precomputed element stiffness matrix
filename = ['Ke' num2str(idEl) '_' num2str(txtQ) '.txt'];
path_to_file = fullfile("unitTest" ,"@testReissnerMindlinPlate", "data", filename);
fl = exist(path_to_file, "file");
if fl == 2
    Ke = readmatrix(path_to_file);
else
    error("%s does not appear to be an existent file", path_to_file)
end
%[text] %[text:anchor:H_9A307628] ### Reading-in the precomputed element load vector
filename = ['Fe' num2str(idEl) '_' num2str(txtQ) '.txt'];
path_to_file = fullfile("unitTest" ,"@testReissnerMindlinPlate", "data", filename);
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
