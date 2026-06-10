%[text] %[text:anchor:T_FC84E213] # Test the Finite Element Formulation of the Reissner Mindlin Plate Formulation
%[text:tableOfContents]{"heading":"**Table of Contents**"}
%[text] %[text:anchor:H_7E81E991] ## Brief summary of this function
%[text] Tests the constant, the linear and the quadratic one-dimensional basis functions which are defined in the unit space $\\xi \\in \\left\[ -1, 1 \\right\]$
function testComputationElStiffMtxLoadVctReissnerMindlin(testCase)
%[text] %[text:anchor:H_2C2CD2D8] ## **Definition of the expected solutions**
%[text] %[text:anchor:H_C59195BB] ### Definition of the expected solution for the stiffness matrix and load vector of the two-dimensional bilinear four-noded element
    path_to_file_Ke_Q1 = fullfile("unitTest" ,"@testReissnerMindlinPlate", "data", "Ke_Q1.txt");
    fl = exist(path_to_file_Ke_Q1, "file");
    if fl == 2
        expKe_Q1 = readmatrix(path_to_file_Ke_Q1);
    else
        error("%s does not appear to be an existent file", path_to_file_Ke_Q1)
    end
    path_to_file_Fe_Q1 = fullfile("unitTest" ,"@testReissnerMindlinPlate", "data", "Fe_Q1.txt");
    fl = exist(path_to_file_Fe_Q1, "file");
    if fl == 2
        expFe_Q1 = readmatrix(path_to_file_Fe_Q1);
    else
        error("%s does not appear to be an existent file", path_to_file_Ke_Q1)
    end
%[text] %[text:anchor:H_2B2C9A6D] ### Definition of the expected solution for the stiffness matrix and load vector of the two-dimensional biquadratic nine-noded element
    path_to_file_Ke_Q2 = fullfile("unitTest" ,"@testReissnerMindlinPlate", "data", "Ke_Q2.txt");
    fl = exist(path_to_file_Ke_Q2, "file");
    if fl == 2
        expKe_Q2 = readmatrix(path_to_file_Ke_Q2);
    else
        error("%s does not appear to be an existent file", path_to_file_Ke_Q2)
    end
    path_to_file_Fe_Q2 = fullfile("unitTest" ,"@testReissnerMindlinPlate", "data", "Fe_Q2.txt");
    fl = exist(path_to_file_Fe_Q2, "file");
    if fl == 2
        expFe_Q2 = readmatrix(path_to_file_Fe_Q2);
    else
        error("%s does not appear to be an existent file", path_to_file_Fe_Q2)
    end
%[text] %[text:anchor:H_B0452F65] ### Definition of the expected solution for the stiffness matrix and load vector of the two-dimensional bilinear nine-noded element using Assumed Natural Strain (ANS)
    path_to_file_Ke_Q1_ANS = fullfile("unitTest" ,"@testReissnerMindlinPlate", "data", "Ke_Q1_ANS.txt");
    fl = exist(path_to_file_Ke_Q1_ANS, "file");
    if fl == 2
        expKe_Q1_ANS = readmatrix(path_to_file_Ke_Q1_ANS);
    else
        error("%s does not appear to be an existent file", path_to_file_Ke_Q1_ANS)
    end
    path_to_file_Fe_Q1_ANS = fullfile("unitTest" ,"@testReissnerMindlinPlate", "data", "Fe_Q1_ANS.txt");
    fl = exist(path_to_file_Fe_Q1_ANS, "file");
    if fl == 2
        expFe_Q1_ANS = readmatrix(path_to_file_Fe_Q1_ANS);
    else
        error("%s does not appear to be an existent file", path_to_file_Fe_Q1_ANS)
    end
%[text] %[text:anchor:H_4BD0AE6F] ## Definition of the tolerances
    tolScale1 = 1e1;
    tolScale2 = 1e1*tolScale1;
    tolKe_Q1 = max(eps(norm(expKe_Q1))*tolScale1, eps);
    tolFe_Q1 = max(eps(norm(expFe_Q1))*tolScale1, eps);
    tolKe_Q2 = max(eps(norm(expKe_Q2))*tolScale1, eps);
    tolFe_Q2 = max(eps(norm(expFe_Q2))*tolScale1, eps);
    tolKe_Q1_ANS = max(eps(norm(expKe_Q1_ANS))*tolScale1, eps)*tolScale2;
    tolFe_Q1_ANS = max(eps(norm(expFe_Q1_ANS))*tolScale1, eps);
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
%[text] %[text:anchor:H_59A41312] ## Discretization
%[text] %[text:anchor:H_F2FA9FA6] ### Generate a two-dimensional quadrilateral mesh
%[text] %[text:anchor:H_BCC3B8B6] #### Choose number of elements in the mesh
    numElx = 12; % number of elements along the X-direction
    numEly = 3; % number of elements along the Y-direction
    numEl = numElx*numEly; % total number of elements
%[text] %[text:anchor:H_5A9A09EB] #### Generation of a quadrilateral mesh
    [nodesX, nodesY] = generateQuadrilateralMesh(X0, XLx, Y0, YLy, numElx, numEly); 
%[text] %[text:anchor:H_39E22DFB] ### Generate a two-dimensional quadrilateral mesh using bilinear three-noded elements
%[text] %[text:anchor:H_B970E577] #### Mesh generation
    numNodesEl_Q1 = 4;
    msh_Q1 = generateBilinearQuadrilateralMesh(numElx, numEly, numEl, nodesX, nodesY, numNodesEl_Q1);
    numNodes_Q1 = numel(msh_Q1.nodes(:, 1));
%[text] %[text:anchor:H_29DC2F44] #### Rendering the mesh non-uniform
    path_to_file_varvec = fullfile("unitTest" ,"@testReissnerMindlinPlate", "data", "varvec.txt");
    fl = exist(path_to_file_varvec, "file");
    if fl == 2
        varvec = readmatrix(path_to_file_varvec);
    else
        error("%s does not appear to be an existent file", path_to_file_varvec)
    end
    for ii = 1:numNodes_Q1
        msh_Q1.nodes(ii, 1:2) = msh_Q1.nodes(ii, 1:2) + varvec(ii, 1:2);
    end
    nodesX = reshape(msh_Q1.nodes(:, 1), height(nodesX), []);
    nodesY = msh_Q1.nodes(:, 2);
%[text] %[text:anchor:H_12DC84DF] ### Generate a two-dimensional quadrilateral mesh using biquadratic nine-noded elements
%[text] %[text:anchor:H_219C8C46] #### **Mesh generation**
    numNodesEl_Q2 = 9;
    [msh_Q2, ~] = generateBiquadraticQuadrilateralMesh ...
        (numElx, numEly, nodesX, nodesY, msh_Q1, numNodesEl_Q1, numNodesEl_Q2);
%[text] %[text:anchor:H_831A355A] ## Computation of the element stiffness matrices and load vectors
%[text] %[text:anchor:H_E0BA4EA3] ### Computation of the element stiffness matrix and load vector using two-dimensional bilinear four-noded elements
    computeBasisFunctionsAndDerivs_Q1 = @computeBilinearBasisFunctionsAndFirstDerivatives;
    el1 = 1;
    id_nodes_el = msh_Q1.elements(el1, :);
    X = msh_Q1.nodes(id_nodes_el, 1:2);
    [Ke_Q1, Fe_Q1] = computeElementStiffMatrixandForceVectorReissnerMindlinPlate ...
        (X, computeBasisFunctionsAndDerivs_Q1, propStr);
%[text] %[text:anchor:H_7346567A] ### Computation of the element stiffness matrix and load vector using two-dimensional biquadratic nine-noded elements
    computeBasisFunctionsAndDerivs_Q2 = @computeBiquadraticBasisFunctionsAndFirstDerivatives;
    el9 = 9;
    id_nodes_el = msh_Q2.elements(el9, :);
    X = msh_Q2.nodes(id_nodes_el, 1:2);
    [Ke_Q2, Fe_Q2] = computeElementStiffMatrixandForceVectorReissnerMindlinPlate ...
        (X, computeBasisFunctionsAndDerivs_Q2, propStr);
%[text] %[text:anchor:H_29C9EACC] ### Computation of the element stiffness matrix and load vector using two-dimensional bilinear four-noded elements with Assumed Natural Strain (ANS)
    computeBasisFunctionsAndDerivs_Q1_ANS = @computeBilinearBasisFunctionsAndFirstDerivatives; %[text:anchor:H_75DF3A39]
    el17 = 17;
    id_nodes_el = msh_Q1.elements(el17, :);
    X = msh_Q1.nodes(id_nodes_el, 1:2);
    [Ke_Q1_ANS, Fe_Q1_ANS] = computeElementStiffMatrixandForceVectorReissnerMindlinPlateANS ...
        (X, computeBasisFunctionsAndDerivs_Q1_ANS, propStr);
%[text] %[text:anchor:H_45A178E9] ### Verification of the stiffness matrices and load vectors
%[text] %[text:anchor:H_DB1DECDF] ### Verification of the stiffness matrix and load vector using two-dimensional bilinear four-noded elements
    testCase.verifyEqual(Ke_Q1, expKe_Q1, "AbsTol", tolKe_Q1);
    testCase.verifyEqual(Fe_Q1, expFe_Q1, "AbsTol", tolFe_Q1);
%[text] %[text:anchor:H_C4B35170] ### Verification of the stiffness matrix and load vector using two-dimensional biquadratic nine-noded elements
    testCase.verifyEqual(Ke_Q2, expKe_Q2, "AbsTol", tolKe_Q2);
    testCase.verifyEqual(Fe_Q2, expFe_Q2, "AbsTol", tolFe_Q2);
%[text] %[text:anchor:H_010A5049] ### Verification of the stiffness matrix and load vector using two-dimensional bilinear four-noded elements and the Assumed Natural Strain (ANS) method
    testCase.verifyEqual(Ke_Q1_ANS, expKe_Q1_ANS, "AbsTol", tolKe_Q1_ANS);
    testCase.verifyEqual(Fe_Q1_ANS, expFe_Q1_ANS, "AbsTol", tolFe_Q1_ANS);
%[text] 
end

%[appendix]{"version":"1.0"}
%---
