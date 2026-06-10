%[text] %[text:anchor:T_0F5C4D7A] # Test the Page-Wise Implementation of the Finite Element Formulation for the Reissner-Mindlin Plate Problem
%[text:tableOfContents]{"heading":"**Table of Contents**"}
%[text] %[text:anchor:H_23F418DD] ## Brief summary of this function
%[text] Test the computation of the master stiffness matrix and consistent nodal force vector for the Finite Element formulation of the Reissner-Mindlin plate problem using the [page-wise](https://blogs.mathworks.com/loren/2021/01/14/paged-matrix-functions/) computation of the element stiffness matrices and the [`sparse`](https://www.mathworks.com/help/matlab/math/computational-advantages-of-sparse-matrices.html)-matrix construction of the master stiffness matrix.
function testMasterStiffMtxLoadVctReissnerMindlinPageWise(testCase)
%[text] %[text:anchor:H_2C2CD2D8] ## **Read-in the expected solutions**
%[text] %[text:anchor:H_C59195BB] ### Read-in the expected solution for the master stiffness matrix when using the two-dimensional bilinear  four-noded Finite Element mesh and the page-wise implementation
    pathToFileKQ1 = fullfile("unitTest" ,"@testReissnerMindlinPlate", "data", "expKQ1PageWise.txt");
    fl = exist(pathToFileKQ1, "file");
    if fl == 2
        expKQ1PageWise = readmatrix(pathToFileKQ1);
    else
        error("%s does not appear to be an existent file", pathToFileKQ1)
    end
    pathToFileFQ1 = fullfile("unitTest" ,"@testReissnerMindlinPlate", "data", "expFQ1PageWise.txt");
    fl = exist(pathToFileFQ1, "file");
    if fl == 2
        expFQ1PageWise = readmatrix(pathToFileFQ1);
    else
        error("%s does not appear to be an existent file", pathToFileFQ1)
    end
%[text] %[text:anchor:H_2B2C9A6D] ### Read-in the expected solution for the master stiffness matrix when using the two-dimensional biquadratic nine-noded Finite Element mesh and the page-wise implementation
    pathToFileKQ2 = fullfile("unitTest" ,"@testReissnerMindlinPlate", "data", "expKQ2PageWise.txt");
    fl = exist(pathToFileKQ2, "file");
    if fl == 2
        expKQ2PageWise = readmatrix(pathToFileKQ2);
    else
        error("%s does not appear to be an existent file", pathToFileKQ2)
    end
    pathToFileFQ2 = fullfile("unitTest" ,"@testReissnerMindlinPlate", "data", "expFQ2PageWise.txt");
    fl = exist(pathToFileFQ2, "file");
    if fl == 2
        expFQ2PageWise = readmatrix(pathToFileFQ2);
    else
        error("%s does not appear to be an existent file", pathToFileFQ2)
    end
%[text] %[text:anchor:H_4BD0AE6F] ## Definition of the tolerances
    scaleTol1 = 1e1;
    scaleTol2 = scaleTol1*1e1;
    tolKQ1 = eps(norm(expKQ1PageWise))*scaleTol2;
    tolFQ1 = eps(norm(expFQ1PageWise))*scaleTol1;
    tolKQ2 = eps(norm(expKQ2PageWise))*scaleTol2;
    tolFQ2 = eps(norm(expFQ2PageWise))*scaleTol1;
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
%[text] %[text:anchor:H_82FF0E1E] ### Computation of the master stiffness matrix using two-dimensional bilinear elements, full-integration and page-wise implementation
computeBasisFunctionsAndDerivsQ1 = @computeBilinearBasisFunctionsAndFirstDerivatives;
computeStiffMatrixandForceVectorQ1 = 'undefined';
[KQ1PageWise, FQ1PageWise] = computeMasterStiffMatrixandForceVctReissnerMindlinPlatePageWise ...
    (mshQ1, computeBasisFunctionsAndDerivsQ1, computeStiffMatrixandForceVectorQ1, propStr);
%[text] %[text:anchor:H_249988C6] ### Computation of the master stiffness matrix using two-dimensional biquadratic elements, full-integration and page-wise implementation
computeBasisFunctionsAndDerivsQ2 = @computeBiquadraticBasisFunctionsAndFirstDerivatives;
computeStiffMatrixandForceVectorQ2 = 'undefined';
[KQ2PageWise, FQ2PageWise] = computeMasterStiffMatrixandForceVctReissnerMindlinPlatePageWise ...
    (mshQ2, computeBasisFunctionsAndDerivsQ2, computeStiffMatrixandForceVectorQ2, propStr);
%[text] %[text:anchor:H_45A178E9] ### Verification of the master stiffness matrices and load vectors
%[text] %[text:anchor:H_DB1DECDF] ### Verification of the master stiffness matrix and load vector of the two-dimensional bilinear four-noded element using the page-wise implementation
    testCase.verifyEqual(KQ1PageWise, sparse(expKQ1PageWise), "AbsTol", tolKQ1);
    testCase.verifyEqual(FQ1PageWise, sparse(expFQ1PageWise), "AbsTol", tolFQ1);
%[text] %[text:anchor:H_C4B35170] ### Verification of the master stiffness matrix and load vector of the one-dimensional quadratic element using full integration using the page-wise implementation
    testCase.verifyEqual(KQ2PageWise, sparse(expKQ2PageWise), "AbsTol", tolKQ2);
    testCase.verifyEqual(FQ2PageWise, sparse(expFQ2PageWise), "AbsTol", tolFQ2);
%[text] 
end

%[appendix]{"version":"1.0"}
%---
