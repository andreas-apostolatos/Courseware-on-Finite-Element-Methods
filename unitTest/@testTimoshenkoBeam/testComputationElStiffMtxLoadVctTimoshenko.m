%[text] %[text:anchor:T_FC84E213] # Test the Finite Element Formulation of the Timoshenko Beam Formulation
%[text:tableOfContents]{"heading":"**Table of Contents**"}
%[text] %[text:anchor:H_7E81E991] ## Brief summary of this function
%[text] Tests the computation of the element stiffness matrices for the one-dimensional linear and quadratic elements both with full and selective-reduced integration
function testComputationElStiffMtxLoadVctTimoshenko(testCase)
%[text] %[text:anchor:H_2C2CD2D8] ## **Definition of the expected solutions**
%[text] %[text:anchor:H_C59195BB] ### Definition of the expected solution for the stiffness matrix and load vector of the one-dimensional linear element using full integration
    expKe_L1 = 1.0e+04 * [     9.960865370257244  -1.602564102564103  -9.960865370257244  -1.602564102564103
                              -1.602564102564103   0.369671822152580   1.602564102564103   0.145988536132287
                              -9.960865370257244   1.602564102564103   9.960865370257244   1.602564102564103
                              -1.602564102564103   0.145988536132287   1.602564102564103   0.369671822152580];
    expFe_L1 = [ -16.088603178487855
                                   0
                 -16.088603178487855
                                   0];
%[text] %[text:anchor:H_0EF1A380] ### Definition of the expected solution for the stiffness matrix and load vector of the one-dimensional linear element using selective-reduced integration
    expKe_L1_RI = 1.0e+04 * [  3.129785572012178  -1.602564102564103  -3.129785572012178  -1.602564102564103
                              -1.602564102564103   0.828708578667642   1.602564102564103   0.812433693693179
                              -3.129785572012178   1.602564102564103   3.129785572012178   1.602564102564103
                              -1.602564102564103   0.812433693693179   1.602564102564103   0.828708578667642];
    expFe_L1_RI = [  -51.203638897657598
                                       0
                     -51.203638897657598
                                       0];
%[text] %[text:anchor:H_2B2C9A6D] ### Definition of the expected solution for the stiffness matrix and load vector of the one-dimensional quadratic element using full integration
    expKe_L2 = 1.0e+05 * [     1.766529828637947  -0.160256410256410   0.263690123285247   0.053418803418803  -2.030219951923193  -0.213675213675214
                              -0.160256410256410   0.020929824817015  -0.053418803418803  -0.008156585550019   0.213675213675214  -0.000282126947003
                               0.263690123285247  -0.053418803418803   0.694492297104473   0.160256410256410  -0.958182420389720   0.213675213675214
                               0.053418803418803  -0.008156585550019   0.160256410256410   0.056206271674398  -0.213675213675214   0.027881000261231
                              -2.030219951923193   0.213675213675214  -0.958182420389720  -0.213675213675214   2.988402372312914  -0.000000000000000
                              -0.213675213675214  -0.000282126947003   0.213675213675214   0.027881000261231                   0   0.149244724096978];
    expFe_L2 = [  -3.897227043837921
                                   0
                 -23.690374152310202
                                   0
                 -55.175202392296242
                                   0];
%[text] %[text:anchor:H_9EBCB9DA] ### Definition of the expected solution for the stiffness matrix and load vector of the one-dimensional quadratic element using selective-reduced integration
    expKe_L2_RI = 1.0e+05 * [  0.552304456985113  -0.160256410256410   0.146526566695298   0.053418803418803  -0.698831023680410  -0.213675213675214
                              -0.160256410256410   0.048017325302445  -0.053418803418803  -0.016209764050433   0.213675213675214   0.058164240519316
                               0.146526566695298  -0.053418803418803   1.499067476749051   0.160256410256410  -1.645594043444349   0.213675213675214
                               0.053418803418803  -0.016209764050433   0.160256410256410   0.023679174220628  -0.213675213675214   0.002103186801523
                              -0.698831023680410   0.213675213675214  -1.645594043444349  -0.213675213675214   2.344425067124760   0.000000000000000
                              -0.213675213675214   0.058164240519316   0.213675213675214   0.002103186801523                   0   0.138821370165252];
    expFe_L2_RI = [  -28.071202152654266
                                       0
                      -2.986650255175949
                                       0
                     -62.115704815660408
                                       0];
%[text] %[text:anchor:H_D85EFAC7] ### Definition of the expected solution for the stiffness matrix and load vector of the one-dimensional linear element Multipliers using the Hellinger-Reissner principle and constant discretization for the Lagrange Multipliers
    expKe_L10_HR = 1.0e+04 * [ 3.581983092986053  -1.602564102564103  -3.581983092986053  -1.602564102564103
                              -1.602564102564103   0.726293564982261   1.602564102564103   0.707667252898734
                              -3.581983092986053   1.602564102564103   3.581983092986053   1.602564102564103
                              -1.602564102564103   0.707667252898734   1.602564102564103   0.726293564982261];
    expFe_L10_HR = [    -44.739577517887042
                                          0
                        -44.739577517887042
                                          0];
%[text] %[text:anchor:H_811442B3] ### Definition of the expected solution for the stiffness matrix and load vector of the one-dimensional quadratic element Multipliers using the Hellinger-Reissner principle and linear discretization for the Lagrange Multipliers
    expKe_L21_HR = 1.0e+05 * [    1.537985749447155  -0.126730535608064   0.155826480847230   0.086944678067150  -1.693812230294385  -0.280726962971906
                                  -0.126730535608064   0.015417590567310  -0.037595934726274  -0.013355050664167   0.164326470334338   0.010428572416851
                                   0.155826480847230  -0.037595934726274   0.643584982414059   0.176079278948940  -0.799411463261289   0.182029476290154
                                   0.086944678067150  -0.013355050664167   0.176079278948940   0.051155893012393  -0.263023957016090   0.038129844037383
                                  -1.693812230294385   0.164326470334338  -0.799411463261289  -0.263023957016090   2.493223693555674   0.098697486681752
                                  -0.280726962971906   0.010428572416851   0.182029476290154   0.038129844037383   0.098697486681752   0.128285180956972];
    expFe_L21_HR = [ -3.897227043837921
                                       0
                     -23.690374152310202
                                       0
                     -55.175202392296242
                                       0];
%[text] %[text:anchor:H_4BD0AE6F] ## Definition of the tolerances
    tolKe_L1 = eps(norm(expKe_L1));
    tolFe_L1 = eps(norm(expFe_L1));
    tolKe_L1_RI = eps(norm(expKe_L1_RI));
    tolFe_L1_RI = eps(norm(expFe_L1_RI));
    tolKe_L2 = eps(norm(expKe_L2));
    tolFe_L2 = eps(norm(expFe_L2));
    tolKe_L2_RI = eps(norm(expKe_L2_RI));
    tolFe_L2_RI = eps(norm(expFe_L2_RI));
    tolKe_L10_HR = eps(norm(expKe_L10_HR));
    tolFe_L10_HR = eps(norm(expFe_L10_HR));
    tolKe_L21_HR = eps(norm(expKe_L21_HR));
    tolFe_L21_HR = eps(norm(expFe_L21_HR));
%[text] %[text:anchor:H_3DF377AE] ## **Problem setup**
%[text] %[text:anchor:H_754843BD] ### Definition of the geometric parameters
    propStr.b = .1; % width [m]
    propStr.t = .1; % height [m]
    propStr.A = propStr.b*propStr.t; % height [m]
    X0 = 0; % X-coordinate at the left end of the beam [m]
    XL = 4; % Y-coordinate at the right end of the beam [m]
%[text] %[text:anchor:H_73F29CAE] ### Definition of the loading and material paramerters
    propStr.qBar = -1e2; % distributed load along the beam's length [N/m]
    propStr.mBar = 0; % distributed moment along the beam's length [N]
    propStr.E = 1e7; % Young's modulus [N/m^2]
    propStr.nu = .3; % Poisson ratio [dimensionless]
    propStr.I = propStr.b*(propStr.t^3)/12; % moment of inertia [m^4]
    propStr.G = propStr.E/(2*(1 + propStr.nu)); % shear modulus (connected to epsilon_12 = E/(1+nu)) [N/m^2]
    propStr.alpha = 5/6; % shear correction factor [dimensionless]
%[text] %[text:anchor:H_59A41312] ## Discretization
%[text] %[text:anchor:H_39E22DFB] ### Generate an one-dimensional mesh using linear two-noded elements
    numEl = 5;
    msh_L1 = generateLinearMeshOnLine(X0, XL, numEl);
    varvec_L1 = [-0.478227936430243, -0.254155158477091, -0.226527122592647, -0.094791550357741];
    counterL1 = 1;
    for ii = 2:numel(msh_L1.nodes(:, 1)) - 1
        msh_L1.nodes(ii, 1) = msh_L1.nodes(ii, 1) + varvec_L1(counterL1);
        counterL1 = counterL1 + 1;
    end
%[text] %[text:anchor:H_12DC84DF] ### Generate an one-dimensional mesh using quadratic three-noded elements
    msh_L2 = generateQuadraticMeshOnLine(msh_L1);
    varvec_L2 = [0.004911427913199, -0.259052626818722, -0.148448603313542, 0.188134139231087, -0.294237326523659];
    counter_L2 = 1;
    for ii = numel(msh_L1.nodes(:, 1)) + 1:numel(msh_L2.nodes(:, 1))
        msh_L2.nodes(ii, 1) = msh_L2.nodes(ii, 1) + varvec_L2(counter_L2);
        counter_L2 = counter_L2 + 1;
    end
%[text] %[text:anchor:H_831A355A] ## Computation of the element stiffness matrices
%[text] %[text:anchor:H_E0BA4EA3] ### Computation of the element stiffness matrix using one-dimensional linear elements and full-integration
    computeBasisFunctionsAndDerivs = @computeLinearBasisFunctionsAndFirstDerivatives;
    el1 = 1;
    id_nodes_el = msh_L1.elements(el1, :);
    X_el1 = msh_L1.nodes(id_nodes_el);
    [Ke_L1, Fe_L1, ~] = computeElementStiffMatrixandForceVectorTimoshenkoBeam ...
        (X_el1, computeBasisFunctionsAndDerivs, 'undefined', propStr);
%[text] %[text:anchor:H_7346567A] ### Computation of the element stiffness matrix using one-dimensional linear elements and reduced-integration
    computeBasisFunctionsAndDerivs = @computeLinearBasisFunctionsAndFirstDerivatives;
    el2 = 2;
    id_nodes_el = msh_L1.elements(el2, :);
    X_el2 = msh_L1.nodes(id_nodes_el);
    [Ke_L1_RI, Fe_L1_RI, ~] = computeElementStiffMatrixandForceVectorTimoshenkoBeamRI ...
        (X_el2, computeBasisFunctionsAndDerivs, 'undefined', propStr);
%[text] %[text:anchor:H_CA49FA87] ### Computation of the element stiffness matrix using one-dimensional quadratic elements and full-integration
    computeBasisFunctionsAndDerivs = @computeQuadraticBasisFunctionsAndFirstDerivatives;
    el3 = 3;
    id_nodes_el = msh_L2.elements(el3, :);
    X_el3 = msh_L2.nodes(id_nodes_el);
    [Ke_L2, Fe_L2, ~] = computeElementStiffMatrixandForceVectorTimoshenkoBeam ...
        (X_el3, computeBasisFunctionsAndDerivs, 'undefined', propStr);
%[text] %[text:anchor:H_7F0B8BAE] ### Computation of the element stiffness matrix using one-dimensional quadratic elements and reduced-integration
    computeBasisFunctionsAndDerivs = @computeQuadraticBasisFunctionsAndFirstDerivatives;
    el4 = 4;
    id_nodes_el = msh_L2.elements(el4, :);
    X_el4 = msh_L2.nodes(id_nodes_el);
    [Ke_L2_RI, Fe_L2_RI, ~] = computeElementStiffMatrixandForceVectorTimoshenkoBeamRI ...
        (X_el4, computeBasisFunctionsAndDerivs, 'undefined', propStr);
%[text] %[text:anchor:H_DFDF3E21] ### Computation of the element stiffness matrix using one-dimensional linear elements and the Hellinger-Reissner principle with constant discretization for the Lagrange Multipliers
    computeBasisFunctionsAndDerivs = @computeLinearBasisFunctionsAndFirstDerivatives;
    computeBasisFunctionsAndDerivs_HR_LM = @computeConstantBasisFunctionAndFirstDerivatives;
    el5 = 5;
    id_nodes_el = msh_L1.elements(el5, :);
    X_el5 = msh_L1.nodes(id_nodes_el);
    [Ke_L10_HR, Fe_L10_HR, ~] = computeElStiffMtxForceVctHellingerReissnerTimoshenkoBeam ...
        (X_el5, computeBasisFunctionsAndDerivs, computeBasisFunctionsAndDerivs_HR_LM, propStr);
%[text] %[text:anchor:H_7EBC639C] ### Computation of the element stiffness matrix using two-dimensional linear elements and the Hellinger-Reissner principle with linear discretization for the Lagrange Multipliers
    computeBasisFunctionsAndDerivs = @computeQuadraticBasisFunctionsAndFirstDerivatives;
    computeBasisFunctionsAndDerivs_HR_LM = @computeLinearBasisFunctionsAndFirstDerivatives;
    el3 = 3;
    id_nodes_el = msh_L2.elements(el3, :);
    X_el3 = msh_L2.nodes(id_nodes_el);
    [Ke_L21_HR, Fe_L21_HR, ~] = computeElStiffMtxForceVctHellingerReissnerTimoshenkoBeam ...
        (X_el3, computeBasisFunctionsAndDerivs, computeBasisFunctionsAndDerivs_HR_LM, propStr);
%[text] %[text:anchor:H_45A178E9] ### Verify the solutions
%[text] %[text:anchor:H_DB1DECDF] ### Verify the solution for the stiffness matrix and load vector of the one-dimensional linear element using full integration
    testCase.verifyEqual(Ke_L1, expKe_L1, "AbsTol", tolKe_L1);
    testCase.verifyEqual(Fe_L1, expFe_L1, "AbsTol", tolFe_L1);
%[text] %[text:anchor:H_20361B14] ### Verify the solution for the stiffness matrix and load vector of the one-dimensional linear element using selective-reduced integration
    testCase.verifyEqual(Ke_L1_RI, expKe_L1_RI, "AbsTol", tolKe_L1_RI);
    testCase.verifyEqual(Fe_L1_RI, expFe_L1_RI, "AbsTol", tolFe_L1_RI);
%[text] %[text:anchor:H_C4B35170] ### Verify the solution for the stiffness matrix and load vector of the one-dimensional quadratic element using full integration
    testCase.verifyEqual(Ke_L2, expKe_L2, "AbsTol", tolKe_L2);
    testCase.verifyEqual(Fe_L2, expFe_L2, "AbsTol", tolFe_L2);
%[text] %[text:anchor:H_C802B293] ### Verify the solution for the stiffness matrix and load vector of the one-dimensional quadratic element using selective-reduced integration
    testCase.verifyEqual(Ke_L2_RI, expKe_L2_RI, "AbsTol", tolKe_L2_RI);
    testCase.verifyEqual(Fe_L2_RI, expFe_L2_RI, "AbsTol", tolFe_L2_RI);
%[text] %[text:anchor:H_77A3BA1A] ### Verify the solution for the stiffness matrix and load vector of the one-dimensional linear element using the Hellinger-Reissner principle and constant discretization for the Lagrange Multipliers fields
    testCase.verifyEqual(Ke_L10_HR, expKe_L10_HR, "AbsTol", tolKe_L10_HR);
    testCase.verifyEqual(Fe_L10_HR, expFe_L10_HR, "AbsTol", tolFe_L10_HR);
%[text] %[text:anchor:H_F00334B6] ### Verify the solution for the stiffness matrix and load vector of the one-dimensional quadratic element using the Hellinger-Reissner principle and linear discretization for the Lagrange Multipliers fields
    testCase.verifyEqual(Ke_L21_HR, expKe_L21_HR, "AbsTol", tolKe_L21_HR);
    testCase.verifyEqual(Fe_L21_HR, expFe_L21_HR, "AbsTol", tolFe_L21_HR);
%[text] 
end

%[appendix]{"version":"1.0"}
%---
