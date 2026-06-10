%[text] %[text:anchor:T_7BB7FF35] # Finite Element Methods
%[text:tableOfContents]{"heading":"Table of Contents"}
%[text] %[text:anchor:H_80A2577E] ## Introduction
%[text] 
%[text] This set of interactive lecture notes was developed in the frame of the lecture [***Advanced Finite Element Methods***](https://www.cee.ed.tum.de/en/st/teaching/) (in German: ***Weiterführende Themen der Finite Elemente Methode***). This is a course taught in the second year of the master programs on [***Computational Mechanics***](https://www.tum.de/studium/studienangebot/detail/computational-mechanics-master-of-science-msc) and [***Civil Engineering***](https://www.tum.de/en/studies/degree-programs/detail/civil-engineering-master-of-science-msc).
%[text] The course addresses topics related to the [*Variational (Weak) Formulation*](https://en.wikipedia.org/wiki/Weak_formulation) and [*Finite Element Approximation*](https://en.wikipedia.org/wiki/Finite_element_method) of shear-deformable lightweight structures of the beam and plate types. For a more introductory courseware regarding the structural analysis of beams using the [Symbolic Math Toolbox](https://www.mathworks.com/products/symbolic.html)™ in [MATLAB](https://www.mathworks.com/products/matlab.html)® for the analytical derivation of the underlying equations please refer to the [courseware on Beam Bending and Deflection](https://www.mathworks.com/matlabcentral/fileexchange/113670-beam-bending-and-deflection/) freely available on [File Exchange](https://www.mathworks.com/matlabcentral/fileexchange/).
%[text] Shear-deformable formulations are typically associated with moderately thick to thick beams and plates. These formulations are among other the [*Timoshenko beam*](https://en.wikipedia.org/wiki/Timoshenko%E2%80%93Ehrenfest_beam_theory) and the [*Reissner-Mindlin plate*](https://en.wikipedia.org/wiki/Mindlin%E2%80%93Reissner_plate_theory) theories. The equivalent theories for thin beams and plates are the [*Euler-Bernoulli beam*](https://en.wikipedia.org/wiki/Euler%E2%80%93Bernoulli_beam_theory) and the [*Kirchhoff plate*](https://en.wikipedia.org/wiki/Kirchhoff%E2%80%93Love_plate_theory) theories, respectively. The latter formulations assume negligble shear deformations allowing for explicitly coupling the vertical deflections with the cross-sectional rotations which in turns raises the variational index of the corresponding variational formulations, namely, the highest derivative order that appears in the variational (weak) formulation of the [*Boundary Value Problem (BVP)*](https://en.wikipedia.org/wiki/Boundary_value_problem).
%[text] It is shown that low-order approximations for the Finite Element discretization of such formulations suffer from a numerical artifact the so-called [*Transverse-Shear Locking*](https://de.wikipedia.org/wiki/Locking_(FEM)) (only available in German on [Wikipedia®](https://en.wikipedia.org/wiki/Main_Page)). This numerical artifact results in spurious reaction transverse-shear forces having the underlying results of the Finite Element discretizations to behave stiffer than they should.
%[text] It is then demonstrated how the transverse-shear locking can be alleviated by various methods such as:
%[text] - The use of *selective-reduced integration*
%[text] - The use of *high-order basis functions*
%[text] - The *Assumed Natural Strain* method
%[text] - The [*Hellinger-Reissner Principle*](https://en.wikiversity.org/wiki/Elasticity/Hellinger-Reissner_principle) \
%[text] Moreover, the reader is introduced to the topic of *weak application* of [*Dirichlet boundary conditions*](https://en.wikipedia.org/wiki/Dirichlet_boundary_condition) in a variational formulation. The reader is introduced to the following methods accordingly:
%[text] - [*The Penalty method*](https://en.wikipedia.org/wiki/Penalty_method)
%[text] - [*The Lagrange Multipliers method*](https://en.wikipedia.org/wiki/Lagrange_multiplier) \
%[text] which are well-known from [*mathematical optimization*](https://en.wikipedia.org/wiki/Mathematical_optimization).
%[text] The accompanying files contain a comprehensive set of [*Live Scripts*](https://www.mathworks.com/help/matlab/matlab_prog/create-live-scripts.html) which allow for experimenting with the different Finite Element Methods (FEMs) in an interactive way.
%[text] 
%[text] **Author: Andreas Apostolatos, PhD (aapostol@mathworks.com)**
%[text] **Date: 20.02.2023**
%%
%[text] %[text:anchor:H_503B025E] ## Quick guide
%[text] 
%[text] %[text:anchor:H_955C077E] **Chapter 1**
%[text] [`main_Chapter1_BasisFunctions.m`](file:.\1_BasisFunctions\main_Chapter1_BasisFunctions.m)
%[text] %[text:anchor:H_94021FBA] **Chapter 2**
%[text] [`main_Chapter2_ALinearStraightTimoshenkoBeamElement.m`](file:.\2_TimoshenkoBeam\main_Chapter2_ALinearStraightTimoshenkoBeamElement.m)
%[text] %[text:anchor:H_C058518E] **Chapter 3**
%[text] [`main_Chapter3_AShearDeformablePlateElement.m`](file:.\3_ReissnerMindlinPlate\main_Chapter3_AShearDeformablePlateElement.m)
%[text] %[text:anchor:H_5A306220] **Chapter 4**
%[text] [`main_Chapter4_TransverseShearLocking_TimoshenkoBeam.m`](file:.\4_TransverseShearLocking\4_TransverseShearLocking_TimoshenkoBeam\main_Chapter4_TransverseShearLocking_TimoshenkoBeam.m)
%[text] [`main_Chapter4_Locking_ReissnerMindlinCantileverPlate.m`](file:.\4_TransverseShearLocking\4_TransverseShearLocking_ReissnerMindlinCantileverPlate\main_Chapter4_Locking_ReissnerMindlinCantileverPlate.m)
%[text] [`main_Chapter4_Locking_ReissnerMindlinSquarePlate.m`](file:.\4_TransverseShearLocking\4_TransverseShearLocking_ReissnerMindlinSquarePlate\main_Chapter4_Locking_ReissnerMindlinSquarePlate.m)
%[text] %[text:anchor:H_DC79D498] **Chapter 5**
%[text] [`main_Chapter5_WeakBoundaryConditions_TimoshenkoBeam.m`](file:.\5_WeakDirichletBoundaryConditions\main_Chapter5_WeakBoundaryConditions_TimoshenkoBeam.m)
%[text] **Chapter 6**
%[text] [`main_Chapter6_HellingerReissnerFormulation_TimoshenkoBeam.m`](file:.\6_HellingerReissnerFormulation\main_Chapter6_HellingerReissnerFormulation_TimoshenkoBeam.m)
%[text] [`main_Chapter6_HellingerReissnerFormulation_TimoshenkoBeam_Study.m`](file:.\6_HellingerReissnerFormulation\main_Chapter6_HellingerReissnerFormulation_TimoshenkoBeam_Study.m)
%[text] **Appendix A**
%[text] [`main_AppexA_EfficiencyCondiderationsTimoshenkoBeam.m`](file:.\AppexA_EfficiencyConsiderations\main_AppexA_EfficiencyCondiderationsTimoshenkoBeam.m)
%[text] [`main_AppexA_EfficiencyCondiderationsReissnerMindlinPlate.m`](file:.\AppexA_EfficiencyConsiderations\main_AppexA_EfficiencyCondiderationsReissnerMindlinPlate.m)
%[text] 
%[text] **Unit tests**
%[text] [`main_unitTests.m`](file:.\main_unitTests.m)
%%
%[text] %[text:anchor:H_F298FD6D] ## Chapter 1: Basis Functions
%[text] 
%[text] In this chapter it is introduced the notion of the [*basis functions*](https://en.wikipedia.org/wiki/Basis_function) typically used in the frame of the [*Finite Element Method*](https://en.wikipedia.org/wiki/Finite_element_method). The basis functions introduced herein are the following ones,
%[text] - one-dimensional two-noded linear,
%[text] - one-dimensional three-noded quadratic,
%[text] - two-dimensional four-noded bilinear, and
%[text] - two-dimensional nine-noded biquadratic \
%[text] It is then shown how these basis functions can be constructed using the interpolatory [*Lagrange Polynomials*](https://en.wikipedia.org/wiki/Lagrange_polynomial)*.* Moreover, it is demonstrated how the geometry can be parametrized using the same basis functions as for the field approximation, known as the *Isoparametric Finite Elements*.
%[text] Lastly, it is also demonstrated how the Jacobian of the geometric transformation is affected for different element shapes. This is an important concept as it adds directly to the numerical approximation of the Finite Element formulation that is used.
%[text] To navigate to the corresponding Live Script please select the hyperlink [`main_Chapter1_BasisFunctions.mlx`](file:./1_BasisFunctions/main_Chapter1_BasisFunctions.mlx).
%[text] This Live Script takes about a minute to execute.
%%
%[text] %[text:anchor:H_8541FA9A] ## Chapter 2: Finite Element Formulation the Timoshenko Beam Problem
%[text] 
%[text] In this chapter it is presented the standard variational formulation of the [*Timoshenko beam*](https://en.wikipedia.org/wiki/Timoshenko%E2%80%93Ehrenfest_beam_theory) problem along with its Finite Element discretization. For the discretization it is used the one-dimensional linear two-noded and the quadratic three-noded basis functions. The transverse-shear locking is highlighted for the linear basis, as the corresponding deformations are much stiffer than the ones when using the quadratic basis.
%[text] It is also provided an analytical solution for the given problem setup. This analytical solution is then used to verify the solution obtained by the Finite Element Method.
%[text] To navigate to the corresponding Live Script select the hyperlink [`main_Chapter2_ALinearStraightTimoshenkoBeamElement.mlx`](file:./2_TimoshenkoBeam/main_Chapter2_ALinearStraightTimoshenkoBeamElement.mlx)`.` This Live Script takes about 40 seconds to execute.
%%
%[text] %[text:anchor:H_9A3CC11F] ## Chapter 3: Finite Element Formulation of the Reissner-Mindlin Plate Problem
%[text] 
%[text] In this chapter it is presented the standard variational formulation of the [*Reissner-Mindlin plate*](https://en.wikipedia.org/wiki/Mindlin%E2%80%93Reissner_plate_theory) along with its Finite Element discretization. For the discretization it is used the two-dimensional bilinear four-noded and the biquadratic nine-noded basis functions. The transverse-shear locking is highlighted for the bilinear basis, as the corresponding deformations are much stiffer than the ones when using the biquadratic basis.
%[text] For the selected setup, namely, the cantilever plate subject to distributed vertical load, it can be used the same analytical solution as for the setup of the Timoshenko beam in [Chapter 3](file:./3_ReissnerMindlinPlate\main_Chapter3_AShearDeformablePlateElement.mlx:T_92F22C48). This solution is then used as in the previous chapter to vertify the solution obtained by the Finite Element Method.
%[text] To navigate to the corresponding Live Script select the hyperlink [`main_Chapter3_AShearDeformablePlateElement.mlx`](file:./3_ReissnerMindlinPlate/main_Chapter3_AShearDeformablePlateElement.mlx)`.` The Live Script takes about 6,5 minutes to execute.
%%
%[text] %[text:anchor:H_81518CB1] ## Chapter 4: Transverse Shear-Locking
%[text] 
%[text] In this chapter it is introduced the transverse-shear locking and detailed studies are carried out. The chapter demonstrates the manifestation of the transverse-shear locking and the methods for its alleviation for both the Timoshenko beam and the Reissner-Mindlin plate problems.
%[text] For the Timoshenko beam problem it is used a cantilever beam as shown in [Chapter 2](file:./2_TimoshenkoBeam/main_Chapter2_ALinearStraightTimoshenkoBeamElement.mlx:H_853F8D92). Moreover, it is also used as a reference the solution of the same setup using the [Euler-Bernoulli beam](https://en.wikipedia.org/wiki/Euler%E2%80%93Bernoulli_beam_theory) formulation discretized by means of the [cubic Hermitian polynomials](https://en.wikipedia.org/wiki/Cubic_Hermite_spline). Emphasis is put on the correspondance of the solutions between the two different beam theories based on the beam's slenderness.
%[text] To navigate to the corresponding Live Script select the hyperlink [`main_Chapter4_TransverseShearLocking_TimoshenkoBeam.mlx`](file:./4_TransverseShearLocking/4_TransverseShearLocking_TimoshenkoBeam/main_Chapter4_TransverseShearLocking_TimoshenkoBeam.mlx). The Live Script takes about 40 seconds to execute.
%[text] For the Reissner-Mindlin plate problem it is used both a [cantilever plate](file:./3_ReissnerMindlinPlate/main_Chapter3_AShearDeformablePlateElement.mlx:H_D949C6AD) as shown in Chapter 3 and a rectangular plate fixed along two of its edges while subject to a distributed vertical load.
%[text] - in what concerns the cantilever plate it is used the [analytical solution](file:./3_ReissnerMindlinPlate/main_Chapter3_AShearDeformablePlateElement.mlx:M_74A0D211) presented at Chapter 3 as a reference solution to evaluate the presence of transverse-shear locking,
%[text] - concerning the rectangular plate it is used a reference solution obtained by the [Kirchhoff plate theory](https://en.wikipedia.org/wiki/Kirchhoff%E2%80%93Love_plate_theory). Emphasis is put on the correspondance of the solutions obtained the two different theories depending on slenderness. Moreover, the Finite Element formulation by means of the Assumed Natural Strain is also used that is based on the four-noded bilinear quadrilateral finite elements. This method is proven to be free of transverse-shear locking. \
%[text] To navigate to the corresponding Live Script regarding the cantilever plate select the hyperlink [`main_Chapter4_Locking_ReissnerMindlinCantileverPlate.mlx`](file:./4_TransverseShearLocking/4_TransverseShearLocking_ReissnerMindlinCantileverPlate/main_Chapter4_Locking_ReissnerMindlinCantileverPlate.mlx). The Live Script takes about 1,5 minutes to execute.
%[text] To navigate to the corresponding Live Script regarding the rectangular plate select the hyperlink [`main_Chapter4_Locking_ReissnerMindlinSquarePlate.mlx`](file:./4_TransverseShearLocking/4_TransverseShearLocking_ReissnerMindlinSquarePlate/main_Chapter4_Locking_ReissnerMindlinSquarePlate.mlx). The Live Script takes about a minute to execute.
%[text] For the detection of the transverse-shear locking three studies are carried out:
%[text] - computation and visualization of the resultant shear-force distribution along the beam for the different sets of basis functions
%[text] - slenderness studies
%[text] - refinement studies \
%%
%[text] %[text:anchor:H_8D66F679] ## Chapter 5: Weak Imposition of Dirichlet boundary Conditions for the Timoshenko Beam Problem
%[text] 
%[text] In this chapter it is introduced the concept of the weak application of the Dirichlet boundary conditions for the Timoshenko beam problem. The variational formulation is extended by additional terms that account for the conditions along the Dirichlet boundary without assuming them to be identically zero and forcing the test functions to vanish along the Dirichlet boundary. In this way some error is allowed on the Dirichlet conditions along the boundary, but the solution can easier adapt within the computational domain. Two methods are used for the application of weak Dirichlet boundary conditions accordingly:
%[text] - The [Penalty](https://en.wikipedia.org/wiki/Penalty_method) method, and
%[text] - the [Lagrange Multipliers](https://en.wikipedia.org/wiki/Lagrange_multiplier) method \
%[text] The methods are quantified by computing the error of the vertical deflection and cross-sectional rotation in the $L^2$-norm within the computational domain.
%[text] To navigate to the corresponding Live Script regarding the weak application of Dirichlet boundary conditions select the hyperlink [`main_Chapter5_WeakBoundaryConditions_TimoshenkoBeam.mlx`](file:./5_WeakDirichletBoundaryConditions/main_Chapter5_WeakBoundaryConditions_TimoshenkoBeam.mlx)`.` The Live Script takes about 2 seconds to execute.
%%
%[text] %[text:anchor:H_354875AF] ## Chapter 6: Hellinger-Reissner Formulation for the Timoshenko Beam Problem
%[text] 
%[text] In this chapter it is shown how the standard variational formulation of the Timoshenko beam element can be extended to account for the [Hellinger-Reissner principle](https://en.wikiversity.org/wiki/Elasticity/Hellinger-Reissner_principle) to obtain a transverse-shear locking-free formulation. The same analytical solution is used as in the previous chapters to assess the different Finite Element formulations based on their accuracy. The Live Script allows for different combination of polynomials orders for the primal unknown fields, namely the vertical displacement and cross-sectional roation fields, and the secondary fields, namely, the reaction force and reaction moment fields, in the frame of the Hellinger-Reissner formulation of the Timoshenko beam problem.
%[text] To navigate to the corresponding Live Scripts regarding the Hellinger-Reissner formulation of the Timoshenko beam problem select the hyperlink [`main_Chapter6_HellingerReissnerFormulation_TimoshenkoBeam.mlx`](file:./6_HellingerReissnerFormulation/main_Chapter6_HellingerReissnerFormulation_TimoshenkoBeam.mlx). The Live script takes about 14 seconds to execute.
%[text] Moreover, it is also provided Live Script [`main_Chapter6_HellingerReissnerFormulation_TimoshenkoBeam_Study.mlx`](file:./6_HellingerReissnerFormulation/main_Chapter6_HellingerReissnerFormulation_TimoshenkoBeam_Study.mlx) which performs similar studies as Live Script [`main_Chapter4_TransverseShearLocking_TimoshenkoBeam.mlx`](file:./4_TransverseShearLocking/4_TransverseShearLocking_TimoshenkoBeam/main_Chapter4_TransverseShearLocking_TimoshenkoBeam.mlx) but it also includes insights with respect to the Hellinger-Reissner principle for the Timoshenko beam formulation and how it compares to the rest of the Finite Element formulations including also the Euler-Bernoulli beam formulation. The scripts takes about 42 seconds to execute.
%%
%[text] %[text:anchor:H_079FDE1D] ## Appendix A: Efficiency Considerations
%[text] 
%[text] This chapter provides hints on how the runtime performance of standard finite element routines can be accelerated. In particular the following finite element routines are dealt with:
%[text] - computation of the element stiffness matrices [page-wise](https://blogs.mathworks.com/loren/2021/01/14/paged-matrix-functions/) using vectorized operations
%[text] - assembly of the master stiffness matrix and consistent nodal force vector using [sparse matrices](https://www.mathworks.com/help/matlab/math/computational-advantages-of-sparse-matrices.html) \
%[text] It is highly recommended to vectorize math-operations in MATLAB for better runtime performance. The most prominent operations to be vectorized are these that involve `for`-loops over the problem's dimensions that can be arbitrarily large. As it is shown in the provided mlx-files the runtime gain increases with the problem's size.
%[text] Operations of [`sparse`](https://www.mathworks.com/help/matlab/ref/sparse.html) matrices and vectors prove to be very efficient as compared to operations involving their full counterparts. This is especially true for problems involving matrices and vectors that are actually sparse, namely, they are dominated by zero components. Moreover there is a size limitation for the the construction of full matrices in MATLAB, see [this MATLAB Answers article](https://www.mathworks.com/matlabcentral/answers/850320-what-is-the-maximum-size-of-a-matrix#:~:text=The%20maximum%20number%20of%20rows%20or%20columns%20in%20MATLAB%20is,array%20must%20fit%20into%20memory.) for more information accordingly. MATLAB [`sparse`](https://www.mathworks.com/help/matlab/ref/sparse.html) matrices allow for significantly larger matrices since the zero entries do not need to be stored. This increases also the efficiency of matrix-vector computations between sparse matrices significantly.
%[text] In Live Script [`main_AppexA_EfficiencyCondiderationsTimoshenkoBeam.mlx`](file:./AppexA_EfficiencyConsiderations/main_AppexA_EfficiencyCondiderationsTimoshenkoBeam.mlx) it is shown how page-wise finite element implementation of the Timoshenko beam problem can be developed and it is also demonstrated the runtime performance gain. The script takes about one minute to execute.
%[text] Similarly, in Live Script [`main_AppexA_EfficiencyCondiderationsReissnerMindlinPlate.mlx`](file:./AppexA_EfficiencyConsiderations/main_AppexA_EfficiencyCondiderationsReissnerMindlinPlate.mlx) it is shown how page-wise finite element implementation of the Reissner-Mindlin plate problem can be developed following the same workflow as the one demonstrated for the Timoshenko beam problem. It is herein also highlighted the substantial runtime benefit obtained by the combination of the page-wise implementation together with the use of the [`sparse`](https://www.mathworks.com/help/matlab/ref/sparse.html) matrices. This script takes about 25 minutes to execute.

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"inline","rightPanelPercent":15.6}
%---
