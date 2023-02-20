function [GP, GW] = getGaussPointsAndWeightsOnCanonicalTriangle(polOrder)
%% GETGAUSSPOINTSANDWEIGHTSONCANONICALTRIANGLE Returns the Gauss points and weights on a canonical triangle based on the polynomial order that the integrand attains
%
% Returns the Gauss Points and Weights for the exact integration of a
% polynomial of degree n = p + q, where p is the polynomial degree in the
% u-direction and q is the polynomial degree in the v-direction. The 
% function works up to n = 8. The Gaussian integration is exact over 
% the canonical element:
%   
%         (0,1)
%           |\
%           | \
%           |  \
%           |   \
%           |    \
%           |     \
%           |      \
%           ---------
%         (0,0)     (1,0)
%
% Source : Quadrature Formulas in Two Dimensions, 
%          Math 5172 - Finite Element Method, Section 001, Spring 2010
%
%    Input :
% polOrder : The polynomial degree of the integrand
%
%   Output :
%       GP : The set of Gauss Points with their three linear dependent 
%            coordinates GPi = [1-lambda1 - lambda2 lambda1 lambda2]
%       GW : The respective set of Gauss weights 
%
% Function layout :
%
% 0. Input validation
%
% 1. Assign the linear independent Gauss Point coordinates and their respective weights
%
% 2. Add the linear dependent coordinate to the Gauss Point array
%
%% Function Implementation

%% 0. Input validation
arguments
    polOrder (1,1) double {mustBeInteger, mustBeNonnegative}
end

%% 1. Assign the linear independent Gauss Point coordinates and their respective weights
if polOrder == 1
    GPLI = [0.33333333333333 0.33333333333333];
    GW = 1.00000000000000;
elseif polOrder == 2
    GPLI  =  [0.16666666666667 0.16666666666667 
              0.16666666666667 0.66666666666667 
              0.66666666666667 0.16666666666667];
    GW = [0.33333333333333
          0.33333333333333
          0.33333333333333];
elseif polOrder == 3
    GPLI  =  [0.33333333333333 0.33333333333333
              0.20000000000000 0.20000000000000 
              0.20000000000000 0.60000000000000 
              0.60000000000000 0.20000000000000];
    GW = [-0.56250000000000
          0.52083333333333
          0.52083333333333
          0.52083333333333];
elseif polOrder == 4
    GPLI  =  [0.44594849091597 0.44594849091597 
              0.44594849091597 0.10810301816807 
              0.10810301816807 0.44594849091597 
              0.09157621350977 0.09157621350977 
              0.09157621350977 0.81684757298046 
              0.81684757298046 0.09157621350977];
    GW = [0.22338158967801
          0.22338158967801
          0.22338158967801
          0.10995174365532
          0.10995174365532
          0.10995174365532];
elseif polOrder == 5
    GPLI  =  [0.33333333333333 0.33333333333333 
              0.47014206410511 0.47014206410511 
              0.47014206410511 0.05971587178977 
              0.05971587178977 0.47014206410511 
              0.10128650732346 0.10128650732346 
              0.10128650732346 0.79742698535309 
              0.79742698535309 0.10128650732346];
    GW = [0.22500000000000
          0.13239415278851
          0.13239415278851
          0.13239415278851
          0.12593918054483
          0.12593918054483
          0.12593918054483];
elseif polOrder == 6
    GPLI  =  [0.24928674517091 0.24928674517091 
              0.24928674517091 0.50142650965818 
              0.50142650965818 0.24928674517091 
              0.06308901449150 0.06308901449150 
              0.06308901449150 0.87382197101700 
              0.87382197101700 0.06308901449150 
              0.31035245103378 0.63650249912140 
              0.63650249912140 0.05314504984482 
              0.05314504984482 0.31035245103378 
              0.63650249912140 0.31035245103378 
              0.31035245103378 0.05314504984482 
              0.05314504984482 0.63650249912140];
    GW = [0.11678627572638
          0.11678627572638
          0.11678627572638
          0.05084490637021
          0.05084490637021
          0.05084490637021
          0.08285107561837
          0.08285107561837
          0.08285107561837
          0.08285107561837
          0.08285107561837
          0.08285107561837];
elseif polOrder == 7
    GPLI  =  [0.33333333333333 0.33333333333333 
              0.26034596607904 0.26034596607904 
              0.26034596607904 0.47930806784192 
              0.47930806784192 0.26034596607904 
              0.06513010290222 0.06513010290222 
              0.06513010290222 0.86973979419557 
              0.86973979419557 0.06513010290222 
              0.31286549600487 0.63844418856981 
              0.63844418856981 0.04869031542532 
              0.04869031542532 0.31286549600487 
              0.63844418856981 0.31286549600487 
              0.31286549600487 0.04869031542532 
              0.04869031542532 0.63844418856981];
    GW = [-0.14957004446768
          0.17561525743321
          0.17561525743321
          0.17561525743321
          0.05334723560884
          0.05334723560884
          0.05334723560884
          0.07711376089026
          0.07711376089026
          0.07711376089026
          0.07711376089026
          0.07711376089026
          0.07711376089026];
elseif polOrder == 8
    GPLI  =  [0.33333333333333 0.33333333333333 
              0.45929258829272 0.45929258829272 
              0.45929258829272 0.08141482341455 
              0.08141482341455 0.45929258829272 
              0.17056930775176 0.17056930775176 
              0.17056930775176 0.65886138449648 
              0.65886138449648 0.17056930775176 
              0.05054722831703 0.05054722831703 
              0.05054722831703 0.89890554336594 
              0.89890554336594 0.05054722831703 
              0.26311282963464 0.72849239295540 
              0.72849239295540 0.00839477740996 
              0.00839477740996 0.26311282963464 
              0.72849239295540 0.26311282963464 
              0.26311282963464 0.00839477740996 
          0.00839477740996 0.72849239295540];
    GW = [0.14431560767779
          0.09509163426728
          0.09509163426728
          0.09509163426728
          0.10321737053472
          0.10321737053472
          0.10321737053472
          0.03245849762320
          0.03245849762320
          0.03245849762320
          0.02723031417443
          0.02723031417443
          0.02723031417443
          0.02723031417443
          0.02723031417443
          0.02723031417443];
else
    error('Gauss Points and weights for the integration of polynomial integrands of order higher than 8 over triangular domain have not been yet implemented');
end

%% 2. Add the linear dependent coordinate to the Gauss Point array
numGP = length(GPLI(:,1));
GP = zeros(numGP,3);
GP(:,2:3) = GPLI;
GP(:,1) = ones(numGP,1) - GP(:,2) - GP(:,3);

%% 3. Multiply all Gauss weights with the area of the canonical triangle
GW = GW.*0.5;

end
