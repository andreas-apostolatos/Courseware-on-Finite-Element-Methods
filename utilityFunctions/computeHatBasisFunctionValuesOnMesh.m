%[text] %[text:anchor:T_586CF7E9] # Evaluation of the one-dimensional Hat-Basis Functions at a Set of Points
%[text:tableOfContents]{"heading":"**Table of Contents**"}
%[text] Returns the values of the hat basis functions (see Eqs. (1) below) at a number of evaluation points on a provided mesh
%[text] $\\begin{array}{l}\nN\_1 \\left(\\xi \\right)=\\frac{1-\\xi }{2},\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\left(1\\ldotp 1\\right)\\\\\nN\_2 \\left(\\xi \\right)=\\frac{1+\\xi }{2}\\ldotp \\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\left(1\\ldotp 2\\right)\n\\end{array}${"editStyle":"visual"}
function N = computeHatBasisFunctionValuesOnMesh(phi1, phi2, msh, numPts)
%[text] %[text:anchor:H_01DCB2B4] ## Function description
%[text] Values of the linear one-dimensional hat-basis functions at a set of points.
%[text]  **Input :**
%[text]  `phi1`, `phi2` : Symbolic expression for the basis functions
%[text]  `msh` : Nodal coordinates of the end vertices of the finite element
%[text]                               .`nodes` : Array containing the coordinates of the nodes in the mesh
%[text]                          .`elements` : Array containing the IDs of the nodes of each element in the mesh
%[text]  `numPts` : Number of evaluation points
%[text]  **Output :**
%[text]  `N` : Array containing the value of each basis function at each evaluation point
%[text] %[text:anchor:H_5530BC53] ## Function Implementation
%[text] %[text:anchor:H_E1129D23] ### Input validation
arguments
    phi1 (1, 1) symfun
    phi2 (1, 1) symfun
    msh (1, 1) {mustHaveNodesAndElements(msh)}
    numPts double {mustBeInteger, mustBePositive}
end
%[text] %[text:anchor:H_96F71D08] ### Declaration of auxiliary symbolic variables
    syms X xi    
%[text] %[text:anchor:H_097EE2CB] ### Generation of evaluation points in the space
    pts = linspace(msh.nodes(1), msh.nodes(end), numPts);
%[text] %[text:anchor:H_49897CB9] ### Number of nodes in the mesh
    numNodes = numel(msh.nodes);
%[text] %[text:anchor:H_4255AA7D] ### Initialization of the matrix containing the values of the basis functions
    N = zeros(numNodes, numPts);
%[text] %[text:anchor:H_77E99C48] ### Loop over all the nodes in the mesh
    for ii = 1:numNodes
%[text] %[text:anchor:H_84BCC713] #### Element that contain the node
        [idElmnts, idy] = find(msh.elements == ii);
        numElmnts = length(idElmnts);
%[text] %[text:anchor:H_0FE9EEEF] #### Loop over all the elements where the node belongs to
        for jj = 1:numElmnts
%[text] Verify whether the node is the first or the last (second for the one-dimensional two-noded linear element)
            if idy(jj) ~= 1 && idy(jj) ~= 2
               error("The node must be either the first or the second for the one-dimensional two-noded linear element")
            end
%[text] Coordinates of the left and right nodes in the element
            xL = msh.nodes(msh.elements(idElmnts(jj), 1));
            xR = msh.nodes(msh.elements(idElmnts(jj), 2));
%[text] Loop over all the evaluation points between the left and the right nodes in the element
            for kk = 1:numPts
                if pts(kk) >= xL && pts(kk) <= xR
%[text] Setup of the geometric map                    
                    x = [phi1 phi2]*[xL; xR];
%[text] Inversion of the geometric map to obtain the parametric coordinate in the element's parametric space
                    Xxi = solve(X == x(xi), xi);
%[text] Computation of the parametric location of the evaluation points in the element's parametric space
                    Xxi_eval = subs(Xxi, pts(kk));
%[text] Evaluation of the corresponding hat shape function depending on the position of the node in the element
                    if idy(jj) == 1
                        N(ii, kk) = subs(phi1, xi, Xxi_eval);
                    elseif idy(jj) == 2
                        N(ii, kk) = subs(phi2, xi, Xxi_eval);
                    end
%[text] 
                end
%[text] 
            end
%[text] 
        end
%[text] 
    end
%[text] 
end

%[appendix]{"version":"1.0"}
%---
