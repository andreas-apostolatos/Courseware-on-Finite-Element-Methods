%[text] %[text:anchor:T_50930D3E] # Plot One-Dimensional Mesh With the Corresponding Geometric Mapping
%[text] Given an one-dimensional mesh, the function generates three graphs concerning an one-dimensional mesh in a tiled format: The one-dimensional basis functions, the provided parametric location mapped onto the physical space and the distribution of the geometric Jacobian in the parametric space of the provided element.
%[text:tableOfContents]{"heading":"**Table of Contents**"}
function t = plot1dMeshAndGeometryMapping(el, xi_, msh, N)
%[text] %[text:anchor:H_7B8A5F7D] ## **Function description**
%[text] Plots a graph displaying the basis functions and their values in the chosen parametric location. Moreover, the mapping of the selected parametric location is displayed in the physical space by means of the geometric mapping.
%[text] %[text:anchor:H_FC688161]  **Input :**
%[text]  `el` : The ID of the element for which the parametric location is sought to be displayed in the physical space
%[text]  `xi_` : The selected parametric location
%[text]  `msh` : An one-dimensional mesh as a structure with the following fields:
%[text]                              .`nodes` : Number of nodes in the one-dimensional mesh
%[text]                        .`elements` : Number of elements in the one-dimensional mesh
%[text]  `N` : Array containing the symbolic expression of the one-dimensional basis functions
%[text] 
%[text]  **Output :**
%[text]  `t` : Handle to the tiledlayout graphical object
%[text] %[text:anchor:H_E47A3A5A] ## **Function implementation**
%[text] %[text:anchor:H_529D13C1] ### Input validation
arguments
    el (1, 1) double {mustBeInteger, mustBePositive}
    xi_ (1, 1) double {mustBeInRange(xi_, -1, +1)}
    msh (1, 1) {mustHaveNodesAndElements(msh)}
    N (:, 1) symfun
end
%[text] %[text:anchor:H_2C657D5C] ### Read input
    if ~ismember(el, msh.elements(:, 1))
        error("The selected element does not belong to the selected mesh")
    end
%[text] %[text:anchor:H_728C6838] ### Node numbering of the corner element nodes
    node_ids = msh.elements(el, :);
    if numel(node_ids) ~= numel(N(xi_))
        error("The number of nodes per element does not match the number of " + ...
            "basis functions");
    end
%[text] %[text:anchor:H_E0841A6D] ### Nodes of the provided element
    nodes_el = msh.nodes(node_ids, :);
    XL = node_ids(1);
    XR = node_ids(2);
    numNodes = numel(nodes_el);
%[text] %[text:anchor:H_BA25D1ED] ### Geometry parametrization
    syms X(xi)
    X(xi) = N(xi)*nodes_el;
%[text] %[text:anchor:H_1C2620F0] ### Physical location of the selected parametric location through the geometric mapping
    X_ = double(X(xi_));
%[text] %[text:anchor:H_71F8722E] ### Geometric Jacobian
    syms Jmtx(xi)
    Jmtx(xi) = diff(X(xi), xi);
%[text] %[text:anchor:H_8072F9D4] ### Initialization of the tiled layout
    t = tiledlayout(3, 1);
%[text] %[text:anchor:H_CC3AFE40] #### Graphical representation of the basis functions with their values at the provided parametric location
    nexttile
    fplot(N, [-1 +1])
    hold on;
    N_value_vct = zeros(numNodes, 1);
    N_value_vct(end) = max(N(xi_));
    plot(ones(numNodes, 1)*xi_, N(xi_), "ko", "MarkerSize", 5);
    plot(ones(numNodes, 1)*xi_, N_value_vct, "-.k", "MarkerSize", 5);
    for ii = 1:numNodes
        filter_mtx = zeros(1, numNodes);
        filter_mtx(ii) = 1;
        text(xi_ + 0.06, filter_mtx*transpose(N(xi_)), ...
            num2str(double(filter_mtx*transpose(N(xi_)))));
    end
    hold off;
    xlabel("\xi");
    ylabel(sprintf("N_i^{L%i}", numNodes));
    title("Parametric space")
%[text] %[text:anchor:H_66CD1AD5] ### Graphical representation of the provided parametric location in the physical space
    nexttile
    plot(transpose(nodes_el), zeros(1, numel(nodes_el)), "-ko", "MarkerSize", 5, ...
        "MarkerFaceColor", "k");
    hold on;
    plot(X_, 0, "ko", "MarkerSize", 5);
    text(X_ - 0.06, -0.3, num2str(double(X_)));
    text(XL - 0.06, 0.3, "X1");
    text(XR - 0.06, 0.3, "X2");
    xlabel("x");
    title("Physical space")
    hold off;
%[text] %[text:anchor:H_1C1A99FA] ### Graphical representation of the geometric Jacobian in the parametric space of the provided element
    nexttile
    J_value_vct = zeros(numNodes, 1);
    J_value_vct(end) = Jmtx(xi_);
    fplot(Jmtx, [-1 +1]);
    hold on;
    plot(xi_, Jmtx(xi_), "ko", "MarkerSize", 5);
    plot(ones(numNodes, 1)*xi_, J_value_vct, "-.k", "MarkerSize", 5);
    plot([-1; +1], [0; 0], "-.k", "MarkerSize", 5);
    xlabel("\xi");
    ylabel("$J = \frac{d x}{d \xi}$", "Interpreter","latex");
    title("Jacobian")
    hold off;
%[text] 
end

%[appendix]{"version":"1.0"}
%---
