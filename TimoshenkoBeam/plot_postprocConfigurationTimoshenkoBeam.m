%[text] %[text:anchor:H_33BBA0DF] # Graphical representation of the deformed configuration for the Timoshenko beam problem
%[text:tableOfContents]{"heading":"**Table of Contents**"}
%[text] The function returns a graphical representation of the current configuration for the Timoshenko beam problem given the Finite Element mesh and the discrete solution vector. Moreover, the rotation field is displayed on the provided Finite Element mesh directly, if selected.
function h = plot_postprocConfigurationTimoshenkoBeam(msh, uh, resultant, choiceShapeFunctionsPrimal)
%[text] %[text:anchor:H_74525E29] ## Function description
%[text] Plots the deformed configuration for a Timoshenko beam.
%[text]  **Input :**
%[text]  `msh` : Struct-variable representing the Finite Element mesh containing the following information:
%[text]                                                                          .nodes : Two-dimensional array containing the Cartesian coordinates of each node in the mesh
%[text]                                                                     .elements : Two-dimensional array containing the global numbering of the nodes that belong to each element
%[text]  `uh` : Solution vector corresponding to the Finite Element solution
%[text]  `resultant` : "displacement" or "rotation"
%[text]  `choiceShapeFunctionsPrimal` : Flag regarding the choice of shape functions for the discretization of the primal field
%[text]  `1` : Linear basis functions
%[text]  `2` : Quadratic basis functions
%[text]  **Output :**
%[text]  `h` : Handle for the generated figure
%[text] %[text:anchor:H_20AA9D71] ## Function implementation
%[text] %[text:anchor:H_927AACB4] ## Check input
if ~strcmp(resultant, "displacement") && ~strcmp(resultant, "rotation")
    error("Input variable 'resultant' should be selected either as 'displacement' or 'rotation'")
end
%[text] %[text:anchor:H_F404E500] ## Read input
numNodes = numel(msh.nodes(:, 1));
numDOFs = 2*numNodes;
if strcmp(resultant, "displacement")
    startIdx = 1;
elseif strcmp(resultant, "rotation")
    startIdx = 2;
end
%[text] %[text:anchor:H_30C67C8C] ## **Plot the deformed configuration taking into account the polynomial order of the chosen interpolation**
switch choiceShapeFunctionsPrimal
    case 1
        h = plot(msh.nodes, uh(startIdx:2:numDOFs, 1), '-o', LineWidth=2);
    case 2
        [nodes, idx] = sort(msh.nodes);
        uSorted = uh(startIdx:2:numDOFs, 1);
        uSorted = uSorted(idx);
        h = plot(nodes, uSorted, '-o');
    otherwise
        error("This function only supports two-noded linear and three-noded quadratic basis functions")
end
xlabel("$x$", Interpreter="latex");
ylabel("$w$", Interpreter="latex");
if strcmp(resultant, "displacement")
    title("Deformation")
elseif strcmp(resultant, "rotation")
    title("Cross-sectional rotation $\beta$", Interpreter="latex")
end
%[text] 
end

%[appendix]{"version":"1.0"}
%---
