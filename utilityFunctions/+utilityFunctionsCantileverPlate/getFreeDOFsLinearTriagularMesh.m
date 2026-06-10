%[text] %[text:anchor:H_1CEB6968] # **Retrieval of the free degrees of freedom corresponding to the linear triangular mesh**
%[text:tableOfContents]{"heading":"**Table of Contents**"}
%[text] It returns the global numbering of the free Degrees of Freedom (DOFs) corresponding to the bilinear quadrilateral mesh created in section [Three-noded linear triangular mesh](file:..\..\2_ReissnerMindlinPlate\main_Chapter2_AShearDeformablePlateElement.mlx:H_FA1680CB) in Live Script [`main_Chapter2_AShearDeformablePlateElement.mlx`](file:..\..\2_ReissnerMindlinPlate\main_Chapter2_AShearDeformablePlateElement.mlx). This function is only relevant for the setup detailed in the latter Live Script.
function [freeDOFsT1, homDOFsT1] = getFreeDOFsLinearTriagularMesh ...
    (mshT1, xCoordDirichlet, yCoordDirichlet)
%[text] %[text:anchor:H_6C4EB3DB] ## **Function Description**
%[text] It returns the global numbering of the DOFs which do not lie on the Dirichlet boundary corresponding to the linear triangular mesh for the problem at hand, namely, the cantilever rectangular plate.
%[text] %[text:anchor:H_C2A6403E]  **Input :**
%[text]  `mshT1` : Structure containing information about the linear triangular mesh:
%[text]   `.nodes` : Two-dimensional array containing the nodes  of the triangular mesh
%[text]  `.elements` : Two-dimensional array containing the elements by means of the corresponding nodal IDs of the triangular mesh
%[text]  `xCoordDirichlet` : Two-dimensional array containing the spans in the form \[xStart xEnd\] regarding the x-coordinates along each Dirichlet boundary
%[text]  `yCoordDirichlet` : Two-dimensional array containing the spans in the form \[yStart yEnd\] regarding the y-coordinates along each Dirichlet boundary
%[text] %[text:anchor:H_CBAFD59D]  **Output :**
%[text]  `freeDOFsT1` : The global numbering of the free DOFs corresponding to the linear triangular mesh
%[text]  `homDOFsQ1` : The global numbering of DOFs along the Dirichlet boundary corresponding to the linear triangular mesh
%[text] %[text:anchor:H_894F0B02] ## **Function Implementation**
%[text] %[text:anchor:H_D88B8E8E] ### Input validation
    arguments
        mshT1 (1, 1) {mustHaveNodesAndElements}
        xCoordDirichlet (1, 2) double
        yCoordDirichlet (1, 2) double
    end
%[text] %[text:anchor:H_CC3575DE] ### Read input
    numNodesT1 = height(mshT1.nodes);
    freeDOFsT1 = 1:3*numNodesT1;
    homNodesT1 = [];
    if height(xCoordDirichlet) ~= height(yCoordDirichlet)
        error("Input arrays 'xCoord_Dirichlet' and 'yCoord_Dirichlet' " + ...
            "must have the same number of rows (same number of " + ...
            "Dirichlet boundaries)");
    end
    numDrchltBndrs = height(xCoordDirichlet);
    tol = 1e-6;
    counter = 1;
%[text] %[text:anchor:H_11A0816A] ### Loop over all the Dirichlet boundaries
    for ii = 1:numDrchltBndrs
        if (xCoordDirichlet(ii, 1) ~= xCoordDirichlet(ii, 2) && ...
                yCoordDirichlet(ii, 1) ~= yCoordDirichlet(ii, 2)) || ...
                (xCoordDirichlet(ii, 1) == xCoordDirichlet(ii, 2) && ...
                yCoordDirichlet(ii, 1) == yCoordDirichlet(ii, 2))
            error("The Dirichlet boundary extension x = [%d, %d] and " + ...
                "y = [%d, %d] for the %d-th boundary is wrong", ...
                xCoordDirichlet(ii, 1), xCoordDirichlet(ii, 2), ...
                yCoordDirichlet(ii, 1), yCoordDirichlet(ii, 2), ii);
        end
        isOnX = false;
        if xCoordDirichlet(ii, 1) ~= xCoordDirichlet(ii, 2)
            isOnX = true;
        end
        if isOnX
            span = xCoordDirichlet;
            fixed = yCoordDirichlet(ii, 1);
        else
            span = yCoordDirichlet;
            fixed = xCoordDirichlet(ii, 1);
        end
    
        for jj = 1:numNodesT1
            coord = mshT1.nodes(jj, :);
            if isOnX
                if (coord(1, 1) > span(1, 1) - tol && ...
                        coord(1, 1) < span(1, 2) + tol) && ...
                        abs(coord(1, 2) - fixed) < tol
                    homNodesT1(counter) = jj;
                    counter = counter + 1;
                end
            else
                if (coord(1, 2) > span(1, 1) - tol && ...
                        coord(1, 2) < span(1, 2) + tol) && ...
                        abs(coord(1, 1) - fixed) < tol
                    homNodesT1(counter) = jj;
                    counter = counter + 1;
                end
            end
        end
    end
    
    homDOFsT1 = [3*homNodesT1 - 2
                  3*homNodesT1 - 1
                  3*homNodesT1];
    homDOFsT1 = reshape(homDOFsT1, [], 3*numel(homNodesT1));
    homDOFsT1 = sort(homDOFsT1);
    [~, idx] = ismember(homDOFsT1, freeDOFsT1);
    freeDOFsT1(idx) = [];
%[text] 
end

%[appendix]{"version":"1.0"}
%---
