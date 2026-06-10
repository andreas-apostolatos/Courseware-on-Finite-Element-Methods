%[text] %[text:anchor:T_6C324AFB] # Plot Two-Dimensional Finite Element mesh and Basis Functions
%[text] Displays the two-dimensional basis functions on a two-dimensional mesh.
%[text] - It is tested and works for the two-dimensional bilinear basis functions
%[text] - The visual representation of the basis functions correspond to the parametric space. However, the parametric space is replaced in an one-on-one correspondance with the physical space for visualization purposes. Therefore, the visual representation of the two-dimensional basis functions in the physical space is not valid for skewed elements
%[text] - This function should be only used for regular uniform meshes to get an accurate representation of the two-dimensional basis functions in the physical space
%[text] - The function can be take long time to compute even for meshes with few elements \
%[text:tableOfContents]{"heading":"**Table of Contents**"}
function [h1, h2, h3] = plotMeshAndBasisFunctions2d(msh, N, isPlotDerivatives)
%[text] %[text:anchor:H_EF551ED2] ## Function Description
%[text] Plot a two-dimensional mesh together together with the basis functions
%[text]  **Input :**
%[text]  `msh` : Nodal coordinates of the end vertices of the finite element
%[text]                                               .`nodes` : Array containing the coordinates of the nodes in the mesh
%[text]                                         .`elements` : Array containing the IDs of the nodes of each element in the mesh
%[text]  `N` : Cell array containing the basis functions in a symbolic form
%[text]  `isPlotDerivatives` : Flag indicating whether the derivatives should also be plotted
%[text] %[text:anchor:H_5361C5B6]  **Output :**
%[text]  `h1` : Figure handle for the figure containing the basis functions
%[text]  `h2` : Figure handle for the figure containing the derivatives of the basis functions along X
%[text]  `h3` : Figure handle for the figure containing the derivatives of the basis functions along Y
%[text] %[text:anchor:H_781C3337] ## Function Implementation
%[text] %[text:anchor:H_05C6349E] ### Input validation
    arguments
        msh (1, 1) {mustHaveNodesAndElements}
        N (1, :) symfun
        isPlotDerivatives (1, 1) logical
    end
%[text] %[text:anchor:H_564604D6] ### Read input
    syms xi eta NPhys(X, Y)
    if isPlotDerivatives
        syms dNPhysdX(X, Y) dN_physdY(X, Y)
    end
    minX = min(msh.nodes(:, 1));
    maxX = max(msh.nodes(:, 1));
    minY = min(msh.nodes(:, 2));
    maxY = max(msh.nodes(:, 2));
    epsil = 1e-6;
%[text] %[text:anchor:H_34D827DE] ### Loop over all the nodes in the mesh
    for ii = 1:length(msh.nodes)
%[text] %[text:anchor:H_775C8358] #### Plot the mesh with the node numberings
        if ~exist('h1', 'var')
            h1 = figure;
        else
            figure(h1.Number);
        end
        trimesh(msh.elements, msh.nodes(:, 1), msh.nodes(:, 2), msh.nodes(:, 3), ...
                    'FaceColor', [217 218 219]/255, EdgeColor='k'); % , visible='off'
        axis equal;
        hold on;
%[text] %[text:anchor:H_8A5A993F] #### Plot the nodes in the mesh
        plot3(msh.nodes(:, 1), msh.nodes(:, 2), msh.nodes(:, 3), 'ok');
%[text] %[text:anchor:H_B09DF157] #### Plot the node numberings on the mesh
        for jj = 1:length(msh.nodes)
            if ii == jj
                textColor = [203 65 84]/255;
            else
                textColor = 'k';
            end
            text(msh.nodes(jj, 1), msh.nodes(jj, 2), msh.nodes(jj, 3) + max(msh.nodes, [], "all", "includenan")/20, ...
                sprintf("X_{%d}", jj), "Color", textColor);
        end
%[text] %[text:anchor:H_A1CA5C47] #### Elements at which the node belongs to
        [idElmnts, idy] = find(msh.elements == ii);
        idElmnts_test = ones(4, 1);
        idyTest = ones(4, 1);
        numElmnts = length(idElmnts);
        idElmnts_test(1:numElmnts) = idElmnts;
        idyTest(1:numElmnts) = idy;
        isElement = false(4, 1);
        isElement(1:numElmnts) = true;
%[text] %[text:anchor:H_5E226581] #### Coordinates of the vertices of each element where the node is found
        nodesEl1 = msh.nodes(msh.elements(idElmnts_test(1), :), 1:2);
        nodesEl2 = msh.nodes(msh.elements(idElmnts_test(2), :), 1:2);
        nodesEl3 = msh.nodes(msh.elements(idElmnts_test(3), :), 1:2);
        nodesEl4 = msh.nodes(msh.elements(idElmnts_test(4), :), 1:2);    
%[text] %[text:anchor:H_E9D4F01A] #### Transformation of the basis functions in the physical space
        NmapXY{1} = subs(N{idyTest(1)}, {xi, eta}, {2/(max(nodesEl1(:, 1)) - min(nodesEl1(:, 1)))*X - (max(nodesEl1(:, 1)) + min(nodesEl1(:, 1)))/(max(nodesEl1(:, 1)) - min(nodesEl1(:, 1))), 2/(max(nodesEl1(:, 2)) - min(nodesEl1(:, 2)))*Y - (max(nodesEl1(:, 2)) + min(nodesEl1(:, 2)))/(max(nodesEl1(:, 2)) - min(nodesEl1(:, 2)))});
        NmapXY{2} = subs(N{idyTest(2)}, {xi, eta}, {2/(max(nodesEl2(:, 1)) - min(nodesEl2(:, 1)))*X - (max(nodesEl2(:, 1)) + min(nodesEl2(:, 1)))/(max(nodesEl2(:, 1)) - min(nodesEl2(:, 1))), 2/(max(nodesEl2(:, 2)) - min(nodesEl2(:, 2)))*Y - (max(nodesEl2(:, 2)) + min(nodesEl2(:, 2)))/(max(nodesEl2(:, 2)) - min(nodesEl2(:, 2)))});
        NmapXY{3} = subs(N{idyTest(3)}, {xi, eta}, {2/(max(nodesEl3(:, 1)) - min(nodesEl3(:, 1)))*X - (max(nodesEl3(:, 1)) + min(nodesEl3(:, 1)))/(max(nodesEl3(:, 1)) - min(nodesEl3(:, 1))), 2/(max(nodesEl3(:, 2)) - min(nodesEl3(:, 2)))*Y - (max(nodesEl3(:, 2)) + min(nodesEl3(:, 2)))/(max(nodesEl3(:, 2)) - min(nodesEl3(:, 2)))});
        NmapXY{4} = subs(N{idyTest(4)}, {xi, eta}, {2/(max(nodesEl4(:, 1)) - min(nodesEl4(:, 1)))*X - (max(nodesEl4(:, 1)) + min(nodesEl4(:, 1)))/(max(nodesEl4(:, 1)) - min(nodesEl4(:, 1))), 2/(max(nodesEl4(:, 2)) - min(nodesEl4(:, 2)))*Y - (max(nodesEl4(:, 2)) + min(nodesEl4(:, 2)))/(max(nodesEl4(:, 2)) - min(nodesEl4(:, 2)))});
%[text] %[text:anchor:H_258C2EAC] #### Symbolic function representation of the basis function at the current node in the physical domain
        NPhys(X, Y) = isElement(1)*(X > min(nodesEl1(:, 1)) + epsil & X < max(nodesEl1(:, 1)) - epsil & Y > min(nodesEl1(:, 2)) + epsil & Y < max(nodesEl1(:, 2)) - epsil)*NmapXY{1} + ...
                isElement(2)*(X > min(nodesEl2(:, 1)) + epsil & X < max(nodesEl2(:, 1)) - epsil & Y > min(nodesEl2(:, 2)) + epsil & Y < max(nodesEl2(:, 2)) - epsil)*NmapXY{2} + ...
                isElement(3)*(X > min(nodesEl3(:, 1)) + epsil & X < max(nodesEl3(:, 1)) - epsil & Y > min(nodesEl3(:, 2)) + epsil & Y < max(nodesEl3(:, 2)) - epsil)*NmapXY{3} + ...
                isElement(4)*(X > min(nodesEl4(:, 1)) + epsil & X < max(nodesEl4(:, 1)) - epsil & Y > min(nodesEl4(:, 2)) + epsil & Y < max(nodesEl4(:, 2)) - epsil)*NmapXY{4};
%[text] %[text:anchor:H_69E71C0B] #### Plot the basis function in the whole computational domain
        fsurf(NPhys, [minX maxX minY maxY], LineStyle='none', FaceAlpha=3/4); % [217 218 219]/255
        plot3([msh.nodes(ii, 1); msh.nodes(ii, 1)], [msh.nodes(ii, 2); msh.nodes(ii, 2)], [msh.nodes(ii, 3); msh.nodes(ii, 3) + 1], ...
            "Color", [203 65 84]/255, "LineWidth", 2);
        xlabel("X");
        ylabel("Y");
        zlabel(sprintf("N_{%d}", ii));
        axis equal;
        hold off;
        pause(1)
%[text] %[text:anchor:H_2A5BECA6] #### Display the node numberings on the mesh for the derivatives of the basis functions with respect to *X*
        if isPlotDerivatives
            if ~exist('h2', 'var')
                h2 = figure;
            else
                figure(h2.Number);
            end
            trimesh(msh.elements, msh.nodes(:, 1), msh.nodes(:, 2), msh.nodes(:, 3), ...
                    'FaceColor', [217 218 219]/255, EdgeColor='k'); % , visible='off'
            axis equal;
            hold on;
%[text] %[text:anchor:H_AAEA1933] #### Display the node numberings on the mesh
            for jj = 1:length(msh.nodes)
                if ii == jj
                    textColor = [203 65 84]/255;
                else
                    textColor = 'k';
                end
                text(msh.nodes(jj, 1), msh.nodes(jj, 2), msh.nodes(jj, 3) + max(msh.nodes, [], "all", "includenan")/20, ...
                    sprintf("X_{%d}", jj), "Color", textColor);
            end
        end
%[text] %[text:anchor:H_22AA7BEE] #### Derivatives of the basis functions with respect to the physical space
        if isPlotDerivatives
            dNmapXYdX = num2cell(diff(NmapXY, X));
            dNmapXYdY = num2cell(diff(NmapXY, Y));
        end
%[text] %[text:anchor:H_82CCAA2F] #### Derivatives of the basis functions with respect to *X* in the whole computational domain
        if isPlotDerivatives
%[text] %[text:anchor:H_CBBFB32D] #### Formulation of the derivative of the basis functions with respect to *X*
            dNPhysdX(X, Y) = isElement(1)*(X > min(nodesEl1(:, 1)) & X < max(nodesEl1(:, 1)) & Y > min(nodesEl1(:, 2)) & Y < max(nodesEl1(:, 2)) )*dNmapXYdX{1} + ...
            isElement(2)*(X > min(nodesEl2(:, 1)) & X < max(nodesEl2(:, 1)) - epsil & Y > min(nodesEl2(:, 2)) & Y < max(nodesEl2(:, 2)) )*dNmapXYdX{2} + ...
            isElement(3)*(X > min(nodesEl3(:, 1)) & X < max(nodesEl3(:, 1)) - epsil & Y > min(nodesEl3(:, 2)) & Y < max(nodesEl3(:, 2)) )*dNmapXYdX{3} + ...
            isElement(4)*(X > min(nodesEl4(:, 1)) & X < max(nodesEl4(:, 1)) - epsil & Y > min(nodesEl4(:, 2)) & Y < max(nodesEl4(:, 2)) )*dNmapXYdX{4};
%[text] %[text:anchor:H_B8C5ED36] #### Plot the basis functions with respect to *X* along the whole computational domain
            fsurf(dNPhysdX, [minX+epsil*1e2 maxX-epsil*1e2 minY+epsil*1e2 maxY-epsil*1e2], ...
                LineStyle='none', FaceAlpha=3/4); % [217 218 219]/255
            xlabel("X");
            ylabel("Y");
            zlabel(sprintf("^{\\partial N_{%d}}/_{\\partial X}", ii));
            axis equal;
            hold off;
            pause(1)
        end
%[text] %[text:anchor:H_577FDD2E] #### Plot the node numbering on the mesh for the derivatives of the basis functions with respect to *Y*
        if isPlotDerivatives
            if ~exist('h3', 'var')
                h3 = figure;
            else
                figure(h3.Number);
            end
            trimesh(msh.elements, msh.nodes(:, 1), msh.nodes(:, 2), msh.nodes(:, 3), ...
                    'FaceColor', [217 218 219]/255, EdgeColor='k'); % , visible='off'
            axis equal;
            hold on;

%[text] %[text:anchor:H_A8A3BF87] #### Plot the node numberings on the mesh
            for jj = 1:length(msh.nodes)
                if ii == jj
                    textColor = [203 65 84]/255;
                else
                    textColor = 'k';
                end
                text(msh.nodes(jj, 1), msh.nodes(jj, 2), msh.nodes(jj, 3) + max(msh.nodes, [], "all", "includenan")/20, ...
                    sprintf("X_{%d}", jj), "Color", textColor);
            end
        end
%[text] %[text:anchor:H_B29ED222] #### Plot the derivatives of the basis functions with respect to *Y* in the whole computational domain
        if isPlotDerivatives            
%[text] %[text:anchor:H_BB682488] #### Formulation of the derivative of the basis functions with respect to *Y*
            dN_physdY(X, Y) = isElement(1)*(X > min(nodesEl1(:, 1)) & X < max(nodesEl1(:, 1)) & Y > min(nodesEl1(:, 2)) & Y < max(nodesEl1(:, 2)) )*dNmapXYdY{1} + ...
            isElement(2)*(X > min(nodesEl2(:, 1)) & X < max(nodesEl2(:, 1)) & Y > min(nodesEl2(:, 2)) & Y < max(nodesEl2(:, 2)) )*dNmapXYdY{2} + ...
            isElement(3)*(X > min(nodesEl3(:, 1)) & X < max(nodesEl3(:, 1)) & Y > min(nodesEl3(:, 2)) & Y < max(nodesEl3(:, 2)) )*dNmapXYdY{3} + ...
            isElement(4)*(X > min(nodesEl4(:, 1)) & X < max(nodesEl4(:, 1)) & Y > min(nodesEl4(:, 2)) & Y < max(nodesEl4(:, 2)) )*dNmapXYdY{4};
%[text] %[text:anchor:H_D84E67E7] #### Plot of the basis functions with respect to *Y* in the whole computational domain
            fsurf(dN_physdY, [minX+epsil maxX-epsil minY+epsil maxY-epsil], ...
                LineStyle='none', FaceAlpha=3/4); % [217 218 219]/255
            xlabel("X");
            ylabel("Y");
            zlabel(sprintf("^{\\partial N_{%d}}/_{\\partial Y}", ii));
            axis equal;
            hold off;
            pause(1)
        end
%[text] %[text:anchor:H_DE62B779] #### Pause the graph for 2 seconds for inspection
        pause(2)
    end
end

%[appendix]{"version":"1.0"}
%---
