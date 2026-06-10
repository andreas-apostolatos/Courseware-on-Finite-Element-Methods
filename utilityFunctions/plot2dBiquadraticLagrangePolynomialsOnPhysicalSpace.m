%[text] %[text:anchor:T_A9279F0F] # Plot Biquadratic Lagrange Polynomials in Physical Space
%[text] Generates three graphs concerning a two-dimensional biquadratic mesh in a tiled format: The two-dimensional basis functions, the provided parametric location mapped onto the physical space and the distribution of the geometric Jacobian in the parametric space of the provided element.
%[text:tableOfContents]{"heading":"**Table of Contents**"}
function t = plot2dBiquadraticLagrangePolynomialsOnPhysicalSpace ...
    (el, msh, numVertices, l1L2, l2L2, l3L2, ij)
%[text] %[text:anchor:H_FAB31DA5] ## **Function description**
%[text] Displays the bilinear basis functions on the two-dimensional biquadratic quadrilateral mesh in the physical space.
%[text]  **Input :**
%[text]  `el` : The element for which the parametric location is sought to be displayed in the physical space
%[text]  `msh` : The one-dimensional biquadratic mesh using four-noded elements as a structure with the following fields:
%[text]                                           .`nodes` : Number of nodes in the one-dimensional mesh
%[text]                                     .`elements` : Number of elements in the one-dimensional mesh
%[text]  `numVertices` : Number of vertices in the bilinear quadratic mesh
%[text]  `l1_L2`, `l2_L2`, `l3_L2` : The three Lagrange polynomials for the one-dimensional quadratic element
%[text]{"align":"center"} ![quadratic_lagange_polynomials.png](text:image:53d3)
%[text]  `ij` : Legend entries for the plot, which should be an array with nine elements
%[text] 
%[text]  **Output :**
%[text]  `t` : Handle to the tiledlayout graphical object
%[text] %[text:anchor:H_E0217A59] ## **Function implementation**
%[text] %[text:anchor:H_EEDCB170] ### **Input validation**
    arguments
        el (1, 1) double {mustBeInteger, mustBeGreaterThanOrEqual(el, 1)}
        msh (1, 1) {mustHaveNodesAndElements}
        numVertices (1, 1) double {mustBeInteger, mustBePositive}
        l1L2 (1, 1) symfun
        l2L2 (1, 1) symfun
        l3L2 (1, 1) symfun
        ij (1, 9)
    end
%[text] %[text:anchor:H_DAC56352] ### Read input
    syms l_Q2(x, y, x1, x2, x3, y1, y2, y3)
    if floor(el) ~= el
        error("Variable 'el' should be an integer");
    end
    if el < 1 || el > numel(msh.elements(:, 1))
        error("Variable 'el' is chosen as %i but it should be an integer in the close interval [1, %i]", el, numel(msh.elements(:, 1)));
    end
    t = tiledlayout(3, 3);
    numNodes_el = numel(msh.elements(el, :));
%[text] %[text:anchor:H_C1A89326] ### Nodes of the selected element in the quadrilateral patch
    Xel_Q2 = msh.nodes(msh.elements(el, :), :); %[text:anchor:H_6203C08B]
%[text] %[text:anchor:H_426FB92C] ### Loop over all the nodes in the quadrilateral patch
    for ii = 1:numNodes_el
%[text] %[text:anchor:H_D722A666] #### Set-up the corresponding Lagrange polynomial for the current node
        switch ii
            case 1
                idx1 = 1;
                idx2 = 7;
                idx3 = 4;
                idy1 = 1;
                idy2 = 5;
                idy3 = 2;
                l_Q2(x, y) = l1L2(x, Xel_Q2(idx1, 1), Xel_Q2(idx2, 1), Xel_Q2(idx3, 1))*l1L2(y, Xel_Q2(idy1, 2), Xel_Q2(idy2, 2), Xel_Q2(idy3, 2));
            case 2
                idx1 = 2;
                idx2 = 8;
                idx3 = 3;
                idy1 = 1;
                idy2 = 5;
                idy3 = 2;
                l_Q2(x, y) = l1L2(x, Xel_Q2(idx1, 1), Xel_Q2(idx2, 1), Xel_Q2(idx3, 1))*l3L2(y, Xel_Q2(idy1, 2), Xel_Q2(idy2, 2), Xel_Q2(idy3, 2));
            case 3
                idx1 = 2;
                idx2 = 8;
                idx3 = 3;
                idy1 = 4;
                idy2 = 6;
                idy3 = 3;
                l_Q2(x, y) = l3L2(x, Xel_Q2(idx1, 1), Xel_Q2(idx2, 1), Xel_Q2(idx3, 1))*l3L2(y, Xel_Q2(idy1, 2), Xel_Q2(idy2, 2), Xel_Q2(idy3, 2));
            case 4
                idx1 = 1;
                idx2 = 7;
                idx3 = 4;
                idy1 = 4;
                idy2 = 6;
                idy3 = 3;
                l_Q2(x, y) = l3L2(x, Xel_Q2(idx1, 1), Xel_Q2(idx2, 1), Xel_Q2(idx3, 1))*l1L2(y, Xel_Q2(idy1, 2), Xel_Q2(idy2, 2), Xel_Q2(idy3, 2));
            case 5
                idx1 = 5;
                idx2 = 9;
                idx3 = 6;
                idy1 = 1;
                idy2 = 5;
                idy3 = 2;
                l_Q2(x, y) = l1L2(x, Xel_Q2(idx1, 1), Xel_Q2(idx2, 1), Xel_Q2(idx3, 1))*l2L2(y, Xel_Q2(idy1, 2), Xel_Q2(idy2, 2), Xel_Q2(idy3, 2));
            case 6
                idx1 = 5;
                idx2 = 9;
                idx3 = 6;
                idy1 = 4;
                idy2 = 6;
                idy3 = 3;
                l_Q2(x, y) = l3L2(x, Xel_Q2(idx1, 1), Xel_Q2(idx2, 1), Xel_Q2(idx3, 1))*l2L2(y, Xel_Q2(idy1, 2), Xel_Q2(idy2, 2), Xel_Q2(idy3, 2));
            case 7
                idx1 = 2;
                idx2 = 8;
                idx3 = 3;
                idy1 = 7;
                idy2 = 9;
                idy3 = 8;
                l_Q2(x, y) = l2L2(x, Xel_Q2(idx1, 1), Xel_Q2(idx2, 1), Xel_Q2(idx3, 1))*l1L2(y, Xel_Q2(idy1, 2), Xel_Q2(idy2, 2), Xel_Q2(idy3, 2));
            case 8
                idx1 = 1;
                idx2 = 7;
                idx3 = 4;
                idy1 = 7;
                idy2 = 9;
                idy3 = 8;
                l_Q2(x, y) = l2L2(x, Xel_Q2(idx1, 1), Xel_Q2(idx2, 1), Xel_Q2(idx3, 1))*l3L2(y, Xel_Q2(idy1, 2), Xel_Q2(idy2, 2), Xel_Q2(idy3, 2));
            case 9
                idx1 = 5;
                idx2 = 9;
                idx3 = 6;
                idy1 = 7;
                idy2 = 9;
                idy3 = 8;
                l_Q2(x, y) = l2L2(x, Xel_Q2(idx1, 1), Xel_Q2(idx2, 1), Xel_Q2(idx3, 1))*l2L2(y, Xel_Q2(idy1, 2), Xel_Q2(idy2, 2), Xel_Q2(idy3, 2));
            otherwise
                error("The biquadratic quadrilateral should only have 9 nodes");
        end
%[text] %[text:anchor:H_1C148C34] #### Plot the mesh along with the corresponding Lagrange polynomial associated to the current node
        nexttile
        plotMesh2d(msh, numVertices, true);
        hold on;
        fsurf(l_Q2(x, y), [min(Xel_Q2([idx1 idx2 idx3], 1)) max(Xel_Q2([idx1 idx2 idx3], 1)) min(Xel_Q2([idy1 idy2 idy3], 2)) max(Xel_Q2([idy1 idy2 idy3], 2))], ...
        "FaceColor",[217, 218, 219]./255, "FaceAlpha", 0.5);
        xlabel("x");
        ylabel("y");
        zlabel(sprintf("l_{%s}", num2str(ij(ii))))
        axis equal;
        hold off;
    end
    title(t, sprintf("Biquadratic basis functions in the physical space of element %i", el));
%[text] 
end

%[appendix]{"version":"1.0"}
%---
%[text:image:53d3]
%   data: {"align":"baseline","height":290,"src":"data:image\/png;base64,iVBORw0KGgoAAAANSUhEUgAAAjAAAAGkCAIAAACgjIjwAAA0iklEQVR42u2dbWhc172v86Uf+uX2wz0fTnpKm9yWglUjIWSaahBGilt8gn0MbnzzchLUoQQjIdRcx4OxMKoSXIKvm8QyxjmBa5zmCjt2xohmeoNiBk2FHbdygo1wE2Ui1ZGZ0olrWRGisaOqGnT\/1kp2xqN53bNf1tr7eRjEjLwtba3Zs579W\/u\/1r5nBQAAQAPuoQkAAAAhAQAAICQAAEBIAAAACAkAABASAAAAQgIAAIQEAACAkAAAACEBAAAgJAAAQEgAAAAICQAAEBIAAABCAgAAhAQAAICQAAAAIQEAACAkAABASAAAAAgJAAAQEgAAAEICAACEBAAAgJAAAAAhAQAAICQAAEBIAAAACAkAABASANTNrfcvLKROq8fN0y9YD+ubsoE8aCiAsAjp93\/+9Gtdb3wfwDUebv7us23fkceZf\/\/WRw\/fK4\/z\/\/Fv8pCXLz\/4bfVQG6iH9R3ZQB6ypfpfavvuB+6XHygPGhZcIvnS\/Uuzf0RI\/vD1n\/2XaMnHHZAjwN8W8H0HAtYIEmsk6GR++VMlEnmiok\/FuFN+H1Si+uTI0+onf9z1Q\/WTHUxRHIo0wmImLkIiIfnGvz74s\/ajl\/gA0AvU+RMsCYkq7Hmi1n0o8FP9cuJQpBEW3nls139+DyH5+faLkHwMSfQCRjeCikG++8CKZbInYilxFYcin8daWZr943yyTYd3IdRCevXdrL8hCYxDen\/LQ7rtm9jIMhNlEVBTPFrMxBGS\/+cj\/oYkwEMu7a0ykzxZ+luGtw\/Kx6O5N+\/TJKeGXUiEJCiPpA1ze3bZbQITVBOPEJIWQhLu+9UFQhIYHYkq9DhfDuUF4G8Bl+IRQtJFSBKSoq9PcmhCfqqQR\/BSRai0ZOJMIO9b6bPLMRWPEJIuQiIkQX4qkl472NddrEHIwAuJHa6IFY8QkkZvPyEJFaky7vCUAARmQBIh1ROPbqcPISQd334JSTNzn9M1h41b718I6gBdNahBPHuzl+jfTd\/h\/HiEkPR6+599+2NCUji745BXoImNlJIDlg4RUsV4JA+EpO\/bT0gKD6pygaqzADcIQqoYj3K3\/oKQ9H37JSERkgKP5KGPHr4XFQVeS1p1rx9++OGLL76ozw6vjUcIScfjVUISvVKAUQuSMku0PGpBcdNH8LTqXpPJ5I9\/\/GN9dnhtPEJIOh6vhKSgspA6TTAKVVRCSKW4nT60Nh4hJE2PV0JSIIMRxQv22s3cqISQaopHCEnT41US0qvvZumMCEZgdFRCSEVZzMSLxiOEpO\/xSkgKBmpJAoJR\/ZgodYRUlPlkW9F4hJD0PV7bj14iJJmOWuWadnA2Khk0fJf\/AX\/27Y99eegmpDLxCCFpfQJFSDK96yQYOY5a1cKUqKS5kKanp733gcSjpdk\/IiTzhERIMhS1RBt3pXMPVTpv+gfcYwqE9MYbb2zfvt3jHS4fjxCS7scrIcnEvpJhOs8yqObW17Z7ffrpp7ds2eK9kMrHI4Sk+\/FKSDIItf5CIBcJ1dlJOje4tt3rH\/7wh5s3b3ospIrxCCEZcLwSkoxALQ\/KRSPv0fmSks7dq\/dCqhiPEFJdZDKZZDKZTqdLbTA9PS0bXLt2rZ7WJyTpj6rt5qIRTkJI9cQjhGSfRCIRiURisVhHR8fg4ODaDQ4ePNja2iobbN68+ciRI\/W0PiFJZ9RSArSDz93r6Rc0fBcQUk3xCCHZZHl5ubm5WRVNzs3NNTU1zczM5G9w5cqV9evXZ7N3ks3i4qJIS75ju\/UJSZybQzVO0i2nIqSa4hFCskkqlRLHWC97e3uHhobyNxgeHu7u7rZeSk769a9\/XU\/rE5J0Q5UwcNFIQyfpU+bA\/ZBqikcIySYFvunr6+vv78\/f4K233tq2bZv1cufOnXv27Cna+haEJLNsRAmDtmR++VNNYitCWqm6uK7KnhAhFSEej\/f09Fgv961y1xnB\/PwDDzxw8ODBixcv\/uY3v1FXm+p8+wlJ2AiqRJOhVIRUUzwiIdkkkUh0dXXlJ6SBgYGCbaanp0VaTzzxxJEjRw4cOFBgLButT0jSAcq7cRJCcjYeIaR6GR8fb2trs16KnERR+Rv8\/e9\/n5iYyN\/gzJkz9bc+IUkHG9EOOAkhuRGPEJJNcrmcCGlsbEyeT01NNTY2zs7OynORkKqsk68NDQ3Xr1+X55cvX96wYcPCwkL9rU9I8hF1wZx2wEkIyaV4hJDqCkmRSKSzs7OlpWVkZER9MxqNxuNx9fy1115rbm5+8skn29vbZWOnWp+Q5Atq6ivtYCI+1jiEXEi1xiOEZN7bT0jyxUZMfTXdSb6seKtV9\/rhhx+++OKLnu2wxKOFdx4LnsIR0l38\/s+fEpI87suwUQCQkOS9k8J8x1gb8QghGRmQCUnYCIxwUmiFZC8eISQjhURIwkZg20leXgsMrZDsxSOEZKSQ7jTKM6OEJLdtxH32cBJCsheP5t68LwAthpAISWQj8MJJ3pxthFNItuMRQjJVSCtcSXINaupCcs7hQS14CIVkY+4RQgqCkFaYk+TOuTM2wklufMA\/evheXx6lhPTBBx\/cunXLcR\/UE48QktlCIiQ5CysDhQ2313HQMyHdvHkzEolEo9H29vaC2+LUucN1xiOEZLaQCEkOotbwph1wUlCFZHHs2DHrHtbf+MY3HNzhOuMRQjJeSNHXJ599+2O6FUdsxBre4STzy5+6dE8\/zbvXbDb77W9\/26kdvp0+VGc8QkjGC4mQ5IiNuPcrOcmNA0DzW5j\/6Ec\/SqVSTu3w3Jv35W79BSGFXUgSkuRBn2IbbAQuOUnb7vXq1as\/\/OEPk8mkUzvsSDxCSEEQEiEJG0H9qJKWpb9lAi+kTCbzgx\/84IMPPnBwhx2JRwgpIEIiJNn8ZLp28QBwkrbda3d399e\/\/vVvfkn9OyzZyJF4hJACIiQVkmbmPqdPqSkb+X6La9ANZxcWCsn9kJyKRwgpOEIiJNWajZgAC6Wc5NTCQmEQkoPxCCEFR0iEpJp6HGwE5c9XHEnPYRCSg\/EIIQVKSISkamACLFTppPqvLwZeSJKNbqcPadslIiSfW5+QVBHK6qBK6i9wCLyQbN9mAiGFQkivvpslJJXvYrARVEn9RXfBFpLj8Qgh1RfqM5lkMplOp0ttMDMzIxtMTk562foSkn7\/50\/pTYq8X57ccQCCRJ2XG4MtJMfjEUKyTyKRiEQisViso6NjcHBw7QbHjx9XG2zevHnfvn2etT4hqSjc5QhsHzm2i+4CLKSFdx5bzMQRkhYsLy83NzdPT0\/fOU2Ym2tqapIwlL9BLpdraGiYmpq6884tLMjzojnJpdYnJBV+eLivBNSB7RXBtepeP\/zwwxdffNGRHV6a\/aMb8Qgh2SSVSkkwsl729vYODQ0VCGndunWZzJ3R56WlpfXr109MTBRtfQtnQ1L70Uv0IwrWToX6kXht4xAK6h1jHY9HLvWEYRHS8PBwd3e39bKvr6+\/v79gm1OnTm3dunVwcPCRRx7Zv3+\/x8crIckCG4EjpzU2ChwCKST34hEJySbxeLynp8d6uW+Vgm1isZio6MSJEzt37oxGo2vvFuxq6xOSrBNbVqsDR7Ax8BtIIbl09Qgh2SeRSHR1deUnpIGBgfwNRkdHN23atLy8rF6KkKxbNHrW+iKkkIck1gcCZ6l1VaHgCUni0Xyyzb2dREh2GB8fb2v76l0ROYmi8jc4c+ZM\/pie6GrPnj0et76EpDDfloIVGcCls5zqCxyCJyRX4xFCskkulxMhjY2NyfOpqanGxsbZ2Vl5PjExkc1m5cnk5KR88+rVqyurVXZbtmwRRXnf+mEOSVw6Apeofnp1\/gf8dvqQLw8HhSTxSITkatsiJPshKRKJdHZ2trS0jIyMqG9Go9F4\/IvTh1OnTsk\/qQ2ef\/55X1pfbBTOkGSvJgqgGqq\/mKSzkP70pz\/985\/\/rKlHEhuJkxBSYPGg9SUkvfpuNmw24tIRuEqVF5P0HLL79NNP29radu7cKV8L7mJeZocXM3G34xFCCr6QVkJ2g3MuHYFn5z0VLybp2b2OjY2dPHlyZXWq7MMPP1zlDs8n29yORwgpFEIKVUji0hF4RsWLSTp3r9PT093d3QVrN5TaYYlHDt6FDyGFWkjhCUlcOgIvqXgxSefu9ciRIw899NChQ4eq2WFv4hFCCouQwnDvPi4dgfeUX3pV\/+71m9\/8Zn5pQ9Edvp0+5E08QkhhEdJK0O\/dx6Uj8As58EqtBqJn93poFXny+eef\/8u\/\/Esulyu\/w87epBwhIaTghyQuHYG\/Tiq6zJ2e3eunn3764IMPPvnkk21tbadOnSq\/w5KNPItHCClEQloJ7oqrLFgH\/lLqYpLO3es\/\/vGP\/GxUaoe9jEcIKVxCCuS9++q8pyeAIxS9mGT6DfrcuEk5QkJIXxGwxYS4dAT6sHbc2GghuXqbCYSEkL4ISUG6LQWXjkC306P8i0lGC8ntdVQREkJabdBnRoMRkmpadxnAAwqWFPq+gfgYjxBSGIUUjNtSMFgH2p4nBaDExoN1VBESQvqCAFxJYrAOtMXGzc61wpt1VBESQvoC029LwWAd6IzplZ+eLRSEkBDSVyHJ0BVXF1KnqfMG\/c+ZDB248zEeIaTwCmnF2BVXGawDIzB04G7uzfv8ikcIKdRCMnExIQbrgCjvHl6uo4qQEFKRkGTQiqsM1gHnT27HIy8XCkJICOkuzFpMiME6MA6DBu68XygIITl6+pPJJJPJdDpd5ERjbu69u5mZmdGz9U1ZcZUVVMFEKt7ETxP8mgmLkJwhkUhEIpFYLNbR0TE4OFjwr2fPnm3OY926dQMDA3q2vhGLCbGCKpiLEQN3viwUhJCcYXl5WTQzPT2twlBTU1PRAKQ4f\/78xo0b5+fntW19zefJ3nr\/wkcP32v0TEMIORKSdB5tlng0n2zTYU8Qkh1SqZQEI+tlb2\/v0NBQ8c701q22trZz587p3Pqaz5Pl0hGYzkLqtBzGOscjH0u9EVK9DA8Pd3d3Wy\/7+vr6+\/uLbjk4OPjUU0+Vaf2CNQ19DEl6zpNlsA6CgbYDd\/7OhNWtJzRSSPF4vKenx3q5b5Ui7\/TiYlNT05UrV4w4HdAzJDFYB4FBz4o7HxcKIiE5QyKR6Orqyk9IRWsWfvvb327dutWU1o++Pvns2x9zUgngEhpOpPN9JixCcoDx8fG2tq+uAYqcRFFrN9u1a9fhw4cNan2t5skyDRaChxzSWl0Q9X0mLEJygFwuJ0IaGxuT51NTU42NjbOzs\/J8YmIim\/3qSkxra6vaxpTW12oxIc0LkwBsoNV9vCQbaRWPEFJdISkSiXR2dra0tIyMjHzRoUej8XjckpY07o0bN8xqfU3myX5y5On8O28CBAZ9jm0dZsIiJFq\/HDrMk+VusBBsdBi402QmLEKi9Svg+zxZ3cbZARyWgd\/rCekzExYh0foV8HeeLBOPICQhyccKUn1mwiIkWr+qkOTXPFkmHkFI8Gtakg4zYRESrV8bvoQklvSG8ODXxAZ\/7wmLkBCSHbwvAWfiEYQN7y+XaljqjZBo\/WpDkpfVDSyiCmHD+4JSDUu9ERKtXxVeloCLjVglCEKIl9OS9Cz1Rki0frV4UwKu7nhE3wThxJuBO21LvRESrV8t3pSAM\/EIwow305K0LfVGSLR+bSHJ1RJwahkA3K4v1bnUGyHR+rXhakhiEVWAFZenJelc6o2QaP3acK8E\/ObpF1hEFcDVz4Lmpd4Iida3E5LcqG5gXQYAC5cupup20yOEROvXixsl4NwQFiAfN66nLrzz2O30IbpEhBQoIa04XQLOPSYA3A5JRpR6IyRa3w5io3ueGdXzgwcQDJw9UTOllgEh0fo2Q5IjJeCUegOUOVdzZCjblFJvhETr2+e+X12Ymfu8zh9CLQNAGRwpATcuHiEkWr9mnn374zpLwD858jS1DABlqL8E3KBSb4RE69cbkmxXN7BsHUCVIcn2RVYJRpqv6o2QnCeTySSTyXQ6XTIyz82Njo6Oj48HTEivvpu1vXYDt+ADqIZ6FrjTf1VvhOQwiUQiEonEYrGOjo7BwcG1G4yNjbW2tu7evXvHjh1PPPFELpcLjJBW7FY3UMsA4PbZm4m1DAipLpaXl5ubm6enp1UMampqmpmZKdhAbHTx4kX1csuWLSMjI0ES0oqtBe5Ytg6g1o9Mrf9lPtlmXC0DQqqLVColwch62dvbOzQ0lL\/B6OioBKNqWt\/CuEOn1gXuWLYOoFZqvX2fubUMpvSEOgppeHi4u7vbetnX19ff35+\/wZkzZ3bt2iXfbGxslCx17NgxQ08HKoak6qsbKPUGsBeSqhxXMLeWgYRUF\/F4vKenx3q5b5X8DQ4cONDQ0HDq1Cl5nk6nN2zYcP78+eAJqfoF7rhDOYA9FlKnqyxMNbeWASHVRSKR6Orqyk9IAwMD+RucPHnyoYceyt9ACJ6QVqqrbmDZOgC3Q5LRtQwIqS7Gx8fb2r5aslDkJIrK3+Ds2bP5QloboQIjpGrucU6pN0CdIaliearRtQwIqS5yuZwIaWxsTJ5PTU01NjbOzs7K84mJiWz2TlxYWlp64IEHUqnUymoZ3saNG4vORgqAkFYqVTfImR2l3gB1Un4x4tvpQybWMiAkJ0NSJBLp7OxsaWmxSrqj0Wg8\/sUY7nvvvdfe3v7oo4\/KBkePHjW09aukTHUDq3oD1E\/5cW+DbsGHkGh91ylV3cBMWAAHQ1LRoW+zbsGHkGh9Lyh6+z5mwgI4yNqQZNwt+BASre8Fa6sbmAkL4HhIKpg+YeI9JhASre9RSHr27Y+tl8yEBXA1JN1OHwpAqTdCovXdwrp9X61LngBANeR\/sgJTy4CQaH1XsKobiEcA7oUk+XAFqZYBIdH6biFCev1\/\/YyFggBcYiF1+q8v\/Chgg3UIidZ3hclzo9wTFsBVbpz4zu3JN+gSERJCqkDmlz\/dtfdATXemAIDq+exy7MZQJJAz\/BASre8k1nzyotOSAKBOrHtMBHINFIRE6zscj9SH5NV3szZuKQsA5Vl45zE18WghdTp4i+gjJFrfuY\/K3QsFVXNnCgConoJ7TAQvJCEkWt8x1i4UZE1LAoD6KViXQT5uASsgQki0vivxSPHs2x9T3QDgzEes2MSjgN1sDCHR+m7FIyskUd0AUCelFlEN2O2YERKt71Y8UlRzS1kAqBiPSi2iGqSQhJBofRfjkaJg0VUAqImKi6gGJiQhJFrfxXhkQXUDgG0qLqIamJCEkGh9d+ORotQtZQGgwglfdYuoBiMkISRav954VOVtJpiWBFArBROPAh+SEJJ9MplMMplMp9PFU\/bc3Ht5LCwsBFJINd2knIE7gJqo6YawAQhJCMkmiUQiEonEYrGOjo7BwcG1Gxw7dqyhoaH5S86fPx88IVUfjxRMSwKo4fNV4x2PAhCSEJIdlpeXxTHT09MqCTU1Nc3MzBRss2vXrhMnTpje+g7GIwWLrgJUQ6mJR8EOSQjJDqlUSoKR9bK3t3doaKhgm82bN4+Pj4uulpaWAimkWuORgmlJANUgNqp+sC4wIQkh2WF4eLi7u9t62dfX19\/fXxCh1q1bt2XLltbWVnmyb9++Uq1vEYZ4pIi+PsnAHUAZPrsck4ftD6aJHjKlJ9RRSPF4vKenx3q5b5X8Df76179KbJKv8vz69esbN248efJkkBKSvXhkwcAdQCmsOx7Zw+iQhJDskEgkurq68hPSwMBAme3379+\/e\/fuIAnJdjxScLckgJJne6VXCQpwSEJI9hkfH29r++p6o8hJFJW\/wbVr1yRFWS\/7+\/v37NkTGCHVGY+skMR6QgAFVFwlqBrk42loSEJIdsjlciKksbExeT41NdXY2Dg7OyvPJyYmstk7cz\/T6XRDQ4Mqw7t+\/XokEglS2becfy39LVP\/z2EhcIACKq4SFOyQhJDshyTRTGdnZ0tLy8jIiPpmNBq1gtGJEyeam5tlA\/l6\/PhxQ1vfpXikYOAO4K4P1zuPLWbijvwo+ZDePP0CQgqLkELS+mv56OF7HYlHCgbuABTVrxIU4JCEkGj9GpBzLqfikQXrCQGs1LhKUJUhyfFPK10iQgpsPFJIQmIhcMBGNa0SFNSQhJBofT\/jkYKFwCHM2FslqMqQZNaVJIRE69dwtuV4PLKgugHCHI+cHawzNyQhJFq\/KhwsrisKC4FDOKl1Se9aMWvhBoRE61d7nlXP0gxVhiSmJUGocG+wztCQhJBoff\/jkYKFwCFsuDpYZ2JIQki0vhbxSMHAHYQHtwfrLOTDa0pIQki0fuV4JGdYnv06Bu4gDHgzWJcfkrw5p0RICClQhzIDdxAGvBmsMy4khVdIP18lm80+cDetra20vl\/xSMHAHQQbzwbrjAtJ4RXSjVWsJ\/nQ+vkHsS+XQ7mDHwQVjwfr\/D25REgmoXnr+xvzGbiDQOLxYJ1ZISngQrpWCVpfw3ikYOAOgocvg3UGhaSAC2nTpk3fL826deto\/TLx6KOH7\/V3H1jjDoKEX4N1+Xg2hQMhFWG5ErR+KTS5CzI3p4DA4ONgXX5I0vmeFFxDovWL43s8UjBwB8HA38E6U0ISQqL1i8cjfVatv+eZUQbuwGgkGEk80kWNGockhETr6xuPLBi4A6PRYbCuICS5dysZhISQHI5Hup1AcVdZMBd9But0\/owjpHrJZDLJZDKdTpffbGJiotRkWz1b3437lNePCEm0RO8GZrGYiYuQNNwxPVcSQkg2SSQSkUgkFot1dHQMDg6W2mx6enr9+vXiLVNa3737lNcP666Cceg2WKd5SEJIdlheXm5ubhbZ3Dng5uaamppmZmbWbra0tLRt27b29naDhKTt4LLw6rtZlm8Ag5BsJAlJ293TMCQhJDukUikJRtbL3t7eoaGhtZsdOHDg8OHDO3fuLCMkCy0+P6nTupUzFMDAHZjC7fQhPQfrLPS5cZ9uPaFhQhoeHu7u7rZe9vX19ff3F2xz8eLF7du3y5PyQtLtjEn\/9YBZdxX0R9V55279RfP91C0kISQ7xOPxnp4e6+W+Ve6KGgsLmzdvVuN4pgjJiMWAVxi4AxOQbKTnpaO1IUmrc1CEZIdEItHV1ZWfkAYGBvI3ED\/t2rVrbJUdO3YMDg4WLcbTqvVNuaekEH19kuUbQFs+uxyThxG7Kh95rU5DEZIdxsfH29q+WiRR5CSKyt9ADLTzS1pbW7dv3378+HGdW9+UG0paUHEHeqLDCqrmnokiJDvkcjkRkqQfeT41NdXY2Dg7O7uyOuUomy1c5MaIITtNllKtHnWnc5ZvAN3Qts67FFqN1SMk+yEpEol0dna2tLSMjIyob0aj0Xg8bqKQNC+uKwrrroJ2nbt+izJUgz7VTAiJ1td3HZGKcMMk0AdtF2WoiD7T4RESra\/pWkFVwsAdaIJxg3UFIYkuESH53\/qa37CrIlSBgw4YOlhnockdZxBS2FvfiMmw5aEKHPxF\/0UZTAlJCCnUrW96PLKgChz8wrg671LosJIQQgp16xs0GbY8VIGDXxh96SgfHSYjIqTwtr5xk2HLw038wIcxBsMvHel2hoqQwiskfdb6dQrWAgcvCcalo7v86vckWYQUXiGZOBm2IlxMAm8wZT3vWvG3ygkhhVRImlR5Og5V4OBRmDBkPW8bIcnHQieEFFIhGT0ZtjwsKQQe2MiU9bzthSSEhJC8Q5+VQlyCJYXAPQJT510KH4dPEFIYhSRnQEGNRxZUgYNLBKbOW8OQhJBCJ6TATIYtDxeTwJWPzzuPLWbigf8z\/SrBRUihE1JgJsNWhCWFwFmCV+ddCr8mKSKkcAlJznoCWe1dCglJXEwCRwj8paMCpKPw\/swVIYVLSMGbDFuNk7iYBPUThktHBSev3k+SRUjhElKo4pGCi0ngQO8crCWCqsT7SbIIKURCMvfOsHXCxSSo00YhuXRUgPfzQxBSiIQUhmrvUnAxCewRtktHazsNhGSGkDKZTDKZTKfTpTaQf5INZmZmdGj9kFR7l3cSy9xBrYTt0lFhL+ftVWeEZJNEIhGJRGKxWEdHx+Dg4NoNXnrppZ\/85Cd79+7dtGnTK6+84nvrh6fauxRcTIKaT+MCumBd9Xhc\/42Q7LC8vNzc3Dw9PX3nBGpurqmpqSAGTU1NrV+\/fn5+Xp7fuHFj3bp1spmPrR+wWx\/ZhmXuoCYbhfPSkY\/nsgjJDqlUSoKR9bK3t3doaCh\/g1wup3QliJakla9fv+5j639y5OmwVXuXgnsmQTUsZuLY6Asxe1j\/jZDsMDw83N3dbb3s6+vr7+8vGqROnTq1bdu2w4cPl2p9C1d3OITV3mXgYhKUR93rKOSDdfm4Xf\/tWU8YTCHF4\/Genh7r5b5V1m5248aN11577amnnnrkkUfU8J0vpwOBX9u7VriYBOWZT7Zho3w8mzGCkOyQSCS6urryE9LAwECZ7Ts7O4sWPnjT+mGu9i4FF5OgFOGcA1tNN4KQNBXS+Ph4W9tXUxNETqKo\/A2uXr2af1Vpz549e\/fu9aX1qfYueWA9M8rFJCggPMun1oo39d8IyQ65XE6ENDY2trJaUNfY2Dg7OyvPJyYmstms+mZDQ4NoSZ7LP0UikdHRUV9an2rvMrQfvcTFJLBQc2Bzt\/5CU6zFm0pdhGQ\/JIlmOjs7W1paRkZG1Dej0Wg8\/sW9Uk6ePNnU1PTzn\/9cvvo1D0mOIcoZysPSq2BBIYPvZ7cIyU\/cbn2qvStCgQMoQnLnvbqayP36b4QUZCFJPKKcoSIUOABzYKvE7QophBRYIVHOUD3Mlg0zFDJUj3QpN0+\/gJAQkp1zGcoZanISBQ4hJOSLedvrWBASQqo5Hnl\/t0fTocAhhFDIUCuuljYgpGAKKYS3Kq8fChzCBot563ayi5CCKSSqve1BgUOobMSlI9vdi0shCSEFUEhuX3gMNhQ4hAEKGerBveUxEVIAhUS1d51wv\/NgQyFD\/bhU2oCQgiYkqr2dchIFDkGFQob6cam0ASEFTUgsXucIFDgEFQoZHEE6GTdKGxBSoITE4nUO8uzbH+Ok4NmIW0s4hRulDQgpUEKinMFZoq9PUnQXJBtRyOAgbpQ2IKRACYlyBseh6C4YUFbnBo6XNiCk4AiJcgaXoOjOdCircwnHr1gjpOAIiXIGV53ESnfm2oiyOpdwvLQBIQVESJQzuIoquqMQ3ESwkas4u4gzQgqIkLgXn9tQdGcilNW5jbOlDQgpIEKinMEbJ1F0Z5aNKGTwJiQhJISU98GjnMErKLozBcrqPMPBewsgpCAIiXIGL7nvVxdwkuZQVuclDpY2IKQ6NJDJJJPJdDpdaoPp6WnZ4NKlS662PuUMvuQkiu50thGFDN6HJEfOiRGSTRKJRCQSicViHR0dg4ODazfYv3+\/\/JNssG3btscff3xxcdGl1qecwa+cRNGdnmAj73GqtAEh2WF5ebm5uVkC0J2jf26uqalpZmYmf4PJycn169fPz8+rl1u3bo3H4y61PuUMvkAhODaCfBwpbUBIdkilUpJ+rJe9vb1DQ0P5G2Sz2QsXLuRvcPjw4aKtb2FvTyhn8BEKwXVj4Z3HFjNx2sEX6iltqL8nDLWQhoeHu7u7rZd9fX39\/f2lNr527ZqkJclMbpwOUM6Ak8CyEVOOfEQ6ovpDEkKyQzwe7+npsV7uW6XoltevX29vb3\/55ZfdaH3KGTRxEpOTdLARRd46hKQ6z48Rkh0SiURXV1d+QhoYGFi72ZUrV1pbW48fP+5S61POoAlMTvKXzy7HsJEWpwV1X0FASHYYHx9va\/tqloPISRRVsM2FCxc2bNhw9uxZ91qfcgacBEyA1Yo6R+0Qkh1yuZwIaWxsTJ5PTU01NjbOzs7K84mJiWz2zn0KMplMc3NzKpVa+pLl5WVnW59yBt1gwqz3MAFWN+pctQEh2Q9JkUiks7OzpaVlZGREfTMajary7gMHDnz\/bp577jlnW9\/ZRXbBqZzEhFmPbUSRt1bUWdqAkPzEdus7UtACLuUknOSNjZhypG1Isn2ujJCMFBLlDNrChFlvwEbaUs\/VBIRkpJAoZ8BJ2Ih20Bbb4zcIyTwhUc6gP2rCLE5y5QSc5Ri0x3ZpA0IyT0iszoCTsBHojO0bUiAk84TE6gxmOYl2cNBGn12O0Q5GYK8MGCEZJqRPjjzNeB1OwkagOfZuSIGQDBOSnHdQzoCTQmgjlmMwMSQhpCALaSF12qlbBYPHTmIBVmwUNmxc7UZIJgmJ6Uc4CRuBMe9d7SfQCMkkIVHOYDTtRy\/hpJpgGW\/TqXXUDiEZIySmHwXDSSzAWr2N5t68j3YwGumybp5+ASEFUEispoqTyEZgFrWuuomQzBASN4fFSdgITKSm0gaEZIaQKGcIGNw8CRuFhJquNSAkM4TEaqo4CRuBoVQ\/aoeQDBAS5QxBhbE7bBQGqh+1Q0gGCKnOuwIDTsJG4CPVr7WKkAwQEuUMOAkbgdFIJ1ZNSEJIuguJ8TqchI3AdG6efqGaCUkIyT6ZTCaZTKbT6fKbnTt3rp7WZ\/oRTsJGEACqKW1ASDZJJBKRSCQWi3V0dAwODpba7OjRo21tbbZbn9VUcRI2gmBQTWkDQrLD8vJyc3Pz9PS0PJ+bm2tqapqZmSnYZn5+fu\/evbJZPUJi+lHYuOeZ0ZCsdycqYmWgUFHN6TVCskMqlZJgZL3s7e0dGhoq2GZgYODgwYMjIyP1CInpRyFEQpJEpcDbiGwUQiqO2iEkOwwPD3d3d1sv+\/r6+vv7C7bJ5XLydWxsrLyQLIqeUFDOEFonBfieftz7NbSUmsFSvidESBWIx+M9PT3Wy32rFN2yopDKv3mUM4TcSTNzn2MjCAwV11pFSHZIJBJdXV35CWlgYMBxITH9CCcFzElio9vpQ7yzIQ9JZc6zEZIdxsfH8zUjchJFOSukT448zXgdBMlJc2\/et5iJ856GnJunXyjTsyEkO+RyOdGMyEaeT01NNTY2zs7OyvOJiYlsNuuIkCTYUs4AwqvvZsVJv\/\/zp+b+CUuzfxQbyVfeTVgpW9qAkOyHpEgk0tnZ2dLSMjIyor4ZjUbj8Xj9QmL6ERRg7tLg2AgKKDNqh5D8pFTrM\/0I1mLitFlsBDWdcCMkHYVEOQMEICfdTh+aT7ZhI1hLqVE7hKSdkJh+BAFw0meXY2Ij3i8oSpkJSQhJLyFx9yMoT\/T1Sc2XcmAhBihPqTskISTthMR4HVRE53Jwpr5CNRS9QxJC0ktIjNeB0U6ae\/M+pr5CNUhHt\/YOSQhJLyFx9yOoHjVFSZNLShTUQU0UXUYIIWkkJHmHGK+DWtGhHHwxE6egDmpl7YQkhKSRkJh+BPbwNyep8u7crb\/wRkBNrF1GCCFpJCSWCwLbRF+fvOeZUe9\/79yb91FQB7YpGLVDSLoIiXIGqBPvyxwoqIM6KRi1Q0i6CInxOnDKSR4M36kSBlbvhnrPae5eRggh6SIkyhnAKdx2EmsCgYPkj9ohJC2ExHgdOIt7qzmwCgM4S\/7aNAhJCyGxXBA4jhuXlLhoBI6Tv4wQQtJCSIzXgXtOcmT4jotG4B7WMkIIyX8hMV4HrlK\/k7hoBK5iLSOEkPwXEssFgduoS0r2hu+4aARuYy0jhJB8FhLLBYE3qOG7V9\/NVv9f1DAdi6WCB6gJSQjJZyEx\/Qi8pPqF7ximAy9RywghJJ+FxHJB4DHR1ycrVt8xTAfeI50hQnIzhGYyyWQynU6X2qD7gfspZwDvKVN9xzAd+NZh\/vKnDzd\/FyG5QiKRiEQisViso6NjcHCw6DYvP\/htxuvAL9Y6iWE68BHpDM\/8+7cQkvMsLy83NzdPT0\/L87m5uaamppmZmbWbUc4AOkQlNXzHMB34zvn\/+DeE5DypVEqCkfWyt7d3aGioYJtb71\/4fz3f+f6XcCyCX07q+z8HGaYDH7G6QUlImg8aGSmk4eHh7u5u62VfX19\/f\/\/azd7\/zf1yTsqdzcBH5AicT7Y9ceS47YlKAI4gR+DXut4gITlPPB7v6emxXu5bpeh5gRq15+QUvGcxE88PRp7dugKggFffzd7zzKgce1TZuUIikejq6spPSAMDA0WFpJ6Ikz67HCMqgcfBaG39gjiJqAQeByM56n7\/509XWKnBJcbHx9va2qyXIidRVBkhrXxZ4ERUAo+D0VqISuBxMCraJSIkx8jlciKksbExeT41NdXY2Dg7O1teSEQl8DEYrUWcFH19kqgEHgQjhORFSIpEIp2dnS0tLSMjI0W3Kdr6RCVwAxszXolK4AYioYJghJC0oEzri5MowAPvg1GRT+Azo0QlcDUYISTdhURUAkeQ46f+OUZEJagfOX5KBSOEZICQ8rXEnTqhViQPOZuzlZZKndsClEKOmSqrNxGS7kKyehaKHaB65GhxY1U61bMwggfVo1aXr\/I8BiEZIKT8qMQIHlQ8TtxeB4gRPKjyOKk4RoeQTBWSq2e+EAy8rIWp6cwXwoa9GdYIyTAhrXw5sZERPCg4U5GjwuNrjWpiIyN4UHCmIkeFHBuedYkIyU8hKawRPLQUcnwfy7VG8NBSyKl\/LBchmSokTToj8FdF+tw2ggtLIVdRrZeLEFIAhYSWwqwiDadOoyVUhJDCLiSFqndAS8FGzQGQh85DtWo2PloKg4qcXRseIQVHSIp6FokBzVVk1ptbfpEYMBc1Hc2NNxchBU1IK6tleKrnYn2HgKUi497QV9\/NKi3Zq7mCsJ1nIKQACgktBQZ1ddD0NxEtBSkVuToSi5ACKyQLpSWuLRmqoiCNvnJtyURUoYo35xMIKfhCUlDyYJCK9Kygcwq1xANaMkJF9zwz6uWFQIQUFiFZWlIzV5hOq62KQrIGh5rPz3RanVXkbAUdQkJIFTo+Li9pgpVfw3aiYHV8XF7SLb\/6cqKAkMIoJIVV9cA4no+RSF0oCvlbYFU9EJh8PDPwoGYBISGkak\/P5SuzlzwOqYydFu0W5SSd2Useh1RNTgUQkn0ymUwymUyn0xW3PHfunP6tr87W6SU9a2Rag8Dku\/sdXPUHIflJIpGIRCKxWKyjo2NwcLDMlkePHm1razOo9VXhAxOYnPWQNTpKDK0edT2DK0zOesiSvYYxFCHZYXl5ubm5eXp6Wp7Pzc01NTXNzMys3Wx+fn7v3r2ypVlCUoiNrKE8zGTbQ1a1Am1YT2BSZpKvmMm2h6yyRp3bECHZIZVKSTCyXvb29g4NDa3dbGBg4ODBgyMjIyYKKd9M6uweM1WJZCB1iYhqBcfNpM7uMVOVSAay5hIZMfELIdlheHi4u7vbetnX19ff3792s1wuJ1\/HxsbKCMnCIDPJV4aeiuYhaRk85KWZ5KueQ0++5yFpGbM8ZEpPqKOQ4vF4T0+P9XLfKqU2Li8kEw93NZpHbLI8ZJVukyC9N5MazaM2ryAMGdoUCKla9u\/f37yK2CWRSHR1deUnpIGBgfAIKeSxSY3I5YchwqJusUmeh0ROBWHI9KJEhFQtV69evbDKxYsXx8fH8x0jchJFhVNI+d20FZvka8CSk\/wtSMis2KTkFLwLTvK35EtIvgZJvQjJDrlcThwjppHnU1NTjY2Ns7Oz6p8mJiay2awpQnJpB6SzVslJdd9KTkV7cB2Ov6L7oGKQVQHvav7TthECsAPSWSs\/qe5byaloD67tu6DG4qw\/wdXLZkE9DAIuJEFCUiQS6ezsbGlpGRkZsb4fjUbj8XjIhbQ2XqjwpDp3S1HyeOyh7+nQCyj95GcgFYO8qU1ASJ7tgBWerM7dUpQm74LST34GUmNx3tQmICRThcTbX2d+UopKvnS\/pSh5qLoA2cClLKJ+teUepZ\/3f3O\/CkB+DcQhJL92wMpP0vt\/resNa4hPHmq+jkvXotTvVZpRv13pR\/ZBBSC\/qhIQUtiFBIpd\/\/k9efzvnu\/+3\/7\/IZaSh3hC2UK+qu+oh2xgPWR762F9M39j+b\/qIT9HHuq\/q+3l10k+o+WhgH998Gfy+O9bd\/23\/\/mrr\/\/sv+QhnlC2kK\/qO+ohG1gP2d56WN\/M31j+r3rIz5GHfEdtLL\/rW61b5UHLU\/YNZqAyk8o3+Q8Va9TD+qa1MXUH4Ea4sR4qQqmHijXqYX0zf2OaLhggJAAAQEgAAAAICQAAEBIAAABCAgAAhAQAAICQAAAAIQEAACCkc+fOefa7MplMMplMp9OlNpiZmZENJicnfdyHubm50dHR8fFxv3ZAMTExcePGDb\/2YXp6Wja4dOmS74dEIP92348B3z8FmvQG\/vaHCKmQo0ePllp61XESiUQkEonFYh0dHYODg2s3OH78uNpg8+bNZW456Oo+jI2Ntba27t69e8eOHU888YS6666XO2D1ievXr5dPoy+NsH\/\/fvkn2WDbtm2PP\/744uKiX4dEIP92348B3z8FmvQG\/vaHCOku5ufn9+7dq+7458GvW15elt8lnzF18tXU1CSnP\/kbyEHf0NAwNTUlzxcWFuS542dGFfdBNpDP4cWLF9XLLVu25K+k7sEOKJaWlqQ3bG9vd6MzqrgP0uzSD8rhoV5u3bq1YP14j1sjYH+778eA758CTXoDf\/tDhFTIwMDAwYMH5VDz5g1IpVJyKmS97O3tHRoaKjgE161bJylefRqlX5iYmPB4H0ZHR+WU0MdGUBw4cODw4cM7d+50Q0gV9yGbzV64cCF\/A9kZH1sjYH+778eA758CTXoDf\/tDhFSIiuFl7pbkLMPDw93d3dbLvr6+\/v7+gm1OnTol56SS3x955JH9+\/d7vw9nzpzZtWuXfLOxsVHOlY4dO+Z9I8iZ6fbt2+WJS0KqZh8srl27Jn2BSyenNe1JYP52348B3z8FmvQG\/vaHCKk4nr0B8Xi8p6fHerlvlYJtYrGYHHwnTpyQz2E0Gr1165bH+yCnpQ0NDfJJkOfpdHrDhg3nz5\/3cgcWFhY2b96shi9cElI1b4Ti+vXr7e3tL7\/8so+HRPD+dt+PAd8\/BZr0Bv72hwjpztXa5lXyW9zVNyD\/NyYSia6urvxzIsnIBQMFmzZtWl5eVi\/lEDxy5IjH+3Dy5MmHHnoofwPByx2Qj6WcnI6tsmPHDjk9dKT8rKZ9UFy5cqW1tfX48ePuHZBV7on3v9HVv92vY6D6HXDjU1DrPrjUGyAkjbh69eqFVazLlW6\/Afm\/cXx8PP8XyeEoB2XBQEF+ipcDdM+ePR7vw9mzZ\/M\/io6cs9e0A9L77PwS6RO3b9\/uSLdY0z4IsqWcF0truHpAVrMn3v9Gt\/92v46B6nfAjU9BrfvgUm+AkBiy+4JcLie\/SH6dPJ+ammpsbJydnV1ZnWmRzWZXVguc5JvSdapRiy1btshB6fE+LC0tPfDAA6lUamW1+Gfjxo3OzsOouAP5uDRkV3EfMpmMZClphKUvsU5UvdkTHw9CD\/52348B3z8FmvQGCCnUQlKnRZFIpLOzs6WlxSoklTBuVdaeOnVK\/klt8Pzzz\/uyD++99157e\/ujjz4qGxw9etT7HXBbSBX34cCBAwU3dX7uuee8PCR8PAi9+dt9PwZ8\/xRo0hsgJFi5detWmXl28k+ygUun5FXug3D79m1X96HiDvj+RgR7T2h\/HT4FmvQGZoGQAAAAIQEAACAkAABASAAAAAgJAAAQEgAAAEICAACEBAAAgJAAfOf999+nEQAQEoD\/DAwMjI2N3bhxg6YAQEgA\/nDp0qUnn3zylVdeee6550ZHR2kQAIQE4APZbLa5udmlO9ICICQAqJbf\/va3+bfeAQCEBOAPly9f3rBhw+zs7MWLF1nUGQAhAfjG\/Pz87t279+zZs2\/fvt\/97nc0CABCAvCBpaWlvr6+xcVFmgIAIQH4ydDQkLqJNQAgJAA\/icVib731Fu0AgJAAfObChQuNjY0HDhy4evUqrQGAkAD8JJ1O9\/X1rV+\/fseOHZcvX6ZBABASgJ8sLi6+9NJLv\/jFL2gKAIQE4BtXrlw5efLk4cOHx8fHaQ0AhATgM7lcjkYAQEgAAICQAAAAEBIAAJjF\/wfyW4CHlVbpGwAAAABJRU5ErkJggg==","width":386}
%---
