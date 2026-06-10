%[text] %[text:anchor:T_809F8C6C] # Plot Bilinear Lagrange Polynomials in Physical Space
%[text] Generates three graphs concerning a two-dimensional bilinear mesh in a tiled format: The two-dimensional basis functions, the provided parametric location mapped onto the physical space and the distribution of the geometric Jacobian in the parametric space of the provided element.
%[text:tableOfContents]{"heading":"**Table of Contents**"}
function t = plot2dBilinearLagrangePolynomialsOnPhysicalSpace ...
    (el, msh, numVertices, l1L1, l2L1, ij)
%[text] %[text:anchor:H_266A335B] ## **Function description**
%[text] Displays the bilinear basis functions on the two-dimensional bilinear quadrilateral mesh in the physical space.
%[text]  **Input :**
%[text]  `el` : The element for which the parametric location is sought to be displayed in the physical space
%[text]  `msh` : The one-dimensional bilinear mesh using four-noded elements as a structure with the following fields:
%[text]                               .`nodes` : Number of nodes in the one-dimensional mesh
%[text]                         .`elements` : Number of elements in the one-dimensional mesh
%[text]  `numVertices` : Number of vertices in the bilinear quadratic mesh
%[text]  `l1_L1`, `l2_L1` : The two Lagrange polynomials for the one-dimensional linear element
%[text]{"align":"center"} ![linear_lagange_polynomials.png](text:image:24c1)
%[text]  `ij` : Legend entries for the plot which should be an array with four elements
%[text] 
%[text]  **Output :**
%[text]  `t` : Handle to the tiledlayout graphical object
%[text] %[text:anchor:H_1760F562] ## **Function Implementation**
%[text] %[text:anchor:H_D68978A7] ### Input validation
    arguments
        el (1, 1) double {mustBeInteger, mustBeGreaterThanOrEqual(el, 1)}
        msh (1, 1) {mustHaveNodesAndElements}
        numVertices (1, 1) double {mustBeInteger, mustBePositive}
        l1L1 (1, 1) symfun
        l2L1 (1, 1) symfun
        ij (1, 4)
    end
%[text] %[text:anchor:H_423B2276] ### Read input
    syms l_Q1(x, y, x1, x2, y1, y2)
%[text] %[text:anchor:H_2875BE80] ### Initialization of the tiled layout
    if floor(el) ~= el
        error("Variable 'el' should be an integer");
    end
    if el < 1 || el > numel(msh.elements(:, 1))
        error("Variable 'el' is chosen as %i but it should be an integer in the close interval [1, %i]", el, numel(msh.elements(:, 1)));
    end
    t = tiledlayout(2, 2);
    numNodes = numel(msh.elements(el, :));
%[text] %[text:anchor:H_F0E23FCC] ### Nodes of the selected element in the quadrilateral patch
    Xel_Q1 =  msh.nodes(msh.elements(el, :), :); %[text:anchor:H_6203C08B]
%[text] %[text:anchor:H_8D2ACC7C] ### Loop over all the nodes in the quadrilateral patch
    for ii = 1:numNodes
%[text] %[text:anchor:H_B42A92CB] #### Lagrange polynomial for the current node
        switch ii
            case 1
                idx1 = 1;
                idx2 = 4;
                idy1 = 1;
                idy2 = 2;
                l_Q1(x, y) = l1L1(x, Xel_Q1(idx1, 1), Xel_Q1(idx2, 1))*l1L1(y, Xel_Q1(idy1, 2), Xel_Q1(idy2, 2));
            case 2
                idx1 = 2;
                idx2 = 3;
                idy1 = 1;
                idy2 = 2;
                l_Q1(x, y) = l1L1(x, Xel_Q1(idx1, 1), Xel_Q1(idx2, 1))*l2L1(y, Xel_Q1(idy1, 2), Xel_Q1(idy2, 2));
            case 3
                idx1 = 2;
                idx2 = 3;
                idy1 = 4;
                idy2 = 3;
                l_Q1(x, y) = l2L1(x, Xel_Q1(idx1, 1), Xel_Q1(idx2, 1))*l2L1(y, Xel_Q1(idy1, 2), Xel_Q1(idy2, 2));
            case 4
                idx1 = 1;
                idx2 = 4;
                idy1 = 4;
                idy2 = 3;
                l_Q1(x, y) = l2L1(x, Xel_Q1(idx1, 1), Xel_Q1(idx2, 1))*l1L1(y, Xel_Q1(idy1, 2), Xel_Q1(idy2, 2));
            otherwise
                error("The bilinear quadrilateral should only have 4 nodes");
        end
%[text] %[text:anchor:H_E4F400ED] #### graphical representation of the mesh along with the corresponding Lagrange polynomial associated to the current node
        nexttile
        plotMesh2d(msh, numVertices, true);
        hold on;
        fsurf(l_Q1(x, y), [min(Xel_Q1([idx1 idx2], 1)) max(Xel_Q1([idx1 idx2], 1)) min(Xel_Q1([idy1 idy2], 2)) max(Xel_Q1([idy1 idy2], 2))], ...
            "FaceColor",[217, 218, 219]./255, "FaceAlpha", 0.5);
        xlabel("x");
        ylabel("y");
        zlabel(sprintf("l_{%s}", num2str(ij(ii))))
        axis equal;
        hold off;
    end
    title(t, sprintf("Bilinear basis functions in the physical space of element %i", el));
%[text] 
end

%[appendix]{"version":"1.0"}
%---
%[text:image:24c1]
%   data: {"align":"baseline","height":255,"src":"data:image\/png;base64,iVBORw0KGgoAAAANSUhEUgAAAjAAAAGkCAIAAACgjIjwAAAlt0lEQVR42u2dX2hbV76o5+U8zNM8nKfeC3NyXgbiMRLCoRkLEeybGUyvcwNpQ0pJ8RVDCTa+nuATYSKMqwkJxaTTImPc04fgTDGJ7VEwRUPRDcaqSPCgpCSIMHVVqRkFwai+dlRjhsw1RsJn1btVVf+RJXn\/WWvv7yOEve3VZvm3916f1vot799PtgAAACTgJ4QAAAAQEgAAAEICAACEBAAAgJAAAAAhAQAAICQAAEBIAAAACAkAABASAAAAQgIAAIQEAACAkAAAACEBAAAgJAAAQEgAAAAICQAAEBIAAABCAgAAhAQAAICQAAAAIQEAACAkAABASAAAAAgJAAAQEgAAAEKqj3v37nEVAQAQksVMTEz4fD6uIgAAQrKMtbW1y5cvezwehAQAgJCsJBQKXb9+PRaL1RDSp199wwUGAEBIxlIul8XfiUSihpB++r\/\/819PDf4CAMDx\/N73b\/f\/139HSAZSW0jiGhy5tvj7\/\/s3q7onOmBtfCzvAEEgCFwFGTqQf\/vVv\/W+LMNVcLSQxN9CSD\/5j4Vc8f\/zADAKEASuggM7IGz09fhFSa6C04VUcZL5UyVGAYJAELgKFnbgxV8Xv3ztpfX4rDxXASF9h7XLdwAAZvJ89g9\/631ZOEmqjwUI6Qf800sdE4+4UwHA3uTfflX8kXCeamchNTFBtjClBABgAmJipCWNEJLsQtqyLqUEAGAoWtKoepkOIckuJA1SSgBgJ3YnjRCSMkLaYvkOAOzCnkkjhKSSkHASANiAL197SUyPdBkSEZKVQtrafuUdKSUAUJH1+GyNpBFCUk9IGqSUAEAtvh6\/WDtphJBUFdIWy3cAoA71JI0QksJCqjjp5sMCtzsAyIm2t7uepBFCUltI3wWFlBIASInwUP1JI4RkByFtkVICAPnQqkhs\/r+8+UMiQrJSSFuklABAMhvt+UIghOQIIW3xkiEAkICmk0YIyVZC0uiYeISTAMASDnwhEEJylpC2tlNK1K0AAJNpYm83QrK\/kLZIKQGA6TY6ZNIIIdlWSFuklADAFHaUHkdICGlf2BEOAMahY9IIIdlfSFvbpdBZvgPlHgflcOBl0jdphJAcIaQtUkrA40CH9UZMjA6\/txshOfSGvvmwQEoJeBzo8OFpqIoEQuIJ3BdSSsDjQIcPg3FJI4TkxBua5TvgcaDDzWFo0gghOfSGxknA44CQGkKvFwIhJFM+OOTz8\/Pz6XR6vwbZbFY0ePbsmSTRpxQ6ML7XyRdffPHee+85WUjr8VlzlukQkg5Eo1Gv1xsIBDo7O8Ph8O4G169fb29vFw26urrGx8fliT5OAoR0IOKj5K9\/\/WvHCunr8YtibtR0FQmEZCqlUsnj8YgJkDguFotutzuXy1U3ePLkSWtra6HwbV3XjY0NIS3xFXmiz\/IdICSEtB9mJo0Qkg7E43HhmMrpwMDA1NRUdYO5ubm+vr7KqZgnvfvuu1JFn1LogJAQ0g7MTxohJB3Y4ZtgMDgyMlLd4JNPPjl9+nTl9MKFC0NDQ3tG39pf+aZuBSAkhKRh2t5uCUdCtYUUiUT6+\/srp8PbVDdYW1s7fvz49evXHzx48Mc\/\/lHLNsl5Qx+5tuifXmI0BITkZCFZuEzHDOmwRKPR3t7e6hlSKBTa0SabzQppnT9\/fnx8fHR0dIexpIo+KSWQTUjiQ5Ilf5wpJN2rSCAkU0kmkz6fr3Iq5CQUVd3gH\/\/4RyqVqm5w584dmaNP3QqQajC6+bBgyZ8aQtI2MdlMSAZVkUBIplIul4WQEomEOM5kMi6Xa3V1VRwLCWk768TfLS0ty8vL4vjx48fHjh1bX1+XP\/q8ZAgYjPYU0p\/+9KczZ87YTEjPZ\/9gwuvpEJJJkySv19vT09PW1haLxb5bZ\/D7I5GIdvzRRx95PJ4333yzo6NDNFYl+v7pJUqhA4NRNRcvXuzu7raZkCRJGiEkon8ApJSAx6Gav\/zlL8+fP7eTkIyuIoGQeAL1dxLLd8DjoGEbIWlJI6mW6RAS0a8LUkrA42AnIVn4m0YIiSdQt6kSy3eAkFQXkpxJI4RE9HES8Dg4SEiWvxAIIfEE6gl1K4DHQdEOW1JFAiHxBBoOKSXgcVCrw1+PXxQ2Mr+KBELiCTQDlu+Ax0GVDiuRNEJIRF8HJ1G3AngcpO2wQkkjhET0dYC6FcDjIGeH5d\/bjZCIvv6QUgKEJFuX8m+\/KmzEPYCQHCekLVJKgJAks5EkVSQQEtG30klMlcBmj8MXX3zx3nvvqdJh2apIICSibyUs34HNHgeFCvQpnTRCSETfEKhbAQjJfMTESLm93QiJ6JsBKSVASGYiJkaqJ40QEtE33Eks3wFCMhTJq0ggJKIvEaSUACEZh52SRgiJ6Js3VWL5Dpp7HNbjs5b8kV9IKr4QCCERfZwECj8OX49ftOTPfkL6\/PPPX7x4Yfnzq+gLgRAS0ZcC6laAoo9DRUjPnz\/3er1+v7+jo+Pdd9+1qsNi9mbLpBFCIvpmQ0oJ1H0cbty4MT4+rh3\/7Gc\/s6TDWhUJ29sIIRF9k2D5DlR\/HAqFws9\/\/nPzO2zvpBFCIvpWOom6FaDi4\/D8+fNf\/epX8XjczA6rW0UCIRF9NaBuBSj3ODx9+vTll1+en583s8M23tuNkIyZSufz4h5Np9P7NcjlcqLB0tISQqqGlBIo9DiIx\/yXv\/zl559\/bmaHtSoSCpUeR0gWE41GvV5vIBDo7OwMh8O7G0xOTmoNurq6hoeHEVI1pJRAlcehr6\/vpz\/96X\/7HhM6bIMqEgjJVEqlksfjyWaz4rhYLLrdbjEZqm5QLpdbWloymYw4Xl9fF8d7zpMcK6QtXjIE6j8OunfYgUkjhKQD8XhcTIwqpwMDA1NTUzuEdPToUTHZF8ebm5utra2pVGrP6Fdw5v1HSgkQkoYzk0bKjYQyCmlubk5M5CunwWBwZGRkR5uZmZlTp06Fw+Fz585dvXrVHk+gERy5tkjdCnC4kBy1t5sZks5EIpH+\/v7K6fA2O9oEAgGholu3bl24cMHv9+9+7whCqkBKCZwsJCcnjRCSDkSj0d7e3uoZUigUqm6wsLBw8uTJUqmknQohVX7ZGyHVcBLLdw4fjJTjkD+ybUqPIyQrSSaTPp+vcirkJBRV3eDOnTvVa3pCV0NDQwjpQNgRDs7ByUkjhKQn5XJZCCmRSIjjTCbjcrlWV1fFcSqVKhS+fRPB0tKS+OLTp0+3tnfZdXd3C0UhpHrwTy+xfAe2h6QRQtJ5kuT1ent6etra2mKx2HeDqd8fiUS045mZGfEtrcE777yjaPQtgZQS2BsxMXLy3m6ERPQV4+bDAiklsB8OqSKBkIi+DSGlBHaCpBFCIvpqw\/Id2AOSRgiJ6NvHSZ9+9Q2hABXhhUAIiejbCkqhg6Ksx2dZpkNIRN+GkFICtdBKjzuwigRCIvqOgJQSqAJJI4RE9J3iJEqhg7SQNEJIRN9ZULcC5IS93QiJ6DuRI9cW\/dNLxAHkgWU6hET0nQspJZCHL197iSoSCIno4yR2hIOVUEUCIRF9+AF2hINVkDRCSEQfduKfXqIUOpgMSSOERPRhb0gpgZlQRQIhEX042Eks34GhaEkjlukQEtGHgyGlBMZB0gghEX1oZqrE8h3oC0kjhET0ASeBxfBCIIZEhASHgroVoAtUkWBIREigD6SU4DBQRYIhESGBnrB8B81B0oghESGBUU6ibgXUCUkjhkSEBMZC3QqoB\/Z2MyQqJqR8Pj8\/P59Op3d\/q1gsfvZjcrkcQpIEUkpwwKP99qvCRsQBISkjpGg06vV6A4FAZ2dnOBze8d27d+96qjh69GgoFEJI8kBKCWrYiCoSCEklIZVKJaGZbDarTYbcbveeEyCN+\/fvnzhxYm1tDSFJ6CSmSlCBKhIISUkhxeNxMTGqnA4MDExNTe19i7944fP57t27t1\/0K3AvWgLLd6BB0shaD6kyEsoopLm5ub6+vsppMBgcGRnZs2U4HH7rrbfU\/TjgBKhbAeztZoaksJAikUh\/f3\/ldHib3c02NjbcbveTJ08QkuSQUnIyVJFASGoLKRqN9vb2Vs+Q9tyz8PHHH586dUrp6DvNSSzfOQqqSCAkOwgpmUz6fL7KqZCTUNTuZoODg2NjYwhJIUgpOQeSRgjJJkIql8tCSIlEQhxnMhmXy7W6uiqOU6lUofDDiwDa29u1NghJuakSy3f2hqQRQrKPkLRJktfr7enpaWtri8Vi2hf9fn8kEqlISwR3ZWUFIeEkkApeCISQ7CYkh0TfsVC3wpasx2dJGjEkIiSEpCSklOyEVkUCGzEkIiSEpCos39kDkkYMiQgJIdnHSdStUBSqSDAkIiSEZDeoW6Ei7O1mSERICMmekFJSC62KBKXHGRIREkKyJ6SUFLIRVSQYEhESQnKEk5gqSQtJI4ZEhISQnAUpJTkhacSQiJAQkhM5cm2RuhVSwd5uhkSEhJCcCyklqWxE0oghESEhJJxESslKKD3OkIiQEBL8ADvCrYKkEUMiQkJIsBP\/9BLLdyZD0oghESEhJNgbUkpmQulxhkSEhJCgFjcfFli+MxqqSDAkIiSEBPWCk4yDpBFDIkJCSNAYLN8ZAUkjhISQEBLgJIvhhUAICSEhJDgUlELXhfX4LMt0CAkhISTQAVJKh0ErPU4VCYSEkBAS6APLd81B0gghISSiD0Y5iVLodULSCCEhJKIPxkLdinpgbzdCQkgNLibk8\/Pz8+l0er8GxWJxYWEhmUwiJKjmyLVF\/\/QScdj3yWKZDiEhpIaIRqNerzcQCHR2dobD4d0NEolEe3v7pUuXzp49e\/78+XK5jJCgAimlGjaiigRCQkgNUCqVPB5PNpvVpkFutzuXy+1oIGz04MED7bS7uzsWiyEk2O0klu8qUEUCEFIzxONxMTGqnA4MDExNTVU3WFhYEBMjG0QfjIYd4RokjQAhNcnc3FxfX1\/lNBgMjoyMVDe4c+fO4OCg+KLL5RJzqRs3buwX\/Qrci47FP73k8FLoJI3wkCojoYxCikQi\/f39ldPhbaobjI6OtrS0zMzMiON0On3s2LH79+8zQ4L9cHJKiSoSwAzpUESj0d7e3uoZUigUqm5w+\/btV155pbqBACHBgU5y1PKdljRimQ4Q0qFIJpM+n69yKuQkFFXd4O7du9VC2j2FQkiw9+3uGCeJWRE2AoSkA+VyWQgpkUiI40wm43K5VldXxXEqlSoUvv09\/M3NzePHj8fj8a3tbXgnTpzY87eREBLsN1Wy9\/IdSSNASDpPkrxeb09PT1tbW2VLt9\/vj0Qi2vFnn33W0dHx+uuviwYTExOKRh9wkr7wQiBASEQfFMOWdSuoIgEIieiDqtjpt5SoIgEIieiD2thj+Y6kESAkog\/2cZKidStIGgFCIvpgN1SsW8ELgQAhEX2wJ2ptc8i\/\/aqwEVcNEBLRB3uiSkqJKhKAkIg+OMVJ0k6VqCIBCInog7OQc0c4SSNASEQfnIhsdSvY2w0IieiDc5EnpSQmRiSNACERfcBJVqaUqCIBCInoA\/yAVSklkkaAkIg+wN5TJTOX70gaAUIi+gDWO4kXAgFCIvoAtTChbsV6fJakESAkog9QF8allLQqEtgIEBLRB6gXI5bvSBoBQiL6AM07SZe6FVSRAIZEhARwWA5ft4K93cCQiJAA9OEwKSWtigSlx4EhESEB6ENzKSWqSABDIkICMMpJdU6VSBoBQyJCAjCWelJKJI2AIREhAZjBkWuLNepWsLcbGBLVE1I+n5+fn0+n03t+t1gsflbF+vo6QgJ52C+lRNIIEJLFQvrtNoVC4fiPaW9v3+8\/iUajXq83EAh0dnaGw+HdDW7cuNHS0uL5nvv37yMkkNBJleU7So8DQpJCSCvbVA6q2bN9qVQSjslms9pMyO1253K5HW0GBwdv3bqlevTB9mg7wkkaAUKSRUiNEo\/HxcSocjowMDA1NbWjTVdXVzKZFLra3NysEf0K3ItgFff\/z\/+cOv8\/ZCg7C471kCojoW5CenYQ9f+v5ubm+vr6KqfBYHBkZGTHFOro0aPd3d3t7e3iYHh4mBkSyImYGInpkTyl0IEZkiOEdPLkyV\/sj9BG\/f+rSCTS399fOR3eprrB3\/\/+dzFtEn+L4+Xl5RMnTty+fRshgVTsqCJx82HB2lLoAA4SUukg6v9fRaPR3t7e6hlSKBSq0f7q1auXLl1CSCAP+yWNrCqFDrBFDqk5ksmkz+ernAo5CUVVN3j27JmYRVVOR0ZGhoaGEBJIQu3fNGL5DhCSSkIql8tCSIlEQhxnMhmXy7W6uiqOU6lUofDtq\/7T6XRLS4u2DW95ednr9bLtG2SgzhcCaU769KtviBggJNmFpE2ShGZ6enra2tpisZj2Rb\/fX5kY3bp1y+PxiAbi78nJSUWjD3ZiPT5b\/95uE0qhAyAkog9ORCs93mgVCVJKwJCIkAD05DCvpyOlBAyJCAlAB3SpIqFjKXQAhET0wYno+0Kgw5dCB0BIRB+ciBFVJI5cW\/RPLxFbQEhEH6ABGxlURYKUEiAkog9QFyZUkWioFDoAQiL64ETMrCLBjnBASEQfYG\/MLz3un16qUQodACERfXAiWhUJ8\/9dUkqAkIg+wHdoSSMLi72SUgKERPQBtuQpPU5KCRAS0QfnYn7SqJ6pEst3gJCIPjgIXV4IhJOAIREhISQ4FA1VkTAf6lYAQiL64AiaqyJhPqSUACERfbAzX772klRJo9qwfAcIieiDDZE2aVSPk6hbAQiJ6INNkGdvd3NQtwIQEtEHO5B\/+1VhI9V\/ClJKgJCIPihvI4OqSJgPKSVASEQflMSEKhJWOYmpEiAkog\/KoHrSqDYs3wFCIvqgBrK9EMgIqFsBCEnvgSOfn5+fT6fTtZulUqmVlRWEBPXw5Wsv2SZpVBtSSoCQdCMajXq93kAg0NnZGQ6H92uWzWZbW1uFtxAS1MbyKhJWOYnlO0BIh6JUKnk8HiEbcVwsFt1udy6X291sc3Pz9OnTHR0dCAlqY++kUW1IKQFCOhTxeFxMjCqnAwMDU1NTu5uNjo6OjY1duHChhpAqcC86FickjeqZKrF852QPqTISyiikubm5vr6+ymkwGBwZGdnR5sGDB2fOnBEHtYXEvehwVHwhEE4CZkgSEYlE+vv7K6fD21Q3WF9f7+rq0tbxEBLsyXp81mlJo9pQtwIQUjNEo9He3t7qGVIoFKpuIPw0ODiY2Obs2bPhcHjPzXgIybFoVSSw0W5IKSEkhNQYyWTS5\/NVToWchKKqGwgDXfie9vb2M2fOTE5OIiTQIGlUG5bvEBJCaoByuSyEJGY\/4jiTybhcrtXV1a3tXzkqFHa+aZ8lO6igaBUJq5xE3QqEhJDqnSR5vd6enp62trZYLKZ90e\/3RyIRhAR74uS93c1B3QqEhJCIPuiPVkVC\/tLjskFKCSEhJKIPOtvIIS8EMgJSSggJIRF90AGSRjo6iakSQkJICAmahKSRvpBSQkgICSFBM7C32wiOXFukbgVDIkJCSNCYjUgaGQQpJYZEhISQoC5sWXpcTiexfMeQiJAQEuwLSSMzYUc4QyJCQkiwNySNzMc\/vcTyHUMiQkJI8CPExIi93ZZASokhESEhJPgOqkhYzs2HBVJKDIkICSE5HZJG8kBKiSERISEk50LSSDZYvmNIREgIyXHwQiCcBAiJ6IP1rMdnWaaTGUqhMyQiJITkCLTS41SRkB9SSgyJCAkh2RmSRmrB8h1DIkJCSDaEpJHSTqIUOkMiQkJINoG93apD3QqGRISEkOwAy3T24Mi1Rf\/0EnFgSERICElhG1FFwjaQUmJIREgISUmoImFjJ7F8x5CIkBCSMpA0sjfsCGdIREgISQ1IGjkB\/\/QSpdAZEm0lpHw+Pz8\/n06n92sgviUa5HI5hKQKVJFwDqSUEJJ9hBSNRr1ebyAQ6OzsDIfDuxu8\/\/77v\/nNby5fvnzy5MkPP\/wQIUmOljRimc6BTmL5DiGpLaRSqeTxeLLZrDguFotut3vHNCiTybS2tq6trYnjlZWVo0ePimYISVpIGjkZUkoISW0hxeNxMTGqnA4MDExNTVU3KJfLmq4EQksiysvLywhJTkgaAct3CElhIc3NzfX19VVOg8HgyMjInhOpmZmZ06dPj42N7Rf9CtyL5sMLgQAnSeIhVUZCGYUUiUT6+\/srp8Pb7G62srLy0UcfvfXWW+fOndOW75ghyQNVJGAHn371Dct3zJDUE1I0Gu3t7a2eIYVCoRrte3p69tz4gJCsgioSsB84CSEpJqRkMunz+SqnQk5CUdUNnj59Wp1VGhoaunz5MkKSBJJGUBuW7xCSSkIql8tCSIlEYmt7Q53L5VpdXRXHqVSqUChoX2xpaRFaEsfiW16vd2FhASFZDkkjaMhJ1K1ASAoISZskCc309PS0tbXFYjHti36\/PxKJaMe3b992u92\/\/e1vxd\/8HpIMsLcbGoW6FQhJDSE5JPq2If\/2q8JGxAEahZQSQyJCAp1tRBUJaBpSSgyJCAl0gCoSoKOTmCoxJCIkaJLns3\/g9XSg52CEkxAS0YcmYG83GAF1KxAS0YfGoIoEGAcpJYRE9KEuqCIBpjmJ5TuERPRhX\/hNIzATdoQjJKIPe0PSCKyaKrF8h5CIPvwALwQCnMSQiJAQksWsx2dJGoG1fPrVN6SUEBLRdzpfj1\/ERiDLUIWTEBLRdywkjUA2WL5DSETfcVBFAiR3EnUrEBLRdwTs7Qb5oW4FQiL69kerIkHpcZAffksJIRF9m9uIKhKgEKSUEBLRtyEkjUBpJzFVQkhE3yaQNALVIaWEkIi+HWBvN9iDI9cWqVuBkIi+2jYiaQS2gZQSQiL6SkLpcbCxk1i+Q0hEXxlIGoG9YUc4QiL6akDSCJyAf3qJ5TuEpOvQmc\/Pz8+n0+n9GmSzWdHg0aNHCKlOKD0OzoGUEkLSjWg06vV6A4FAZ2dnOBze3eDq1aviW6LB6dOn33jjjY2NDYRUA6pIgAO5+bBASgkhHZZSqeTxeMQESBwXi0W3253L5aobLC0ttba2rq2taaenTp2KRCIIaT9IGoGTIaWEkA5FPB4Xs5\/K6cDAwNTUVHWDQqGwuLhY3WBsbAwh7QlJIwCW7xBS88zNzfX19VVOg8HgyMjIfo2fPXsmZktizrRn9Cs48ObjhUAAOEmtkVBGIUUikf7+\/srp8DZ7tlxeXu7o6Pjggw8U\/ThgHOvxWZbpAKqhFDpCaoZoNNrb21s9QwqFQrubPXnypL29fXJyUt3oG8TX4xepIgGwJ05OKSGkZkgmkz6fr3Iq5CQUtaPN4uLisWPH7t69q3T0jYCkEUBtHLt8h5CaoVwuCyElEglxnMlkXC7X6uqqOE6lUoXCt6WL8\/m8x+OJx+Ob31MqlRASSSOAhpzktFLoCKn5SZLX6+3p6Wlra4vFYtoX\/X6\/tr17dHT0Fz\/mypUrDhcSe7sBGsVpdSsQEtE3A5bpAJrjyLVF\/\/QSQkJIRF83G1FFAqBpnJNSQkhE30CoIgGgo5Nsv3yHkIi+UZA0AtAX2+8IR0hE3xBIGgEYgX96ycal0BES0dcfqkgAGIeNU0oIiejriZY0YpkOwAQn2W\/5DiERfd0gaQRgJvZLKSEkoq8PJI0ALBgf\/2PBTst3CInoHxZeCARgIXZKKSEkon8oqCIBYDm2qVuBkIh+81BFAkAebJBSQkhEv0lIGgHIhurLdwiJ6DcMSSMAyZ2kaN0KhET0G4O93QDyo2jdCoRE9Bsg\/\/arwkY87QDyo2JKCSER\/QZsRBUJAIVQLqWEkIj+wVBFAkBpJ6kyVUJIRP8ASBoBqI4qy3cIiejXgr3dAPZAiboVCIno7wtVJADshPwpJYRE9PeAKhIANnaStMt3CIno74SkEYC9kTalhJCI\/o8gaQTgnKmSbMt3CIno\/wAvBALASQyJSgopn8\/Pz8+n0+naze7duyd\/9NfjsySNAJyGbHUrEFKTRKNRr9cbCAQ6OzvD4fB+zSYmJnw+n+TR16pIYCMAZyJPSgkhNUOpVPJ4PNlsVhwXi0W3253L5Xa0WVtbu3z5smgmuZBIGgGAJMt3CKkZ4vG4mBhVTgcGBqampna0CYVC169fj8VitYVUwfyfgioSALDDSZbUrbB2JFReSHNzc319fZXTYDA4MjKyo025XBZ\/JxIJOWdI7O0GgN1YW7cCITVDJBLp7++vnA5vs2dLOYWkVZGg9DgA7MbClBJCaoZoNNrb21s9QwqFQqoI6cvXXqKKBADUwKqUEkJqhmQyWa0ZISehKPmFRNIIABpykslTJYTUDOVyWWhGyEYcZzIZl8u1uroqjlOpVKFQkFNIJI0AoFFMTikhpOYnSV6vt6enp62tLRaLaV\/0+\/2RSERCIbG3GwCadpJpdSsQkpWYE31KjwPAYTAtpYSQ7CwkSo8DgI5OMnr5DiHZVkgkjQBAX4zeEY6Q7CkkkkYAYAT+6SXjlu8Qkg2FROlxADAO41JKCMlWQqKKBACYwM2HBSNSSgjJPkIiaQQAZqJ7Sgkh2URIJI0AwHz0Xb5DSMoLiRcCAYA9nISQ1BbSenyWZToAsBa9SqEjJIWFpJUep4oEAMjA4VNKCElVIZE0AgDZOOTyHUJST0gkjQBAcic1VwodISkmJPZ2A4D8NFe3AiGpJCSW6QBAFY5cW\/RPLyEkewqJKhIAoBaNppQQkgJCoooEACjtpDqX7xCS7EIiaQQAqlPnjnCEJLWQSBoBgD3wTy8dWAodIckrJKpIAICdODClhJBkFJKWNGKZDgBs6aT9lu8QknRCImkEAPZmv5QSQpJLSCSNAMA5U6Udy3cISRYh8UIgAHC4kxCSgeTz+fn5+XQ6faCQrKoiYfnll+H+IwgEgatgVQd21K1ASEYRjUa9Xm8gEOjs7AyHwzUuv4VVJBgFCAJB4CpY3oFKSgkhGUKpVPJ4PNlsVhwXi0W3253L5fa8\/NYmjRgFCAJB4CrI0AFt+e7f3e0ISX\/i8biYGFVOBwYGpqamdrQRU6IvX3vp975\/+wUAgOP511OD\/9L7J4SkP3Nzc319fZXTYDA4MjKyuxnFXgEAFEJJIUUikf7+\/srp8DZcSwAAhGQ20Wi0t7e3eoYUCoW4lgAACMlsksmkz+ernAo5CUVxLQEAEJLZlMtlIaREIiGOM5mMy+VaXV3lWgIAICRrJkler7enp6etrS0Wi3EhAQAQEgAAAEICAACEBAAAgJAa4N69e6b9Wwe+7zWXy4kGS0tLFvahWCwuLCwkk0mrOqCRSqVWVlas6kM2mxUNHj16ZPktYcuf3fJ7wPKnQJLRwNrxECHtZGJionqDuKEc+L7XyclJrUFXV5dBv8l7YB8SiUR7e\/ulS5fOnj17\/vz5crlscgcqY2Jra6t4Gi0JwtWrV8W3RIPTp0+\/8cYbGxsbVt0StvzZLb8HLH8KJBkNrB0PEdKPWFtbu3z5ssfjMecCHPi+V3HTt7S0ZDIZcby+vi6Odf9kdGAfRAPxHD548EA77e7u1neDYp0vvd3c3BSjYUdHhxGD0YF9EGEX46C4PbTTU6dORSIRS24JW\/7slt8Dlj8FkowG1o6HCGknoVDo+vXr4lYz5wIc+L5XcQsePXpUzOK1p1GMC6lUyuQ+LCwsiI+EFgZBY3R0dGxs7MKFC0YI6cA+FAqFxcXF6gaiMxZGw2Y\/u+X3gOVPgSSjgbXjIULaiTYNF3Nzcy5APe97nZmZEZ9Jxfz93LlzV69eNb8Pd+7cGRwcFF90uVzis9KNGzfMD4L4ZHrmzBlxYJCQ6nzxrsazZ8\/EWGDQh9OGemKbn93ye8Dyp0CS0cDa8RAh7Y1pF6Ce970GAgFx8926dUs8h36\/\/8WLFyb3QXwsbWlpEU+COE6n08eOHbt\/\/76ZHVhfX+\/q6tKWLwwSUv0v3l1eXu7o6Pjggw8svCXs97Nbfg9Y\/hRIMhpYOx4ipG+ztZ5tqiNu6AWo\/hcPfN\/rwsLCyZMnS6WSdipuwfHxcZP7cPv27VdeeaW6gcDMDojHUnw4TWxz9uxZ8fFQl+1nDfVB48mTJ+3t7ZOTk8bdkOa\/AliGn92qe6D+DhjxFDTaB4NGA4QkEU+fPl3cppKuNPoCVP+LB77v9c6dO9WzeHGDDg0NmdyHu3fvVj+Kunxmb6gDYvS58D1iTDxz5owuw2JDfRCIluJzsYiGoTek+a8AluFnt+oeqL8DRjwFjfbBoNEAIbFk9x37ve81lUoVCoWt7Q1O4oti6NRWLbq7u8VNaXIfNjc3jx8\/Ho\/Ht7Y3\/5w4cULf38M4sAPVGLRkd2Af8vm8mEuJIGx+T+WDqjk9sfAmNOFnt\/wesPwpkGQ0QEiOFtLWPu97FZPxys7amZkZ8S2twTvvvGNJHz777LOOjo7XX39dNJiYmDC\/A0YL6cA+jI6O7ijtfOXKFTNvCQtvQnN+dsvvAcufAklGA4QEWy9evKjxe3biW6KBQR\/J6+yD4J\/\/\/KehfTiwA5ZfCHv3hPjL8BRIMhqoBUICAACEBAAAgJAAAAAhAQAAICQAAEBIAAAACAkAABASAAAAQgKwnL\/+9a8EAQAhAVhPKBRKJBIrKyuEAgAhAVjDo0eP3nzzzQ8\/\/PDKlSsLCwsEBAAhAVhAoVDweDwGVaQFQEgAUC8ff\/xxdekdAEBIANbw+PHjY8eOra6uPnjwgJc6AyAkAMtYW1u7dOnS0NDQ8PDwn\/\/8ZwICgJAALGBzczMYDG5sbBAKAIQEYCVTU1NaEWsAQEgAVhIIBD755BPiAICQACxmcXHR5XKNjo4+ffqUaAAgJAArSafTwWCwtbX17Nmzjx8\/JiAACAnASjY2Nt5\/\/\/3f\/e53hAIAIQFYxpMnT27fvj02NpZMJokGAEICsJhyuUwQABASAAAgJAAAAIQEAABq8V+QqKoTDszZ0wAAAABJRU5ErkJggg==","width":339}
%---
