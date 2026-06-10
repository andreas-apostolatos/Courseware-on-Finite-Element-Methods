%[text] %[text:anchor:T_5395A80D] # Assembly of Sparse-Matrices for Finite Element Formulations
%[text:tableOfContents]{"heading":"**Table of Contents**"}
%[text] It receives as input an array of page-wise matrices corresponding to the element stiffness or element consistent nodal force vectors and returns their assembled master stiffness or master consistent nodal force vectors using sparse-matrices
function varargout = assembleSparseMatrices ...
    (EFT, numDOFs, numDOFsEl, numEl, varargin)
%[text] %[text:anchor:H_2D911AB8] ## Function Description
%[text] Returns the assembled matrices and vectors as sparse-matrices using the corresponding element constituents provided in the form of page-wise arrays.
%[text]  **Input** **:**
%[text]  `EFT` : Array `[numDOFsEl x numEl]` corresponding to the Element Freedom Table for each element in the mesh
%[text]  `numDOFs` : Number of Degrees of Freedom (DOFs) in the mesh
%[text]  `numDOFsEl` : Number of DOFs per element
%[text]  `numEl` : Number of elements in the mesh
%[text]  `varargin` : Array of input page-wise matrices and/or vectors
%[text]  **Output :**
%[text]  `varargout` : Array of output sparse matrices and/or vectors
%[text] %[text:anchor:H_006BB52A] ## Function Implementation
%[text] %[text:anchor:H_C12AEF4C] ## Input validation
arguments
    EFT (:, :) double
    numDOFs double {mustBeInteger, mustBePositive}
    numDOFsEl double {mustBeInteger, mustBePositive}
    numEl double {mustBeInteger, mustBePositive}
end
arguments (Repeating)
    varargin
end
%[text] %[text:anchor:H_F9A09E8E] ### Check input
    szEFT = size(EFT);
    if numel(szEFT) ~= 2
        error("The Element Freedom Tables should be organized in a matrix format " + ...
            "[numDOFsEl x numEl] but the provided array EFT has %i dimensions", szEFT);
    else
        if szEFT(1) ~= numDOFsEl
            error("The Element Freedom Tables' array has %i rows which does not correspond " + ...
                "to the number of DOFs per element which is %i", szEFT(1), numDOFsEl);
        end
        if szEFT(2) ~= numEl
            error("The Element Freedom Tables' array has %i columns which does not " + ...
                "correspond to the number of elements in the mesh which is %i", szEFT(1), numEl);
        end
    end
%[text] %[text:anchor:H_4D567000] ### Initialize output
    varargout = cell(nargin - 4, 1);
%[text] %[text:anchor:H_729A67C3] ### Loop over all input page-wise arrays
    for ii = 1:nargin - 4
%[text] %[text:anchor:H_C9306546] #### Verification that the input is a page-wise matrix/vector
        sz = size(varargin{ii});
        if numel(sz) ~= 3 && numEl ~= 1
            error("The %i input is not a page-wise matrix", ii + 4)
        end
%[text] %[text:anchor:H_A24A052D] #### Verification that the input page-wise matrix/vector has the expected number of rows
        if sz(1) ~= numDOFsEl
            error("The input page-wise matrix/vector has %i number of rows which does " + ...
                "not correspond to the number of DOFs per element which is %i", sz(1), numDOFsEl)
        end
%[text] %[text:anchor:H_0141A805] #### Verification that the input page-wise matrix/vector has the expected number of rows
        if sz(2) ~= numDOFsEl && sz(2) ~= 1
            error("The page-wise matrix has %i number of columns neither corresponds to " + ...
                "the number of DOFs per element which is %i (page-wise matrix) nor to 1 (page-wise vector)", sz(2), numDOFsEl)
        end
%[text] %[text:anchor:H_E09DBFD3] #### Verification that the input page-wise matrix/vector has the expected number of pages
        if numEl ~= 1
            if sz(3) ~= numEl
                error("The input page-wise matrix/vector has %i number of pages but the " + ...
                    "number of elements in the mesh is %i", sz(3), numEl)
            end
        end
%[text] %[text:anchor:H_87E15ABA] #### Assembly of the global constituent depending on whether the input is a page-wise matrix or page-wise array
        if sz(2) == 1
%[text] %[text:anchor:H_9404E619] #### Assembly of the global consistent nodal force vector using sparse matrices
            iIdxLoadVct = EFT(:);
            sLoadVct = permute(reshape(reshape(varargin{ii}, ...
                [numDOFsEl,1,numEl]), [1, 1, numDOFsEl*numEl]),[3, 2, 1]); % Non-zero entries of the element load vectors
            varargout{ii} = sparse(iIdxLoadVct, ones(numEl*numDOFsEl, 1), sLoadVct); % Assembly to the global consistent nodal force vector
%[text] 
        else
%[text] %[text:anchor:H_FFCC8862] #### Assembly of the global stiffness matrix using sparse matrices
            iIdx = permute(reshape(reshape(permute(pagetranspose(reshape(repmat(EFT, ...
                numDOFsEl, 1, 1), [numDOFsEl numDOFsEl numEl])), [2, 1, 3]), [numDOFsEl, 1, ...
                numDOFsEl*numEl]), [1, 1, numDOFsEl^2*numEl]), [3, 2, 1]); % Indices of the non-zeros entries row-wise
            jIdx = permute(reshape(reshape(pagetranspose(reshape(repmat(EFT, numDOFsEl, ...
                1, 1),[numDOFsEl numDOFsEl numEl])), [numDOFsEl, 1, numDOFsEl*numEl]), ...
                [1, 1, numDOFsEl^2*numEl]),[3, 2, 1]); % Indices of the non-zeros entries column-wise
            sStiffMtx = permute(reshape(reshape(varargin{ii}, ...
                [numDOFsEl, 1, numDOFsEl*numEl]),[1, 1, numDOFsEl^2*numEl]),[3, 2, 1]); % Non-zero entries of the element stiffness matrices
            varargout{ii} = sparse(iIdx, jIdx, sStiffMtx, numDOFs, numDOFs); % Assembly to the global stiffness matrix
        end
%[text] 
    end
%[text] 
end

%[appendix]{"version":"1.0"}
%---
