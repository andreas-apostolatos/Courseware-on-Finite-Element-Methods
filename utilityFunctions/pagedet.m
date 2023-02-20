function pdet = pagedet(M)
%% PAGEDET Page-wise computation of the determinant for 2x2 or 3x3 page-wise matrices
%
% Returns the determinants of 2x2 or 3x3 matrices page-wise (by rule of 
% Sarrus).
%
%  Input :
%      M : The 3d array of matrices
% 
% Output :
%  pwdet : The vector of determinants. The size( det ) is [n, 1, 1].
% 
%% Function Implementation
if size(M, 1) == 3 && size(M, 2) == 3 && (length(size(M)) == 3 || length(size(M)) == 2) % length(size(M)) == 2 in case of one single element
	is3d = true;
elseif size(M, 1) == 2 && size(M, 2) == 2 && (length(size(M)) == 3 || length(size(M)) == 2) % length(size(M)) == 2 in case of one single element
	is3d = false;
else
	error("The array dimensions match neither a 2x2 nor a 3x3 page-wise matrix")
end
if is3d
	pdet = M(1, 1, :) .* M(2, 2, :) .* M(3, 3, :) ...
		+ M(1, 2, :) .* M(2, 3, :) .* M(3, 1, :) ...
		+ M(1, 3, :) .* M(2, 1, :) .* M(3, 2, :) ...
		- M(3, 1, :) .* M(2, 2, :) .* M(1, 3, :) ...
		- M(3, 2, :) .* M(2, 3, :) .* M(1, 1, :) ...
		- M(3, 3, :) .* M(2, 1, :) .* M(1, 2, :);
else
	pdet = M(1, 1, :) .* M(2, 2, :) - M(1, 2, :) .* M(2, 1, :);
end

end