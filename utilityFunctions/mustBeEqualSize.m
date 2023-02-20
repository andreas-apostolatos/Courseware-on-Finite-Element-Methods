function mustBeEqualSize(a, b)
%% MUSTBEEQUALSIZE Checks whether the input variables have equal size
%
% An input validation function that verifies whether the input arrays have 
% equal size, see also here:
% https://www.mathworks.com/help/matlab/matlab_prog/function-argument-validation-1.html#mw_d3c46813-2f36-4bfe-8cfa-ba9157619c06
%
%% Function Implementation

if ~isequal(size(a),size(b))
    eid = 'Size:notEqual';
    msg = 'Size of first input must equal size of second input.';
    throwAsCaller(MException(eid, msg))
end

end