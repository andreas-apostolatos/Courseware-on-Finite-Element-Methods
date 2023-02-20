function mustHaveReissnerMindlinPlateProperties(propStr)
%% MUSTHAVEREISSNERMINDLINPLATEPROPERTIES Checks whether the input variable is a struct having all properties needed for the Reissner-Mindlin plate
%
% An input validation function that verifies whether the input variable is
% a MATLAB®-struct containing all fields needed to sufficiently define the
% properties of a Reissner-Mindlin plate
%
%% Function Implementation
if ~isstruct(propStr)
    error("Input variable 'propStr' should be a MATLAB-struct " + ...
        "containing all properties needed to fully define a " + ...
        "Reissner-Mindlin plate")
else
    % Thickness
    if ~isfield(propStr, 't')
        error("struct-variable 'propStr' should have a field 't' representing the " + ...
            "plate's thickness")
    else
        if ~isnumeric(propStr.t)
            error("The plate's thickness must be a numeric variable")
        else
            if ~isscalar(propStr.t)
                error("The plate's thickness must be a scalar variable")
            end
        end
    end

    % Vertical distributed load
    if ~isfield(propStr, 'pBar')
        error("struct-variable 'propStr' should have a field 'pBar' representing the " + ...
            "vertical distributed over the plate's surface")
    else
        if ~isnumeric(propStr.pBar)
            error("The vertical distributed load must be a numeric variable")
        else
            if ~isscalar(propStr.pBar)
                error("The vertical distributed load must be a scalar variable")
            end
        end
    end

    % Moment around the X-axis
    if ~isfield(propStr, 'mxBar')
        error("struct-variable 'propStr' should have a field 'mxBar' representing the " + ...
            "distributed moment around the X-axis over the plate's surface")
    else
        if ~isnumeric(propStr.mxBar)
            error("The distributed moment around the X-axis  must be a numeric variable")
        else
            if ~isscalar(propStr.mxBar)
                error("The distributed moment around the X-axis must be a scalar variable")
            end
        end
    end

    % Moment around the Y-axis
    if ~isfield(propStr, 'myBar')
        error("struct-variable 'propStr' should have a field 'myBar' representing the " + ...
            "distributed moment around the Y-axis over the plate's surface")
    else
        if ~isnumeric(propStr.myBar)
            error("The distributed moment around the Y-axis load must be a numeric variable")
        else
            if ~isscalar(propStr.myBar)
                error("The distributed moment around the Y-axis load must be a scalar variable")
            end
        end
    end

    % Plate's Young's modulus
    if ~isfield(propStr, 'E')
        error("struct-variable 'propStr' should have a field 'E' representing the " + ...
            "plate's Young's modulus")
    else
        if ~isnumeric(propStr.E)
            error("The plate's Young's modulus must be a numeric variable")
        else
            if ~isscalar(propStr.E)
                error("The plate's Young's modulus must be a scalar variable")
            else
                if propStr.E <= 0
                    error("The plate's Young's modulus must be strictly positive")
                end
            end
        end
    end

    % Plate's Poisson ratio
    if ~isfield(propStr, 'nu')
        error("struct-variable 'propStr' should have a field 'nu' representing the " + ...
            "plate's Poisson ratio")
    else
        if ~isnumeric(propStr.nu)
            error("The plate's Poisson ratio must be a numeric variable")
        else
            if ~isscalar(propStr.nu)
                error("The plate's Poisson ratio must be a scalar variable")
            else
                if propStr.nu < 0 || propStr.nu > 1
                    error("The plate's Poisson ratio must have a value in the closed interval [0, 1]")
                end
            end
        end
    end

    % Plate's shear modulus
    if ~isfield(propStr, 'G')
        error("struct-variable 'propStr' should have a field 'G' representing the " + ...
            "plate's shear modulus")
    else
        if ~isnumeric(propStr.G)
            error("The plate's shear modulus must be a numeric variable")
        else
            if ~isscalar(propStr.G)
                error("The plate's shear modulus must be a scalar variable")
            else
                if propStr.G <= 0
                    error("The plate's shear modulus must be strictly positive")
                end
            end
        end
    end

    % Plate's stiffness
    if ~isfield(propStr, 'D')
        error("struct-variable 'propStr' should have a field 'D' representing the " + ...
            "plate's stiffness")
    else
        if ~isnumeric(propStr.D)
            error("The plate's stiffness must be a numeric variable")
        else
            if ~isscalar(propStr.D)
                error("The plate's stiffness must be a scalar variable")
            else
                if propStr.D <= 0
                    error("The plate's stiffness must be strictly positive")
                end
            end
        end
    end

    % Shear correction factor
    if ~isfield(propStr, 'alpha')
        error("struct-variable 'propStr' should have a field 'alpha' representing the " + ...
            "plate's shear correction factor")
    else
        if ~isnumeric(propStr.alpha)
            error("The plate's shear correction factor must be a numeric variable")
        else
            if ~isscalar(propStr.alpha)
                error("The plate's shear correction factor must be a scalar variable")
            else
                if propStr.alpha < 2/3 || propStr.alpha > 1
                    error("The plate's shear correction factor must lie in the closed interval [2/3 1]")
                end
            end
        end
    end

end

end