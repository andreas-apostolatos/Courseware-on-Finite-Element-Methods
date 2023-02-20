function mustHaveTimoshenkoBeamProperties(propStr)
%% MUSTHAVETIMOSHENKOBEAMPROPERTIES Checks whether the input variable is a struct having all properties needed for the Timoshenko beam
%
% An input validation function that verifies whether the input variable is
% a MATLAB®-struct containing all fields needed to sufficiently define the
% properties of a Timoshenko beam
%
%% Function Implementation
if ~isstruct(propStr)
    error("Input variable 'propStr' should be a MATLAB-struct " + ...
        "containing all properties needed to fully define a " + ...
        "Timoshenko beam")
else
    % Thickness
    if ~isfield(propStr, 'A')
        error("struct-variable 'propStr' should have a field 'A' representing the " + ...
            "beam's cross sectional area")
    else
        if ~isnumeric(propStr.A)
            error("The beam's cross sectional area must be a numeric variable")
        else
            if ~isscalar(propStr.A)
                error("The beam's cross sectional area must be a scalar variable")
            end
        end
    end

    % Vertical distributed load
    if ~isfield(propStr, 'qBar')
        error("struct-variable 'propStr' should have a field 'qBar' representing the " + ...
            "vertical distributed across the beam's length")
    else
        if ~isnumeric(propStr.qBar)
            error("The vertical distributed load must be a numeric variable")
        else
            if ~isscalar(propStr.qBar)
                error("The vertical distributed load must be a scalar variable")
            end
        end
    end

    % Moment
    if ~isfield(propStr, 'mBar')
        error("struct-variable 'propStr' should have a field 'mBar' representing the " + ...
            "distributed moment across the beam's length")
    else
        if ~isnumeric(propStr.mBar)
            error("The distributed moment must be a numeric variable")
        else
            if ~isscalar(propStr.mBar)
                error("The distributed moment must be a scalar variable")
            end
        end
    end

    % Beam's Young's modulus
    if ~isfield(propStr, 'E')
        error("struct-variable 'propStr' should have a field 'E' representing the " + ...
            "beam's Young's modulus")
    else
        if ~isnumeric(propStr.E)
            error("The beam's Young's modulus must be a numeric variable")
        else
            if ~isscalar(propStr.E)
                error("The beam's Young's modulus must be a scalar variable")
            else
                if propStr.E <= 0
                    error("The beam's Young's modulus must be strictly positive")
                end
            end
        end
    end

    % Beam's Poisson ratio
    if ~isfield(propStr, 'nu')
        error("struct-variable 'propStr' should have a field 'nu' representing the " + ...
            "beam's Poisson ratio")
    else
        if ~isnumeric(propStr.nu)
            error("The beam's Poisson ratio must be a numeric variable")
        else
            if ~isscalar(propStr.nu)
                error("The beam's Poisson ratio must be a scalar variable")
            else
                if propStr.nu < 0 || propStr.nu > 1
                    error("The beam's Poisson ratio must have a value in the closed interval [0, 1]")
                end
            end
        end
    end

    % Beam's moment of inertia
    if ~isfield(propStr, 'I')
        error("struct-variable 'propStr' should have a field 'I' representing the " + ...
            "beam's moment of inertia")
    else
        if ~isnumeric(propStr.I)
            error("The beam's moment of inertia must be a numeric variable")
        else
            if ~isscalar(propStr.I)
                error("The beam's moment of inertia must be a scalar variable")
            else
                if propStr.I <= 0
                    error("The beam's moment of inertia must be strictly positive")
                end
            end
        end
    end

    % Beam's shear modulus
    if ~isfield(propStr, 'G')
        error("struct-variable 'propStr' should have a field 'G' representing the " + ...
            "beam's shear modulus")
    else
        if ~isnumeric(propStr.G)
            error("The beam's shear modulus must be a numeric variable")
        else
            if ~isscalar(propStr.G)
                error("The beam's shear modulus must be a scalar variable")
            else
                if propStr.G <= 0
                    error("The beam's shear modulus must be strictly positive")
                end
            end
        end
    end

    % Shear correction factor
    if ~isfield(propStr, 'alpha')
        error("struct-variable 'propStr' should have a field 'alpha' representing the " + ...
            "beam's shear correction factor")
    else
        if ~isnumeric(propStr.alpha)
            error("The beam's shear correction factor must be a numeric variable")
        else
            if ~isscalar(propStr.alpha)
                error("The beam's shear correction factor must be a scalar variable")
            else
                if propStr.alpha < 2/3 || propStr.alpha > 1
                    error("The beam's shear correction factor must lie in the closed interval [2/3 1]")
                end
            end
        end
    end

end

end