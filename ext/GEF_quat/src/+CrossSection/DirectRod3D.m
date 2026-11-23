classdef DirectRod3D < matlab.mixin.SetGetExactNames
    %DIRECTROD3D Section properties supplied as beam resultants.
    % Use this class for composite or experimentally identified sections
    % where EA, GA, GJ, EI, and mass resultants are known directly.

    properties
        rhoA    % Linear density
        Jrho    % Section rotary inertia matrix per unit length
        CA_mB   % Stretch-shear stiffness matrix
        CI_mB   % Torsion-bending stiffness matrix
    end

    methods
        function obj = DirectRod3D(rhoA, Jrho, CA_mB, CI_mB)
            if rhoA <= 0
                error('DirectRod3D:InvalidMass', ...
                    'rhoA must be positive.');
            end

            obj.rhoA = rhoA;
            obj.Jrho = obj.asMatrix3(Jrho, 'Jrho');
            obj.CA_mB = obj.asMatrix3(CA_mB, 'CA_mB');
            obj.CI_mB = obj.asMatrix3(CI_mB, 'CI_mB');
        end
    end

    methods (Static, Access = private)
        function value = asMatrix3(value, name)
            if isvector(value) && numel(value) == 3
                value = diag(value(:));
            end

            if ~isequal(size(value), [3, 3])
                error('DirectRod3D:InvalidMatrix', ...
                    '%s must be a 3-by-3 matrix or a 3-vector.', name);
            end

            if any(~isfinite(value(:)))
                error('DirectRod3D:InvalidMatrix', ...
                    '%s contains non-finite values.', name);
            end
        end
    end
end
