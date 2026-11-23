classdef IsotropicLE < matlab.mixin.SetGetExactNames
    % ISOTROPICLE 
    % define isotropic linear elastic material
    
    properties
        E       % Young's modulus 
        nu      % Poisson's ratio 
        G       % Shear modulus 
        Lame    % Lame constant 
        density % Density 
    end
    
    methods
        function obj = IsotropicLE(E, nu, density)
            
            if (E < 0.0) || (nu > 0.5) || (nu < -1.0) || (density < 0.0)
                error('Material:InvalidProperties', ...
                      'Isotropic linear elastic constant is incorrect: E=%.2f, nu=%.2f, Rho=%.2f', ...
                      E, nu, density);
            end
            
            obj.E = E;
            obj.nu = nu;
            obj.density = density;
            
            obj.G = E / (2 * (1 + nu));

            obj.Lame = (nu * E) / ((1 + nu) * (1 - 2*nu));
        end
    end
end

