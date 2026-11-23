classdef Rod3D < matlab.mixin.SetGetExactNames
    % ROD3D : define the properties of cross section of beam/rod
    
    properties
        k   % Shear correction factor, [2 x 1] vector
        A   % Area of the cross-section
        I1  % Polar moment of inertia         
        I2  % Moment of inertia at 2-direction
        I3  % Moment of inertia at 3-direction 
    end
    
    methods
        function obj = Rod3D(A, I1, I2, I3, k)
            
            if (A < 0) || (I2 < 0) || (I3 < 0)
                error('CrossSection:InvalidGeometry', ...
                      'The geometric parameters of the cross-section are incorrect: A=%.2e, I1=%.2e, I2=%.2e', ...
                      A, I2, I3);
            end

            obj.k  = k;
            obj.A  = A;
            obj.I1 = I1;             
            obj.I2 = I2;
            obj.I3 = I3;

        end
    end
end