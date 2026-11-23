function d2A_de2 = genElmd2A_de2()
% GENELMD2A_DE2 Generates the constant 2nd-order derivatives of rotation matrix
% Output is a 3x3x4x4 tensor, where d2A_de2(:,:,alpha,beta) = d^2 A / (de_alpha de_beta)

    d2A_de2 = zeros(3,3,4,4);
    
    % e_alpha = e0 (index 1)
    d2A_de2(:,:,1,1) =  2 * [ 1  0  0;  0  1  0;  0  0  1];
    d2A_de2(:,:,1,2) =  2 * [ 0  0  0;  0  0 -1;  0  1  0];
    d2A_de2(:,:,1,3) =  2 * [ 0  0  1;  0  0  0; -1  0  0];
    d2A_de2(:,:,1,4) =  2 * [ 0 -1  0;  1  0  0;  0  0  0];
    
    % e_alpha = e1 (index 2)
    d2A_de2(:,:,2,1) =  2 * [ 0  0  0;  0  0 -1;  0  1  0];
    d2A_de2(:,:,2,2) =  2 * [ 1  0  0;  0 -1  0;  0  0 -1];
    d2A_de2(:,:,2,3) =  2 * [ 0  1  0;  1  0  0;  0  0  0];
    d2A_de2(:,:,2,4) =  2 * [ 0  0  1;  0  0  0;  1  0  0];

    % e_alpha = e2 (index 3)
    d2A_de2(:,:,3,1) =  2 * [ 0  0  1;  0  0  0; -1  0  0];
    d2A_de2(:,:,3,2) =  2 * [ 0  1  0;  1  0  0;  0  0  0];
    d2A_de2(:,:,3,3) =  2 * [-1  0  0;  0  1  0;  0  0 -1];
    d2A_de2(:,:,3,4) =  2 * [ 0  0  0;  0  0  1;  0  1  0];

    % e_alpha = e3 (index 4)
    d2A_de2(:,:,4,1) =  2 * [ 0 -1  0;  1  0  0;  0  0  0];
    d2A_de2(:,:,4,2) =  2 * [ 0  0  1;  0  0  0;  1  0  0];
    d2A_de2(:,:,4,3) =  2 * [ 0  0  0;  0  0  1;  0  1  0];
    d2A_de2(:,:,4,4) =  2 * [-1  0  0;  0 -1  0;  0  0  1];
end