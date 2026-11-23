function Kt_check = genElmKt_check(e, r_prim, N, S_prim, N_prim, Fc, Mc, numDOF)
% GENELMKT_CHECK Computes the geometric stiffness matrix (K_check) 
% using highly optimized vector outer products.
%
% Inputs:
% e      : 4x1 quaternion vector
% r_prim : 3x1 derivative of position vector w.r.t arc length
% N      : 4x21 shape function matrix for quaternions
% S_prim : 3x21 shape function derivative matrix for position
% N_prim : 4x21 shape function derivative matrix for quaternions
% Fc     : 3x1 internal force vector (CA_mB * gam)
% Mc     : 3x1 internal moment vector (CI_mB * kap)
%
% Output:
% Kt_check : 21x21 geometric stiffness matrix

    dA_de = Base.genElmdA_de(e);
    dG_de = Base.genElmdG_de();
    d2A_de2 = Base.genElmd2A_de2(); % Constant 3x3x4x4 tensor

    Kt_check = zeros(numDOF, numDOF);

    % ---------------------------------------------------------------------
    % 1. Geometric stiffness from Curvature (Moment term)
    % ---------------------------------------------------------------------
    for alpha = 1:4
        % V2_alpha is a 1x21 row vector
        V2_alpha = Mc.' * dG_de(:,:,alpha) * N_prim;
        N_alpha = N(alpha, :); % 1x21 row vector
        
        % Fast outer product for the symmetric cross term
        Kt_check = Kt_check + 2 * (V2_alpha.' * N_alpha + N_alpha.' * V2_alpha);
    end

    % ---------------------------------------------------------------------
    % 2. Geometric stiffness from Strain (Force term)
    % ---------------------------------------------------------------------
    for alpha = 1:4
        % V1_alpha is a 1x21 row vector
        V1_alpha = Fc.' * dA_de(:,:,alpha).' * S_prim;
        N_alpha = N(alpha, :);
        
        % Fast outer product for the symmetric cross term
        Kt_check = Kt_check + (V1_alpha.' * N_alpha + N_alpha.' * V1_alpha);

        % Term involving the 2nd derivative of the rotation matrix A
        for beta = 1:4
            % S_ab is a scalar coefficient
            S_ab = Fc.' * d2A_de2(:,:,alpha,beta).' * r_prim;
            N_beta = N(beta, :);
            
            % Rank-1 update using outer product of shape functions
            Kt_check = Kt_check + S_ab * (N_alpha.' * N_beta);
        end
    end
end