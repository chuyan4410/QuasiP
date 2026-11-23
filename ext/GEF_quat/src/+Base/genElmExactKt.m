function Kt = genElmExactKt(c, e, G, r_prim, N, e_prim, ...
                            lambd, S_prim, N_prim, CA_mB, ...
                            CI_mB, gam, gam0, kap, kap0, numDOF)

    % 1. Material Stiffness (K_hat)
    [dgam_dq, dkap_dq] = Base.genElmdgam_dq(e, G, r_prim, N, e_prim, ...
                                            lambd, S_prim, N_prim);
    Kt1 = dgam_dq.' * CA_mB * dgam_dq + dkap_dq.' * CI_mB * dkap_dq;

    % 2. Geometric Stiffness (K_check)
    % Calculate the internal forces and moments based on current strains
    Fc = CA_mB * (gam - gam0); % Adjust to (gam - gam0) if your constitutive law requires
    Mc = CI_mB * (kap - kap0); % Adjust to (kap - kap0) if your constitutive law requires
    
    Kt2 = Base.genElmKt_check(e, r_prim, N, S_prim, N_prim, Fc, Mc, numDOF);

    % 3. Total Tangent Stiffness Matrix
    Kt = c^2 * (Kt1 + Kt2);
end