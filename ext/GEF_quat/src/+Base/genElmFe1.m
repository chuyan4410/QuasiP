function Fe1 = genElmFe1(c, e, lambd, G, S_prim, r_prim, N, e_prim, ...
                            N_prim, CI_mB, kap0)

[~, dkap_dq] = Base.genElmdgam_dq(e, G, r_prim, N, e_prim, ...
    lambd, S_prim, N_prim);

% ---------------- vectorization ----------------
Fe1 = -c^2 * dkap_dq.' * CI_mB * kap0;

end