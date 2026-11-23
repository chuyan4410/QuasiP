function Fe3 = genElmFe3(c, e, lambd, G, S_prim, r_prim, N, e_prim, ...
                            N_prim, CI_mB, kap)

[~, dkap_dq] = Base.genElmdgam_dq(e, G, r_prim, N, e_prim, ...
    lambd, S_prim, N_prim);

% ---------------- vectorization ----------------
Fe3 = c^2 * dkap_dq.' * CI_mB * kap;

end