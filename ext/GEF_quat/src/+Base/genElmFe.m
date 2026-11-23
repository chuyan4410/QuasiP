function Fe = genElmFe(c, e, lambd, G, S_prim, r_prim, N, e_prim, ...
                            N_prim, CA_mB, gam, gam0, CI_mB, kap, kap0)

[dgam_dq, dkap_dq] = Base.genElmdgam_dq(e, G, r_prim, N, e_prim, ...
    lambd, S_prim, N_prim);

% ---------------- vectorization ----------------
Fe = c^2 * (dgam_dq.' * CA_mB * (gam - gam0) + dkap_dq.' * CI_mB * (kap - kap0));

end