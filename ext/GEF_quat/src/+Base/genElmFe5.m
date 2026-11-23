function Fe5 = genElmFe5(c, e, lambd, G, S_prim, r_prim, N, e_prim, ...
                            N_prim, CA_mB, gam)

[dgam_dq, ~] = Base.genElmdgam_dq(e, G, r_prim, N, e_prim, ...
    lambd, S_prim, N_prim);

% ---------------- vectorization ----------------
Fe5 = c^2 * dgam_dq.' * CA_mB * gam;

end