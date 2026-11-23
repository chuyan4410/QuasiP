function Fe2 = genElmFe2(c, e, lambd, G, S_prim, r_prim, N, e_prim, ...
                            N_prim, CA_mB, gam0)

[dgam_dq, ~] = Base.genElmdgam_dq(e, G, r_prim, N, e_prim, ...
    lambd, S_prim, N_prim);

% ---------------- vectorization ----------------
Fe2 = -c^2 * dgam_dq.' * CA_mB * gam0;

end