function [dgam_dq, dkap_dq] = genElmdgam_dq(e, G, r_prim, N, e_prim, ...
    lambd, S_prim, N_prim)

dA_de = Base.genElmdA_de(e);
dG_de = Base.genElmdG_de();
term1 = zeros(3,21);
term2 = zeros(3,21);
for alpha = 1:4
    term1 = term1 + dA_de(:,:,alpha).' * r_prim * N(alpha,:);
    term2 = term2 + 2 * dG_de(:,:,alpha) * e_prim * N(alpha,:);
end
dgam_dq = lambd.' * S_prim + term1;
dkap_dq = 2 * G * N_prim + term2;
end

