function Ct = genElmCt(dq_elm, N, Jrho, G, numDof)

dM_dq = Base.genElmdM_dq(N, Jrho, G, numDof);     
dM    = Base.genElmdM(dM_dq, dq_elm, numDof); 

tmp = reshape(dM_dq, numDof, []);     
T   = reshape(dq_elm.' * tmp, numDof, numDof);

Ct = dM + T - T.';

end
