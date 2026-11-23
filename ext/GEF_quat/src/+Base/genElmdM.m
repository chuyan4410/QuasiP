function dM = genElmdM(dM_dq, dq_elm, numDof)
tmp = reshape(dM_dq, [], numDof);   
dM = reshape(tmp * dq_elm, numDof, numDof);
end