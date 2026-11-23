function Damp = genDamp(elements, nodes, constraints, dampRatio)
[q0_sys, dq0_sys, ddq0_sys] = Base.sysInitCond(nodes);
dof_sys = length(q0_sys);
myAssembly = Assembly.Assembly(elements, nodes, constraints, ...
             q0_sys, dq0_sys, ddq0_sys);
myAssembly.genSysMatVec(2);

M  = myAssembly.M;
Kt = myAssembly.Kt;
Cq = myAssembly.Cq;
[U1,la] = eigs(Cq.' * Cq, dof_sys);
idx = diag(la) == 0;
U = U1(:,idx);
M_hat  = U.' * M  * U;
Kt_hat = U.' * Kt * U;
E_max = 10;
if dof_sys < E_max
    E_max = dof_sys;
end
[~,D2_hat] = eigs(Kt_hat, M_hat, E_max, 'smallestabs');
D2 = sqrt(diag(D2_hat));

alpha = (2 * D2(1) * D2(2) * (D2(2) * dampRatio - D2(1) * dampRatio)) / (D2(2)^2-D2(1)^2);
beta = 2 * (D2(2) * dampRatio - D2(1) * dampRatio) / (D2(2)^2 - D2(1)^2);
Damp = alpha * M + beta * Kt;
end