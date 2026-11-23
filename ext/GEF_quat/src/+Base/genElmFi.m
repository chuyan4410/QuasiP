function F_i = genElmFi(dq_elm, N, Jrho, G, numDof)

dM_dq = Base.genElmdM_dq(N, Jrho, G, numDof);  % 21 × 21 × 21
dM    = Base.genElmdM(dM_dq, dq_elm, numDof);  % 21 × 21
Fi2   = dM * dq_elm;                   % 21 × 1

dq = full(dq_elm(:));                        % keep column vector

% ---- dq^T * dM_dq(:,:,i) * dq，for all i ----
tmp = pagemtimes(dM_dq, dq);           % 21 × 1 × 21
quad = squeeze(pagemtimes(dq.', tmp)); % 1  × 21

F_i = 0.5 * quad - Fi2;              % 21 × 1
end