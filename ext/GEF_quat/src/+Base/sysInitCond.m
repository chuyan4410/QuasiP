function [q0_sys, dq0_sys, ddq0_sys] = sysInitCond(nodes)

q_cells   = cellfun(@(nd) nd.q(:,1), nodes, 'UniformOutput', false);
q0_sys    = vertcat(q_cells{:});
dq_cells  = cellfun(@(nd) nd.dq(:,1), nodes, 'UniformOutput', false);
dq0_sys   = vertcat(dq_cells{:});
ddq_cells = cellfun(@(nd) nd.ddq(:,1), nodes, 'UniformOutput', false);
ddq0_sys  = vertcat(ddq_cells{:});

end