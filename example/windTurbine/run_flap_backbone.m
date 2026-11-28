function result = run_flap_backbone()
% Autonomous EPMC backbone of the 12.1 rpm blade, with 0.48% structural damping.
% Run from MATLAB with run_flap_backbone. CSV, FIG and PNG are saved beside this file.
    folder = fileparts(mfilename('fullpath'));
    previous = pwd;
    restoreFolder = onCleanup(@() cd(previous));
    root = fileparts(fileparts(folder));
    if isempty(which('coco_prob'))
        startupFile = fullfile(fileparts(root), 'coco_original', 'coco', 'startup.m');
        assert(isfile(startupFile), 'Run the COCO startup before calling run_flap_backbone.');
        run(startupFile);
    end
    cd(root);
    install;
    cd(folder);

    order = 5;
    etaRange = log10([0.05, 10]);
    m = build_blade('flap', 1, 12.1);
    m.F(:) = 0;
    n = m.n;
    total = 3 * n + m.nc;
    v = n + (1:n);
    b = 2 * n + (1:n);
    omega = m.omega;
    beta = 2 * 0.0048 * omega;

    % A bordered linear solve fixes unit flap-tip displacement for the seed.
    tipRow = sparse(1, m.tip, 1, 1, n);
    constraints = [m.H; tipRow];
    dynamic = m.K - omega^2 * m.M + 1i * omega * (m.G + m.D - beta * m.M);
    bordered = [dynamic, constraints.'; constraints, sparse(m.nc + 1, m.nc + 1)];
    rowScale = max(abs(bordered), [], 2);
    balanced = spdiags(1 ./ rowScale, 0, numel(rowScale), numel(rowScale)) * bordered;
    colScale = max(abs(balanced), [], 1).';
    rhs = [zeros(n + m.nc, 1); 1] ./ rowScale;
    mode = (balanced * spdiags(1 ./ colScale, 0, numel(colScale), numel(colScale)) ...
        \ rhs) ./ colScale;
    q = mode(1:n);
    velocity = 1i * omega * q;
    effectiveAcceleration = -omega^2 * q - beta * velocity;
    referenceSquared = real(q' * m.M * q);

    % b = dv/dt - beta*v introduces the exact feedback beta*M(q)*v.
    % The structural kernel receives b; physical viscous damping still acts on v.
    m.scale = 5 * m.scale;
    m.rows = [1 ./ m.scale(v); 1 ./ m.scale(b); 1e-5 * ones(n, 1); ones(m.nc, 1)];
    B = spdiags([ones(2 * n, 1); zeros(n + m.nc, 1)], 0, total, total);
    C = sparse(total, total);
    C(1:n, v) = -speye(n);
    C(v, b) = -speye(n);
    C(b, v) = m.D;
    left = spdiags(m.rows, 0, total, total);
    right = spdiags(m.scale, 0, total, total);
    B = left * B * right;
    C = left * C * right;
    m.feedbackScale = m.rows(v) .* m.scale(v);
    seed.k_wave = 1;
    seed.Xc = (10^etaRange(1) * [q; velocity; effectiveAcceleration; ...
        mode(n + (1:m.nc))] ./ m.scale).';

    phi = struct('param', 'MP', 'method', 'mVCF', 'isvar', true, ...
        'sparse_jacobian', true, 'type', {{'HB'}}, 'k', {{order}}, ...
        'P', {{1}}, 'S', {{32}}, 'S_wave', {{1024}});
    ode = struct('n', total, 'd', 1, 'e', 0, 'isfwithu', false);
    amp = struct('isamp', true, 'outdof', m.tip, 'Nsampling', 1024);
    prob = coco_prob();
    prob = coco_set(prob, 'QuasiP', 'ode_opts', ode, 'amp_opts', amp, ...
        'sta_opts', struct('issta', false), 'PHI_INITIAL', phi);
    prob = coco_set(prob, 'cont', 'PtMX', [0, 300], 'h0', 0.04, 'h_max', 0.20, ...
        'h_min', 1e-5, 'NSV', 10, 'TOL', 1e-8);
    prob = coco_set(prob, 'corr', 'ItMX', 20, 'TOL', 1e-8);
    prob = ode_isol2QuasiP(prob, '', @(x, p, t) epmc_rhs(x, p, t, m), ...
        @(x, p, t) state_jacobian(x, p, t, m), @(x, p, t) parameter_jacobian(x, m), ...
        @(p) operators(B, C), [], seed, omega, {'beta', 'eta'}, [beta; etaRange(1)]);

    [data, idx] = coco_get_func_data(prob, 'QuasiP', 'data', 'uidx');
    harmonics = 2 * order + 1;
    coordScale = spdiags(m.scale(1:n), 0, n, n);
    metric.mass = coordScale * m.M * coordScale / referenceSquared;
    metric.harmonics = harmonics;
    selected = idx([1:(n * harmonics), data.p_idx(2)]);
    prob = coco_add_func(prob, 'modal_amplitude', @amplitude_condition, metric, ...
        'zero', 'uidx', selected, 'F+dF');
    name = 'backbone_flap_12p1rpm_H5';
    bd = coco(prob, name, [], 1, {'eta', 'beta', 'u_om1'}, etaRange);

    eta = coco_bd_col(bd, 'eta');
    frequency = coco_bd_col(bd, 'u_om1');
    feedback = coco_bd_col(bd, 'beta');
    tip = blade_branch_amplitude(bd, m);
    assert(max(eta) >= etaRange(2) - 1e-6, 'Continuation stopped below target amplitude.');
    coeffs = coco_bd_col(bd, 'C', false);
    workError = zeros(numel(eta), 1);
    normError = workError;
    theta = (0:127).' * (2 * pi / 128);
    basis = ones(numel(theta), harmonics);
    basis(:, 2:2:end) = cos(theta * (1:order));
    basis(:, 3:2:end) = sin(theta * (1:order));
    for k = 1:numel(eta)
        displacement = coeffs{k}(:, 1:n) .* m.scale(1:n).';
        normError(k) = abs(sum((displacement * m.M) .* displacement, 'all') ...
            / (referenceSquared * 10^(2 * eta(k))) - 1);
        position = m.q + (basis * displacement).';
        speed = m.scale(v) .* (basis * coeffs{k}(:, v)).';
        % The force is affine in acceleration. Its imaginary part yields M(q)*v.
        momentum = imag(blade_force(m, position, speed, 1i * speed));
        dissipation = sum(speed .* (m.D * speed), 'all');
        workError(k) = abs(feedback(k) * sum(speed .* momentum, 'all') / dissipation - 1);
    end
    result = table(eta(:), sqrt(referenceSquared) * 10.^eta(:), frequency(:), ...
        tip(:), feedback(:) ./ (2 * frequency(:)), normError, workError, ...
        'VariableNames', {'log10_a_over_ref', 'modal_amplitude', 'omega_rad_s', ...
        'tip_flap_m', 'xi', 'normalization_error', 'relative_work_error'});
    writetable(result, [name '.csv']);
    fig = figure('Visible', 'off');
    plot(frequency, tip, 'LineWidth', 1.5);
    xlabel('Modal angular frequency [rad/s]');
    ylabel('Tip flap amplitude [m]');
    grid on;
    exportgraphics(fig, [name '.png'], 'Resolution', 180);
    set(fig, 'Visible', 'on');
    savefig(fig, [name '.fig']);
    close(fig);
end

function [f, J] = epmc_rhs(x, p, theta, m)
    v = m.n + (1:m.n);
    if nargout > 1
        [f, J] = blade_rhs(x, theta, m);
        feedback = diag(m.feedbackScale * p(1, 1));
        J(v, v, :) = J(v, v, :) + feedback;
    else
        f = blade_rhs(x, theta, m);
    end
    f(v, :) = m.feedbackScale .* x(v, :) .* p(1, :);
end

function J = state_jacobian(x, p, theta, m)
    [~, J] = epmc_rhs(x, p, theta, m);
end

function J = parameter_jacobian(x, m)
    v = m.n + (1:m.n);
    J = zeros(size(x, 1), 2, size(x, 2));
    J(v, 1, :) = reshape(m.feedbackScale .* x(v, :), m.n, 1, []);
end

function [flags, depends, A, B, C, Ap, Bp, Cp] = operators(B, C)
    flags = [0, 1, 1];
    depends = [0, 0, 0];
    A = [];
    Ap = [];
    Bp = [];
    Cp = [];
end

function [metric, residual, J] = amplitude_condition(~, metric, u)
    coeffs = reshape(u(1:end-1), metric.harmonics, []);
    weighted = (coeffs * metric.mass) / 10^(2 * u(end));
    ratio = sum(coeffs .* weighted, 'all');
    residual = ratio - 1;
    J = [2 * weighted(:).', -2 * log(10) * ratio];
end
