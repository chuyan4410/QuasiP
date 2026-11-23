%% Save Aditya Case-II numerical results without class objects

saveFile = fullfile(exampleDir, 'Aditya-Case-II-results(new14).mat');

%% Case/config data
caseInfo = struct();
caseInfo.name = 'Jeong NREL 5MW aeroelastic reproduction';
caseInfo.savedAt = char(datetime('now'));
caseInfo.script = mfilename('fullpath');

caseInfo.dt = dt;
caseInfo.endTime = endTime;
caseInfo.time = t(:);
caseInfo.gravity = gravity;
caseInfo.hubWind = hubWind;
caseInfo.rotSpeedRPM = rotSpeedRPM;
caseInfo.OmegaHub = OmegaHub;
caseInfo.pitchDeg = pitchDeg;
caseInfo.loadRampTime = loadRampTime;
caseInfo.elmnum = elmnum;
caseInfo.nodenum = nodenum;
caseInfo.dof_sys = dof_sys;
caseInfo.reduce = reduce;

caseInfo.blade = blade;
caseInfo.radius = radius(:);
caseInfo.aeroData = aeroData;
caseInfo.aeroOptions = aeroOptions;

%% Solver settings
solverInfo = struct();
solverInfo.rho_s = solver.rho_s;
solverInfo.alpha_f = solver.alpha_f;
solverInfo.alpha_m = solver.alpha_m;
solverInfo.gamma = solver.gamma;
solverInfo.beta = solver.beta;
solverInfo.gamma_hat = solver.gamma_hat;
solverInfo.beta_hat = solver.beta_hat;
solverInfo.tol = solver.tol;
solverInfo.max_it = solver.max_it;
solverInfo.dampRatio = solver.dampRatio;
solverInfo.aeroForceUpdateMode = solver.aeroForceUpdateMode;
solverInfo.itNum = solver.itNum;

solverInfo.rotFrame = rotFrame;
if isfield(rotFrame, 'e')
    solverInfo.rotFrame.azimuth = unwrap(2 * atan2(rotFrame.e(4,:), rotFrame.e(1,:))).';
end

%% Structural time histories
result = struct();
result.t = solver.Solution.t(:);
result.q = solver.Solution.q_sys;
result.dq = solver.Solution.dq_sys;
result.ddq = solver.Solution.ddq_sys;
result.a = solver.Solution.a;
result.lm = solver.Solution.lm;

tipNode = nodenum;
tipDof = (tipNode - 1) * 7 + (1:7);
result.tip.q = result.q(tipDof, :);
result.tip.dq = result.dq(tipDof, :);
result.tip.ddq = result.ddq(tipDof, :);
result.tip.position = result.q(tipDof(1:3), :);
result.tip.velocity = result.dq(tipDof(1:3), :);
result.tip.acceleration = result.ddq(tipDof(1:3), :);
result.tip.quaternion = result.q(tipDof(4:7), :);

%% Rotor integrated loads
[thrust, torque] = rotorAeroLoadsFromHistory(solver.aeroForceHistory, radius, 3);
result.rotor.thrust = thrust(:);
result.rotor.torque = torque(:);
result.rotor.power = torque(:) * OmegaHub;

%% Aero histories: station and Gauss-point outputs
aeroHistory = solver.aeroForceHistory;
result.aero = packAeroHistory(aeroHistory);

result.aero.station.alphaDeg = rad2deg(result.aero.station.alpha);
result.aero.station.phiDeg = rad2deg(result.aero.station.phi);
result.aero.station.thetaElasticDeg = rad2deg(result.aero.station.thetaElastic);

result.aero.gauss.alphaDeg = rad2deg(result.aero.gauss.alpha);
result.aero.gauss.phiDeg = rad2deg(result.aero.gauss.phi);
result.aero.gauss.thetaElasticDeg = rad2deg(result.aero.gauss.thetaElastic);

%% Save
save(saveFile, 'caseInfo', 'solverInfo', 'result', '-v7.3');
fprintf('Saved numerical results to:\n  %s\n', saveFile);

%% Local helper
function [thrust, torque] = rotorAeroLoadsFromHistory(aeroHistory, radius, bladeCount)
    thrust = zeros(numel(aeroHistory), 1);
    torque = zeros(numel(aeroHistory), 1);

    valid = find(~cellfun(@isempty, aeroHistory), 1, 'first');
    if isempty(valid)
        return
    end

    [weights, ~] = Base.gaussQuadrature(3);
    elemnum = (numel(radius) - 1) / 2;
    for step = 1:numel(aeroHistory)
        if isempty(aeroHistory{step})
            continue
        end

        g = aeroHistory{step}.gauss;
        scale = aeroHistory{step}.loadScale;
        bladeThrust = 0;
        bladeTorque = 0;
        for e = 1:elemnum
            jac = (radius(2*e+1) - radius(2*e-1)) / 2;
            ids = find(g.elementID == e);
            for k = 1:numel(ids)
                id = ids(k);
                w = weights(g.gaussID(id));
                bladeThrust = bladeThrust + scale * g.normalForce(id) * w * jac;
                % Positive tangential force follows the rotating-frame +Y direction.
                bladeTorque = bladeTorque + scale * g.tangentialForce(id) * g.radius(id) * w * jac;
            end
        end
        thrust(step) = bladeCount * bladeThrust;
        torque(step) = bladeCount * bladeTorque;
    end
end

function aero = packAeroHistory(aeroHistory)
    valid = find(~cellfun(@isempty, aeroHistory));
    if isempty(valid)
        aero = struct();
        return
    end

    first = aeroHistory{valid(1)};
    pointSets = {'station', 'gauss'};

    aero = struct();
    aero.validStep = valid(:);

    for p = 1:numel(pointSets)
        setName = pointSets{p};
        if ~isfield(first, setName)
            continue
        end

        fields = fieldnames(first.(setName));
        nPoint = numel(first.(setName).radius);
        nStep = numel(aeroHistory);

        for f = 1:numel(fields)
            field = fields{f};
            sample = first.(setName).(field);

            if isnumeric(sample) && isvector(sample)
                aero.(setName).(field) = nan(nPoint, nStep);
            end
        end

        for s = valid(:).'
            data = aeroHistory{s}.(setName);
            for f = 1:numel(fields)
                field = fields{f};
                if isfield(aero.(setName), field)
                    value = data.(field);
                    aero.(setName).(field)(1:numel(value), s) = value(:);
                end
            end
        end
    end

    aero.loadScale = nan(1, numel(aeroHistory));
    aero.priming = false(1, numel(aeroHistory));
    for s = valid(:).'
        aero.loadScale(s) = aeroHistory{s}.loadScale;
        aero.priming(s) = aeroHistory{s}.priming;
    end
end
