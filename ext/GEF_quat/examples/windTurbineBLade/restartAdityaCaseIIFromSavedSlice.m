%% Restart Aditya Case-II from a saved result time slice
% This script reloads the pure-struct result file written by saveResults.m,
% rebuilds the beam/aero model with the original reference configuration,
% and starts a new time integration from any saved time slice.

clear
clc

%% User settings
thisFile = mfilename('fullpath');
exampleDir = fileparts(thisFile);
repoRoot = fileparts(fileparts(exampleDir));
addpath(genpath(fullfile(repoRoot, 'src')));
addpath(genpath(fullfile(repoRoot, 'ext')));

restartFile = fullfile(exampleDir, 'Aditya-Case-II-results(new11).mat');
% restartFile = fullfile(exampleDir, 'Jeong-results(new8).mat');

% Pick either restartIndex or restartTime. If restartTime is non-empty, the
% nearest saved time slice is used.
restartIndex = 4400;
restartTime = [];

% Number of new integration increments. Use 1 to recompute only the next
% step, which is usually the fastest way to debug a failed time slice.
numStepsToRun = 1;

saveRestartResult = false;

%% Load saved pure data
data = load(restartFile, 'caseInfo', 'solverInfo', 'result');
caseInfo = data.caseInfo;
solverInfo = data.solverInfo;
result = data.result;

savedTime = result.t(:);
restartIndex = selectRestartIndex(savedTime, restartIndex, restartTime);
restartTimeAbs = savedTime(restartIndex);

if numStepsToRun < 1 || fix(numStepsToRun) ~= numStepsToRun
    error('restart:numStepsToRun', 'numStepsToRun must be a positive integer.');
end

dt = caseInfo.dt;
timeLocal = restartTimeAbs + (0:numStepsToRun).' * dt;
timenum = numel(timeLocal);

fprintf('\nRestart source file:\n  %s\n', restartFile);
fprintf('Restarting from saved index %d at t = %.6f s\n', ...
    restartIndex, restartTimeAbs);
fprintf('Running %d new increment(s), final local t = %.6f s\n', ...
    numStepsToRun, timeLocal(end));

%% Rebuild the model with the original reference configuration
referenceQ = result.q(:, 1);
restartQ = result.q(:, restartIndex);
restartDQ = result.dq(:, restartIndex);
restartDDQ = result.ddq(:, restartIndex);

[myNodes, myElements, myConstraints, modelInfo] = buildRestartBeamModel( ...
    caseInfo, referenceQ, timenum);
setNodeInitialState(myNodes, restartQ, restartDQ, restartDDQ);

rotFrame = restartRotationFrame(solverInfo, caseInfo, restartIndex, timeLocal);
aeroOptions = restartAeroOptions(caseInfo, restartIndex, timeLocal);
aeroModel = Aero.AdityaBEMT(caseInfo.aeroData, aeroOptions);
% aeroModel = Aero.JeongBEMT(caseInfo.aeroData, aeroOptions);

%% Time integration from the selected slice
solver = Solver.G_alpha(myElements, myNodes, myConstraints, dt, modelInfo.dof_sys);
solver.rho_s = solverInfo.rho_s;
solver.rotationFrame = rotFrame;
solver.aeroForceModel = aeroModel;
solver.aeroForceUpdateMode = solverInfo.aeroForceUpdateMode;
solver.tol = solverInfo.tol;
solver.max_it = solverInfo.max_it;
solver.dampRatio = solverInfo.dampRatio;
solver.Damp = zeros(modelInfo.dof_sys);

if isfield(result, 'a') && size(result.a, 2) >= restartIndex
    solver.initialAuxiliaryAcceleration = result.a(:, restartIndex);
end

tic
solver.Integrate(timeLocal);
toc

%% Pack lightweight output for inspection
[thrust, torque] = rotorAeroLoadsFromHistory( ...
    solver.aeroForceHistory, caseInfo.radius, aeroOptions.bladeCount, aeroOptions);

restartInfo = struct();
restartInfo.sourceFile = restartFile;
restartInfo.sourceIndex = restartIndex;
restartInfo.sourceTime = restartTimeAbs;
restartInfo.numStepsToRun = numStepsToRun;
restartInfo.savedAt = char(datetime('now'));
restartInfo.caseInfo = caseInfo;
restartInfo.solverInfo = solverInfo;

restartResult = struct();
restartResult.t = solver.Solution.t(:);
restartResult.q = solver.Solution.q_sys;
restartResult.dq = solver.Solution.dq_sys;
restartResult.ddq = solver.Solution.ddq_sys;
restartResult.a = solver.Solution.a;
restartResult.lm = solver.Solution.lm;
restartResult.rotor.thrust = thrust(:);
restartResult.rotor.torque = torque(:);
restartResult.rotor.power = torque(:) * caseInfo.OmegaHub;
restartResult.aero = packAeroHistory(solver.aeroForceHistory);

fprintf('\nRestart run finished.\n');
if numel(solver.itNum) > 0
    fprintf('  Newton iterations in recomputed increment(s): %s\n', ...
        mat2str(solver.itNum(:).'));
end
if numel(thrust) >= 2
    fprintf('  Last computed 3-blade thrust = %.6e N\n', thrust(end));
    fprintf('  Last computed 3-blade torque = %.6e N-m\n', torque(end));
    fprintf('  Last computed 3-blade power  = %.6e W\n', torque(end) * caseInfo.OmegaHub);
end

if saveRestartResult
    outputFile = fullfile(exampleDir, sprintf( ...
        'Aditya-Case-II-restart-from-%06d.mat', restartIndex));
    save(outputFile, 'restartInfo', 'restartResult', '-v7.3');
    fprintf('Saved restart result to:\n  %s\n', outputFile);
end

%% Local helpers
function restartIndex = selectRestartIndex(savedTime, restartIndex, restartTime)
    if ~isempty(restartTime)
        [~, restartIndex] = min(abs(savedTime - restartTime));
    end

    nStep = numel(savedTime);
    if restartIndex < 1 || restartIndex > nStep
        error('restart:InvalidIndex', ...
            'restartIndex must be between 1 and %d.', nStep);
    end
end

function [nodes, elements, constraints, modelInfo] = buildRestartBeamModel( ...
    caseInfo, referenceQ, timenum)

    blade = caseInfo.blade;
    radius = caseInfo.radius(:);
    elmnum = caseInfo.elmnum;
    nodenum = caseInfo.nodenum;
    dof_sys = caseInfo.dof_sys;

    if numel(referenceQ) ~= dof_sys
        error('restart:ReferenceStateSize', ...
            'The saved reference state size does not match caseInfo.dof_sys.');
    end

    AMList = repmat({{zeros(3,1), zeros(3,3)}}, timenum, 1);
    gList = repmat({{caseInfo.gravity}}, timenum, 1);
    fList = repmat({{zeros(6,1)}}, timenum, 1);
    RM = [0; 0; 0];

    nodes = cell(nodenum, 1);
    dq0 = zeros(7,1);
    ddq0 = zeros(7,1);
    dummyMaterial = Material.IsotropicLE(0, 0, 0);

    for n = 1:nodenum
        sec = directRodSection(blade, n);
        mySection = CrossSection.DirectRod3D( ...
            sec.rhoA, sec.Jrho, sec.CA_mB, sec.CI_mB);

        dofID = transpose(7 * n - 6 : 7 * n);
        q0 = referenceQ(dofID);

        nodes{n} = Node.HM3DCurvRodNode_quat(dofID, dummyMaterial, mySection, ...
            RM, q0, dq0, ddq0, timenum, AMList, gList, fList, 1);
    end

    elements = cell(elmnum, 1);
    for e = 1:elmnum
        elemNodes = {nodes{2*e-1}; nodes{2*e}; nodes{2*e+1}};
        elements{e} = Element.HM3DCurvRodL3_quat(elemNodes, caseInfo.reduce);
    end

    constraints = buildRootConstraints(nodes, referenceQ, timenum, dof_sys, nodenum);

    modelInfo = struct();
    modelInfo.elmnum = elmnum;
    modelInfo.nodenum = nodenum;
    modelInfo.dof_sys = dof_sys;
    modelInfo.radius = radius;
end

function constraints = buildRootConstraints(nodes, referenceQ, timenum, dof_sys, nodenum)
    dof_constraint = 6;
    constrainedDOF = {{1, [1, 2, 3, 5, 6, 7]}};
    rootPosition0 = referenceQ(1:3);
    rootQuat0 = referenceQ(4:7);
    dispSet = [rootPosition0; rootQuat0(2:4)];

    CqSet = zeros(dof_constraint, dof_sys);
    CqSet(1:3, 1:3) = eye(3);
    CqSet(4:6, 5:7) = eye(3);
    KbSet = zeros(dof_sys, dof_sys);
    dispList = repmat({{dispSet, CqSet, KbSet}}, timenum, 1);
    constrainedType = {{'displacement boundary'}; {[]}};
    myDispConstraint = Constraint.Constraint(constrainedDOF, ...
        constrainedType, dispList, dof_constraint);

    forcList = repmat({{zeros(6,1)}}, timenum, 1);
    dof_applied = {{nodenum, 2}};
    myForcConstraint = Constraint.Constraint(dof_applied, ...
        {{'force boundary'}}, forcList, dof_constraint);

    constraints.Dconstraints = myDispConstraint;
    constraints.Fconstraints = myForcConstraint;

    if ~isequal(nodes{1}.dofID(:), (1:7).')
        error('restart:RootDofLayout', ...
            'The restart script assumes node 1 owns global DOFs 1:7.');
    end
end

function setNodeInitialState(nodes, q, dq, ddq)
    for n = 1:numel(nodes)
        id = nodes{n}.dofID;
        nodes{n}.q(:,1) = q(id);
        nodes{n}.dq(:,1) = dq(id);
        nodes{n}.ddq(:,1) = ddq(id);
        nodes{n}.updateNode(1, q(id), dq(id), ddq(id));
    end
end

function aeroOptions = restartAeroOptions(caseInfo, restartIndex, timeLocal)
    aeroOptions = caseInfo.aeroOptions;
    aeroOptions.time = timeLocal(:);

    if isfield(caseInfo.aeroOptions, 'loadScale')
        loadScale0 = caseInfo.aeroOptions.loadScale;
        lastIndex = restartIndex + numel(timeLocal) - 1;
        if ~isscalar(loadScale0) && numel(loadScale0) >= lastIndex
            aeroOptions.loadScale = loadScale0(restartIndex:lastIndex);
        elseif isfield(caseInfo, 'loadRampTime')
            aeroOptions.loadScale = min(timeLocal(:) / caseInfo.loadRampTime, 1);
        else
            aeroOptions.loadScale = 1;
        end
    end
end

function rotFrame = restartRotationFrame(solverInfo, caseInfo, restartIndex, timeLocal)
    if ~isfield(solverInfo, 'rotFrame') || ~isfield(solverInfo.rotFrame, 'enabled') ...
            || ~solverInfo.rotFrame.enabled
        rotFrame = struct('enabled', false);
        return
    end

    lastIndex = restartIndex + numel(timeLocal) - 1;
    source = solverInfo.rotFrame;
    canSlice = all(isfield(source, {'e', 'de', 'dde'})) ...
        && size(source.e, 2) >= lastIndex ...
        && size(source.de, 2) >= lastIndex ...
        && size(source.dde, 2) >= lastIndex;

    if canSlice
        idx = restartIndex:lastIndex;
        rotFrame = struct('enabled', true, ...
            'e', source.e(:, idx), ...
            'de', source.de(:, idx), ...
            'dde', source.dde(:, idx));
    elseif isfield(caseInfo, 'OmegaHub') && isfield(caseInfo, 'loadRampTime')
        rotFrame = prescribedPrimingRotatingFrame( ...
            timeLocal(:).', caseInfo.OmegaHub, caseInfo.loadRampTime);
    else
        error('restart:RotationFrame', ...
            'Cannot slice or regenerate the rotating frame for the restart run.');
    end
end

function sec = directRodSection(blade, n)
    sec.rhoA = blade.mass(n);
    sec.Jrho = diag([blade.flpIner(n) + blade.edgIner(n); ...
                     blade.flpIner(n); ...
                     blade.edgIner(n)]);
    sec.CA_mB = diag([blade.EA(n); blade.GA(n); blade.GA(n)]);
    sec.CI_mB = diag([blade.GJ(n); blade.flpStff(n); blade.edgStff(n)]);
end

function rotFrame = prescribedPrimingRotatingFrame(t, omegaRated, rampTime)
    e = zeros(4, numel(t));
    de = zeros(4, numel(t));
    dde = zeros(4, numel(t));
    for k = 1:numel(t)
        if t(k) < rampTime
            omega = omegaRated * t(k) / rampTime;
            dotOmega = omegaRated / rampTime;
            phi = 0.5 * omegaRated / rampTime * t(k)^2;
        else
            omega = omegaRated;
            dotOmega = 0;
            phi = 0.5 * omegaRated * rampTime + omegaRated * (t(k) - rampTime);
        end
        e(:,k) = [cos(phi/2); 0; 0; sin(phi/2)];
        de(:,k) = 0.5 * omega * [-sin(phi/2); 0; 0; cos(phi/2)];
        dde(:,k) = -0.25 * omega^2 * [cos(phi/2); 0; 0; sin(phi/2)] ...
            + 0.5 * dotOmega * [-sin(phi/2); 0; 0; cos(phi/2)];
    end
    rotFrame = struct('enabled', true, 'e', e, 'de', de, 'dde', dde);
end

function [thrust, torque] = rotorAeroLoadsFromHistory( ...
    aeroHistory, radius, bladeCount, aeroOptions)

    thrust = zeros(numel(aeroHistory), 1);
    torque = zeros(numel(aeroHistory), 1);

    valid = find(~cellfun(@isempty, aeroHistory), 1, 'first');
    if isempty(valid)
        return
    end

    [weights, ~] = Base.gaussQuadrature(aeroOptions.numGauss);
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
                bladeTorque = bladeTorque - scale * g.tangentialForce(id) ...
                    * g.radius(id) * w * jac;
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
