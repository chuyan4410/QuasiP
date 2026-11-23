%% Full-order time-integration validation of selected rotating-beam pSSM FRC points
clearvars
clc

exampleDir = fileparts(mfilename('fullpath'));
repoRoot = fileparts(fileparts(fileparts(fileparts(exampleDir))));
rotBeamPssmDir = fullfile(repoRoot,'examples','ParSSM', ...
    'TaylorParm_NL_iner','rotBeam_GEF');
addpath(fullfile(repoRoot,'ext','GEF_quat','src'),'-begin');
addpath(rotBeamPssmDir,'-begin');
clear genParam

%% User settings
frcFile = fullfile(exampleDir,'FRCdata_rotBeam1.mat');

% Frequencies are matched to the nearest saved stable FRC points. Replace
% these values with the desired excitation frequencies [rad/s]. Set this
% vector empty and fill selectedFrcIndices to select columns directly.
targetFrequencies = 150;
selectedFrcIndices = [];

settings.stepsPerPeriod = 100;
settings.initialPeriods = 50;
settings.extensionPeriods = 10;
settings.maximumPeriods = 200;
settings.savedPeriods = 3;

% The final three periods provide two consecutive period-pair checks.
settings.requiredConsecutivePeriodPairs = 2;
settings.displacementTolerance = 1e-5;
settings.velocityTolerance = 1e-5;
settings.amplitudeTolerance = 1e-4;
settings.amplitudeAbsoluteScale = 1e-10;

settings.rho_s = 0.9;
settings.newtonTolerance = 1e-6;
settings.maximumNewtonIterations = 100;
settings.suppressTimeStepOutput = true;
settings.requireStablePssmPoint = true;
settings.projectInitialQuaternionConstraints = true;
settings.overwriteExistingResults = false;
settings.continueAfterFailure = true;

resultsDirectory = fullfile(exampleDir,'full-order-results');
directoryOverride = getenv('ROTBEAM_FRC_OUTPUT_DIR');
if ~isempty(directoryOverride)
    resultsDirectory = directoryOverride;
end

%% Load the saved pSSM periodic-orbit initial conditions
loaded = load(frcFile,'FRC');
assert(isfield(loaded,'FRC') && isstruct(loaded.FRC), ...
    'The selected MAT-file does not contain an FRC structure.');
FRC = loaded.FRC;

requiredFields = {'labels','frequency','parameters','stable','Z0_frc', ...
    'dof_sys','forcingPhase'};
for k = 1:numel(requiredFields)
    assert(isfield(FRC,requiredFields{k}), ...
        'FRC is missing the field "%s".',requiredFields{k});
end

assert(size(FRC.parameters,1)==3, ...
    'This script expects parameters [omegaR;epsilon;Omega].');
assert(size(FRC.Z0_frc,1)==2*FRC.dof_sys, ...
    'FRC.Z0_frc must contain [q;dq] in each column.');
assert(size(FRC.Z0_frc,2)==numel(FRC.frequency), ...
    'FRC.Z0_frc and FRC.frequency contain different numbers of points.');

param = genParam();
expectedDof = 7*(2*param.elmnum+1);
assert(FRC.dof_sys==expectedDof, ...
    ['FRC.dof_sys=%d does not match the rotBeam_GEF model ' ...
     'dimension %d.'],FRC.dof_sys,expectedDof);

% Optional bounded smoke-test settings used for interface verification.
if strcmp(getenv('ROTBEAM_FRC_SMOKE'),'1')
    stableIndices = find(FRC.stable);
    selectedFrcIndices = stableIndices(round((numel(stableIndices)+1)/2));
    targetFrequencies = [];
    settings.stepsPerPeriod = 20;
    settings.initialPeriods = 2;
    settings.extensionPeriods = 2;
    settings.maximumPeriods = 4;
    settings.savedPeriods = 2;
    settings.requiredConsecutivePeriodPairs = 1;
    settings.overwriteExistingResults = true;
end

validate_settings(settings);
selectedFrcIndices = resolve_frc_indices( ...
    FRC,targetFrequencies,selectedFrcIndices, ...
    settings.requireStablePssmPoint);

if ~isfolder(resultsDirectory)
    mkdir(resultsDirectory);
end

settings.sourceFrcFile = frcFile;
settings.resultsDirectory = resultsDirectory;
settings.targetFrequencies = targetFrequencies;
settings.selectedFrcIndices = selectedFrcIndices;
settings.modelParameters = param;

fprintf('Loaded %d pSSM FRC points from:\n  %s\n', ...
    numel(FRC.frequency),frcFile);
fprintf('Selected FRC indices: %s\n',mat2str(selectedFrcIndices));
fprintf(['Time integration: %d steps/period, %d initial periods, ' ...
    '+%d periods/check, maximum %d periods.\n'], ...
    settings.stepsPerPeriod,settings.initialPeriods, ...
    settings.extensionPeriods,settings.maximumPeriods);

%% Fixed Rayleigh damping used by the pSSM model
% Assembly is used once here only to recover the reference M and K needed
% for the identical Rayleigh matrix. Time stepping itself is performed by
% Solver.G_alpha and its element-level assembly; no PDIM Assembly force
% evaluator is called inside the frequency or time-step loops.
[Damp,dampingInfo] = reference_rayleigh_damping(param);
settings.dampingInfo = dampingInfo;

%% Integrate the selected frequencies one by one and save immediately
nSelected = numel(selectedFrcIndices);
summary = struct();
summary.createdAt = char(datetime('now'));
summary.sourceFrcFile = frcFile;
summary.settings = settings;
summary.points = cell(nSelected,1);
summaryFile = fullfile(resultsDirectory,'RotBeam_FOM_summary.mat');

for pointNumber = 1:nSelected
    frcIndex = selectedFrcIndices(pointNumber);
    omega = FRC.frequency(frcIndex);
    label = FRC.labels(frcIndex);
    omegaTag = strrep(sprintf('%.6f',omega),'.','p');
    pointFile = fullfile(resultsDirectory,sprintf( ...
        'RotBeam_FOM_idx%03d_lab%04d_Omega%s.mat', ...
        frcIndex,round(label),omegaTag));

    fprintf('\n============================================================\n');
    fprintf('FRC point %d/%d: index=%d, label=%d, Omega=%.10g rad/s\n', ...
        pointNumber,nSelected,frcIndex,round(label),omega);
    fprintf('Result file:\n  %s\n',pointFile);

    if isfile(pointFile) && ~settings.overwriteExistingResults
        fprintf('Existing result retained; skipping this point.\n');
        existing = load(pointFile,'result');
        if isfield(existing,'result')
            summary.points{pointNumber} = summarize_result(existing.result);
        end
        save(summaryFile,'summary','-v7.3');
        continue
    end

    try
        result = integrate_frc_point( ...
            FRC,frcIndex,param,Damp,settings);
    catch ME
        result = failed_result(FRC,frcIndex,ME);
        warning('RotBeamFRC:PointFailed', ...
            'FRC index %d failed: %s',frcIndex,ME.message);
    end

    result.savedAt = char(datetime('now'));
    save(pointFile,'result','settings','-v7.3');
    summary.points{pointNumber} = summarize_result(result);
    save(summaryFile,'summary','-v7.3');

    if isfield(result,'failed') && result.failed && ...
            ~settings.continueAfterFailure
        error('RotBeamFRC:PointFailed', ...
            'Stopped after failure at FRC index %d.',frcIndex);
    end
end

fprintf('\nFull-order validation loop finished. Summary saved to:\n  %s\n', ...
    summaryFile);

%% Local functions
function validate_settings(settings)
validateattributes(settings.stepsPerPeriod,{'numeric'}, ...
    {'scalar','integer','>=',4});
validateattributes(settings.initialPeriods,{'numeric'}, ...
    {'scalar','integer','positive'});
validateattributes(settings.extensionPeriods,{'numeric'}, ...
    {'scalar','integer','positive'});
validateattributes(settings.maximumPeriods,{'numeric'}, ...
    {'scalar','integer','>=',settings.initialPeriods});
validateattributes(settings.savedPeriods,{'numeric'}, ...
    {'scalar','integer','positive'});
validateattributes(settings.requiredConsecutivePeriodPairs,{'numeric'}, ...
    {'scalar','integer','positive'});
assert(settings.savedPeriods>= ...
    settings.requiredConsecutivePeriodPairs+1, ...
    ['savedPeriods must be at least ' ...
     'requiredConsecutivePeriodPairs+1.']);
assert(settings.initialPeriods>= ...
    settings.requiredConsecutivePeriodPairs+1, ...
    ['initialPeriods must contain enough periods for the requested ' ...
     'consecutive periodicity checks.']);
validateattributes(settings.displacementTolerance,{'numeric'}, ...
    {'scalar','real','positive','finite'});
validateattributes(settings.velocityTolerance,{'numeric'}, ...
    {'scalar','real','positive','finite'});
validateattributes(settings.amplitudeTolerance,{'numeric'}, ...
    {'scalar','real','positive','finite'});
validateattributes(settings.rho_s,{'numeric'}, ...
    {'scalar','real','>=',0,'<=',1});
end

function indices = resolve_frc_indices( ...
        FRC,targetFrequencies,selectedFrcIndices,requireStable)
if ~isempty(selectedFrcIndices)
    validateattributes(selectedFrcIndices,{'numeric'}, ...
        {'vector','integer','positive','<=',numel(FRC.frequency)});
    indices = selectedFrcIndices(:).';
else
    validateattributes(targetFrequencies,{'numeric'}, ...
        {'vector','real','positive','finite','nonempty'});
    if requireStable
        candidates = find(FRC.stable);
    else
        candidates = 1:numel(FRC.frequency);
    end
    assert(~isempty(candidates),'No eligible FRC points were found.');

    indices = zeros(1,numel(targetFrequencies));
    for k = 1:numel(targetFrequencies)
        [mismatch,localIndex] = min(abs( ...
            FRC.frequency(candidates)-targetFrequencies(k)));
        indices(k) = candidates(localIndex);
        fprintf(['Target Omega %.10g matched to FRC index %d at ' ...
            '%.10g rad/s (difference %.3e).\n'], ...
            targetFrequencies(k),indices(k), ...
            FRC.frequency(indices(k)),mismatch);
    end
end

indices = unique(indices,'stable');
if requireStable
    assert(all(FRC.stable(indices)), ...
        'At least one selected pSSM point is not stable.');
end
end

function result = integrate_frc_point(FRC,frcIndex,param,Damp,settings)
dof_sys = FRC.dof_sys;
qInitial = FRC.Z0_frc(1:dof_sys,frcIndex);
dqInitial = FRC.Z0_frc(dof_sys+1:2*dof_sys,frcIndex);
rawConstraintInfo = constraint_diagnostics( ...
    qInitial,dqInitial,param.elmnum);
if settings.projectInitialQuaternionConstraints
    [qInitial,dqInitial] = project_initial_constraints( ...
        qInitial,dqInitial,param.elmnum);
end
qIntegrationInitial = qInitial;
dqIntegrationInitial = dqInitial;
ddqInitial = zeros(dof_sys,1);
aInitial = zeros(dof_sys,1);

rotationSpeed = FRC.parameters(1,frcIndex);
epsilon = FRC.parameters(2,frcIndex);
omega = FRC.parameters(3,frcIndex);
assert(abs(omega-FRC.frequency(frcIndex))<= ...
    1e-10*max(1,abs(omega)), ...
    'The stored parameter frequency does not match FRC.frequency.');

forcingPhase = FRC.forcingPhase(frcIndex);
period = 2*pi/omega;
dt = period/settings.stepsPerPeriod;
outputDof = dof_sys-5;

cyclesIntegrated = 0;
converged = false;
history = struct('cycles',{},'displacementError',{}, ...
    'velocityError',{},'amplitudeError',{},'tipAmplitude',{}, ...
    'meanNewtonIterations',{},'maximumNewtonIterations',{}, ...
    'restartDisplacementMismatch',{},'restartVelocityMismatch',{}, ...
    'restartAccelerationMismatch',{},'blockRuntime',{});
lastSolution = struct();

integrationConstraintInfo = constraint_diagnostics( ...
    qInitial,dqInitial,param.elmnum);

fprintf(['omegaR=%.10g rad/s, epsilon=%.6g, force amplitude=%.6g N, ' ...
    'period=%.6g s, dt=%.6g s.\n'], ...
    rotationSpeed,epsilon,epsilon*param.forcAmp,period,dt);
fprintf(['Raw pSSM quaternion norm residual %.3e; ' ...
    'tangent residual %.3e.\n'], ...
    rawConstraintInfo.quaternionNormResidual, ...
    rawConstraintInfo.quaternionTangentResidual);
if settings.projectInitialQuaternionConstraints
    fprintf(['Projected integration initial residuals: norm %.3e; ' ...
        'tangent %.3e.\n'], ...
        integrationConstraintInfo.quaternionNormResidual, ...
        integrationConstraintInfo.quaternionTangentResidual);
end

while cyclesIntegrated<settings.maximumPeriods
    if cyclesIntegrated==0
        blockPeriods = settings.initialPeriods;
    else
        blockPeriods = settings.extensionPeriods;
    end
    blockPeriods = min(blockPeriods, ...
        settings.maximumPeriods-cyclesIntegrated);

    nBlockSteps = blockPeriods*settings.stepsPerPeriod;
    startTime = cyclesIntegrated*period;
    time = startTime+(0:nBlockSteps)*dt;

    fprintf('Integrating periods %d--%d (%d time steps)...\n', ...
        cyclesIntegrated+1,cyclesIntegrated+blockPeriods,nBlockSteps);

    qBlockInitial = qInitial;
    dqBlockInitial = dqInitial;
    ddqBlockInitial = ddqInitial;
    aBlockInitial = aInitial;

    [nodes,elements,constraints,rotationFrame] = build_beam_block( ...
        param,time,rotationSpeed,omega,epsilon,forcingPhase, ...
        qInitial,dqInitial,ddqInitial);

    solver = Solver.G_alpha( ...
        elements,nodes,constraints,dt,dof_sys);
    solver.rotationFrame = rotationFrame;
    solver.rho_s = settings.rho_s;
    solver.tol = settings.newtonTolerance;
    solver.max_it = settings.maximumNewtonIterations;
    solver.Damp = Damp;
    solver.initialAuxiliaryAcceleration = aInitial;

    blockTimer = tic;
    if settings.suppressTimeStepOutput
        solverLog = evalc('solver.Integrate(time);'); %#ok<NASGU>
    else
        solver.Integrate(time);
    end
    blockRuntime = toc(blockTimer);
    cyclesIntegrated = cyclesIntegrated+blockPeriods;

    restartDisplacementMismatch = norm( ...
        solver.Solution.q_sys(:,1)-qBlockInitial);
    restartVelocityMismatch = norm( ...
        solver.Solution.dq_sys(:,1)-dqBlockInitial);
    restartAccelerationMismatch = max( ...
        norm(solver.Solution.ddq_sys(:,1)-ddqBlockInitial), ...
        norm(solver.Solution.a(:,1)-aBlockInitial));

    metrics = periodicity_metrics( ...
        solver.Solution,settings.stepsPerPeriod,period,outputDof, ...
        settings.requiredConsecutivePeriodPairs,settings);

    history(end+1) = struct( ...
        'cycles',cyclesIntegrated, ...
        'displacementError',metrics.displacementError, ...
        'velocityError',metrics.velocityError, ...
        'amplitudeError',metrics.amplitudeError, ...
        'tipAmplitude',metrics.tipAmplitude, ...
        'meanNewtonIterations',mean(solver.itNum), ...
        'maximumNewtonIterations',max(solver.itNum), ...
        'restartDisplacementMismatch',restartDisplacementMismatch, ...
        'restartVelocityMismatch',restartVelocityMismatch, ...
        'restartAccelerationMismatch',restartAccelerationMismatch, ...
        'blockRuntime',blockRuntime); %#ok<AGROW>

    fprintf(['Completed %d periods in %.2f s: eq=%.3e, ' ...
        'edq=%.3e, eAmp=%.3e, tip amplitude=%.6g m.\n'], ...
        cyclesIntegrated,blockRuntime,metrics.displacementError, ...
        metrics.velocityError,metrics.amplitudeError, ...
        metrics.tipAmplitude);
    fprintf('Restart transfer mismatches: q=%.3e, dq=%.3e, acc=%.3e.\n', ...
        restartDisplacementMismatch,restartVelocityMismatch, ...
        restartAccelerationMismatch);

    lastSolution = retain_last_periods( ...
        solver.Solution,settings.stepsPerPeriod,settings.savedPeriods);

    qInitial = solver.Solution.q_sys(:,end);
    dqInitial = solver.Solution.dq_sys(:,end);
    ddqInitial = solver.Solution.ddq_sys(:,end);
    aInitial = solver.Solution.a(:,end);

    converged = metrics.hasEnoughPeriods && ...
        metrics.displacementError<=settings.displacementTolerance && ...
        metrics.velocityError<=settings.velocityTolerance && ...
        metrics.amplitudeError<=settings.amplitudeTolerance;
    if converged
        fprintf('Periodic convergence achieved after %d periods.\n', ...
            cyclesIntegrated);
        break
    end
end

if ~converged
    fprintf(['Periodic convergence was not achieved within the maximum ' ...
        '%d periods.\n'],settings.maximumPeriods);
end

finalConstraintInfo = constraint_diagnostics( ...
    qInitial,dqInitial,param.elmnum);

result = struct();
result.failed = false;
result.converged = converged;
result.frcIndex = frcIndex;
result.frcLabel = FRC.labels(frcIndex);
result.pssmStable = FRC.stable(frcIndex);
result.frequency = omega;
result.rotationSpeed = rotationSpeed;
result.rotationAcceleration = 0;
result.epsilon = epsilon;
result.forceAmplitude = epsilon*param.forcAmp;
result.forcingPhase = forcingPhase;
result.period = period;
result.timeStep = dt;
result.stepsPerPeriod = settings.stepsPerPeriod;
result.cyclesIntegrated = cyclesIntegrated;
result.outputDof = outputDof;
result.tipAmplitude = history(end).tipAmplitude;
result.periodicityHistory = history;
result.pssmConstraintInfoRaw = rawConstraintInfo;
result.initialConstraintInfo = integrationConstraintInfo;
result.finalConstraintInfo = finalConstraintInfo;
result.pssmInitialCondition = FRC.Z0_frc(:,frcIndex);
result.integrationInitialCondition = [ ...
    qIntegrationInitial;dqIntegrationInitial];
result.finalRestartState = struct( ...
    'q',qInitial,'dq',dqInitial,'ddq',ddqInitial,'a',aInitial);
result.lastPeriods = lastSolution;
end

function [nodes,elements,constraints,rotationFrame] = build_beam_block( ...
        param,time,rotationSpeed,omega,epsilon,forcingPhase, ...
        qInitial,dqInitial,ddqInitial)
timenum = numel(time);
elmnum = param.elmnum;
nodenum = 2*elmnum+1;
dof_sys = 7*nodenum;
elmlen = param.len/elmnum;

material = Material.IsotropicLE(param.E,param.nu,param.rho);
section = CrossSection.Rod3D(param.Area,param.I1,param.I2, ...
    param.I3,param.shearFactor);

AMList = repmat({{zeros(3,1),zeros(3,3)}},timenum,1);
gList = repmat({{zeros(3,1)}},timenum,1);
fList = repmat({{zeros(6,1)}},timenum,1);

nodes = cell(nodenum,1);
qNode = zeros(7,1);
dqNode = zeros(7,1);
ddqNode = zeros(7,1);
for n = 1:nodenum
    dofID = transpose(7*n-6:7*n);
    qNode(:) = 0;
    qNode(1) = elmlen/2*(n-1);
    qNode(4) = 1;
    nodes{n} = Node.HM3DCurvRodNode_quat( ...
        dofID,material,section,param.RM,qNode,dqNode,ddqNode, ...
        timenum,AMList,gList,fList,param.KI);
end

elements = cell(elmnum,1);
for e = 1:elmnum
    elementNodes = {nodes{2*e-1};nodes{2*e};nodes{2*e+1}};
    elements{e} = Element.HM3DCurvRodL3_quat(elementNodes,1);
end

dof_constraint = 6;
constrainedDOF = {{1,[1,2,3,5,6,7]}};
dispSet = zeros(dof_constraint,1);
CqSet = zeros(dof_constraint,dof_sys);
CqSet(1:3,1:3) = eye(3);
CqSet(4:6,5:7) = eye(3);
KbSet = zeros(dof_sys);
dispList = repmat({{dispSet,CqSet,KbSet}},timenum,1);
displacementConstraint = Constraint.Constraint( ...
    constrainedDOF,{{'displacement boundary'};{[]}}, ...
    dispList,dof_constraint);

forceList = cell(timenum,1);
forceAmplitude = epsilon*param.forcAmp;
for k = 1:timenum
    forceValue = forceAmplitude*cos(omega*time(k)+forcingPhase);
    forceList{k} = {[0;forceValue;0;0;0;0]};
end
dofApplied = {{nodenum,2}};
forceConstraint = Constraint.Constraint( ...
    dofApplied,{{'force boundary'}},forceList,dof_constraint);

constraints.Dconstraints = displacementConstraint;
constraints.Fconstraints = forceConstraint;

set_node_initial_state(nodes,qInitial,dqInitial,ddqInitial);
rotationFrame = constant_speed_rotation_frame( ...
    param.rotationAxis,rotationSpeed,time);
end

function set_node_initial_state(nodes,q,dq,ddq)
for n = 1:numel(nodes)
    id = nodes{n}.dofID;
    nodes{n}.q(:,1) = q(id);
    nodes{n}.dq(:,1) = dq(id);
    nodes{n}.ddq(:,1) = ddq(id);
    nodes{n}.updateNode(1,q(id),dq(id),ddq(id));
end
end

function rotationFrame = constant_speed_rotation_frame(axisR,omega,time)
axisR = axisR(:)/norm(axisR);
phi = omega*time(:).';
c = cos(phi/2);
s = sin(phi/2);

e = [c;axisR*s];
de = 0.5*omega*[-s;axisR*c];
dde = -0.25*omega^2*[c;axisR*s];

rotationFrame = struct( ...
    'enabled',true,'e',e,'de',de,'dde',dde);
end

function metrics = periodicity_metrics( ...
        solution,stepsPerPeriod,period,outputDof,requiredPairs,settings)
nTime = size(solution.q_sys,2);
availablePeriods = floor((nTime-1)/stepsPerPeriod);
nPairs = min(requiredPairs,availablePeriods-1);

metrics.hasEnoughPeriods = nPairs>=requiredPairs;
metrics.displacementError = inf;
metrics.velocityError = inf;
metrics.amplitudeError = inf;
metrics.tipAmplitude = max(abs(solution.q_sys(outputDof, ...
    max(1,nTime-stepsPerPeriod):nTime)));

if nPairs<1
    return
end

qErrors = zeros(nPairs,1);
dqErrors = zeros(nPairs,1);
amplitudeErrors = zeros(nPairs,1);
amplitudes = zeros(nPairs,1);

for pairIndex = 1:nPairs
    currentEnd = nTime-(pairIndex-1)*stepsPerPeriod;
    currentIndices = currentEnd-stepsPerPeriod:currentEnd;
    previousIndices = currentEnd-2*stepsPerPeriod: ...
        currentEnd-stepsPerPeriod;

    qCurrent = solution.q_sys(:,currentIndices);
    qPrevious = solution.q_sys(:,previousIndices);
    dqCurrent = solution.dq_sys(:,currentIndices);
    dqPrevious = solution.dq_sys(:,previousIndices);

    qScale = max(1,norm(qCurrent,'fro'));
    qErrors(pairIndex) = ...
        norm(qCurrent-qPrevious,'fro')/qScale;
    dqErrors(pairIndex) = period* ...
        norm(dqCurrent-dqPrevious,'fro')/qScale;

    currentAmplitude = max(abs(qCurrent(outputDof,:)));
    previousAmplitude = max(abs(qPrevious(outputDof,:)));
    amplitudeScale = max([currentAmplitude,previousAmplitude, ...
        settings.amplitudeAbsoluteScale]);
    amplitudeErrors(pairIndex) = ...
        abs(currentAmplitude-previousAmplitude)/amplitudeScale;
    amplitudes(pairIndex) = currentAmplitude;
end

metrics.displacementError = max(qErrors);
metrics.velocityError = max(dqErrors);
metrics.amplitudeError = max(amplitudeErrors);
metrics.tipAmplitude = amplitudes(1);
metrics.pairDisplacementErrors = qErrors;
metrics.pairVelocityErrors = dqErrors;
metrics.pairAmplitudeErrors = amplitudeErrors;
end

function retained = retain_last_periods(solution,stepsPerPeriod,nPeriods)
nTime = size(solution.q_sys,2);
firstIndex = max(1,nTime-nPeriods*stepsPerPeriod);
indices = firstIndex:nTime;

retained = struct();
retained.t = solution.t(indices);
retained.q = solution.q_sys(:,indices);
retained.dq = solution.dq_sys(:,indices);
retained.ddq = solution.ddq_sys(:,indices);
retained.a = solution.a(:,indices);
retained.lambda = solution.lm(:,indices);
end

function info = constraint_diagnostics(q,dq,elmnum)
nodenum = 2*elmnum+1;
normResidual = zeros(nodenum,1);
tangentResidual = zeros(nodenum,1);
for n = 1:nodenum
    id = 7*n-3:7*n;
    quat = q(id);
    dquat = dq(id);
    normResidual(n) = abs(quat.'*quat-1);
    tangentResidual(n) = abs(quat.'*dquat);
end
info.quaternionNormResidual = max(normResidual);
info.quaternionTangentResidual = max(tangentResidual);
info.rootDisplacementResidual = norm(q([1,2,3,5,6,7]));
info.rootVelocityResidual = norm(dq([1,2,3,5,6,7]));
end

function [qProjected,dqProjected] = project_initial_constraints(q,dq,elmnum)
qProjected = q;
dqProjected = dq;

% Enforce the clamped-root coordinates exactly.
rootConstrainedDofs = [1,2,3,5,6,7];
qProjected(rootConstrainedDofs) = 0;
dqProjected(rootConstrainedDofs) = 0;

nodenum = 2*elmnum+1;
for n = 1:nodenum
    id = 7*n-3:7*n;
    quat = qProjected(id);
    quatNorm = norm(quat);
    assert(quatNorm>eps, ...
        'A nodal quaternion in the pSSM initial state has zero norm.');
    quat = quat/quatNorm;
    dquat = dqProjected(id);
    dquat = dquat-quat*(quat.'*dquat);
    qProjected(id) = quat;
    dqProjected(id) = dquat;
end
end

function [Damp,info] = reference_rayleigh_damping(param)
nodenum = 2*param.elmnum+1;
dof_sys = 7*nodenum;
[qReference,~,~] = reference_state(param);
zerosState = zeros(dof_sys,1);

[nodes,elements,constraints,~] = build_beam_block( ...
    param,[0,1],0,1,0,0,qReference,zerosState,zerosState);
assembly = Assembly.Assembly( ...
    elements,nodes,constraints,qReference,zerosState,zerosState);
assembly.genSysMatVec(2);

M = assembly.M;
K = assembly.Kt;
U = null(full(assembly.Cq));
Mhat = U.'*M*U;
Khat = U.'*K*U;
nEig = min(10,size(Mhat,1)-1);
[~,D2] = eigs(sparse(Khat),sparse(Mhat),nEig,'smallestabs');
frequencies = sort(real(sqrt(diag(D2))));
frequencies = frequencies( ...
    isfinite(frequencies) & frequencies>0);
assert(numel(frequencies)>=2, ...
    'At least two elastic modes are required for Rayleigh damping.');

omega1 = frequencies(1);
omega2 = frequencies(2);
alphaPerRatio = ...
    2*omega1*omega2*(omega2-omega1)/(omega2^2-omega1^2);
betaPerRatio = ...
    2*(omega2-omega1)/(omega2^2-omega1^2);
Damp = param.dampRatio*(alphaPerRatio*M+betaPerRatio*K);

info = struct();
info.fixedDampingRatio = param.dampRatio;
info.alpha = param.dampRatio*alphaPerRatio;
info.beta = param.dampRatio*betaPerRatio;
info.referenceFrequencies = frequencies;
end

function [qReference,dqReference,ddqReference] = reference_state(param)
nodenum = 2*param.elmnum+1;
dof_sys = 7*nodenum;
elmlen = param.len/param.elmnum;
qReference = zeros(dof_sys,1);
for n = 1:nodenum
    id = 7*n-6:7*n;
    qReference(id(1)) = elmlen/2*(n-1);
    qReference(id(4)) = 1;
end
dqReference = zeros(dof_sys,1);
ddqReference = zeros(dof_sys,1);
end

function compact = summarize_result(result)
fields = {'failed','converged','frcIndex','frcLabel','frequency', ...
    'rotationSpeed','epsilon','forceAmplitude','cyclesIntegrated', ...
    'tipAmplitude','savedAt'};
compact = struct();
for k = 1:numel(fields)
    if isfield(result,fields{k})
        compact.(fields{k}) = result.(fields{k});
    end
end
if isfield(result,'failureMessage')
    compact.failureMessage = result.failureMessage;
end
end

function result = failed_result(FRC,frcIndex,ME)
result = struct();
result.failed = true;
result.converged = false;
result.frcIndex = frcIndex;
result.frcLabel = FRC.labels(frcIndex);
result.frequency = FRC.frequency(frcIndex);
result.rotationSpeed = FRC.parameters(1,frcIndex);
result.epsilon = FRC.parameters(2,frcIndex);
result.failureIdentifier = ME.identifier;
result.failureMessage = ME.message;
result.failureReport = getReport(ME,'extended','hyperlinks','off');
end
