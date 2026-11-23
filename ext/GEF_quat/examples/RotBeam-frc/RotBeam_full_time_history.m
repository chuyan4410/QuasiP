%% Full time-history integration of one rotating-beam FRC point
% This script performs one uninterrupted full-order time integration.  It
% does not test periodic convergence and does not discard the transient.
% The complete time histories of the selected output DOFs are saved.
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

% targetFrequency is matched to the nearest saved stable FRC point.  Set
% selectedFrcIndex to a column number to select an FRC point directly.
targetFrequency = 150;             % excitation frequency [rad/s]
selectedFrcIndex = [];

% For the present 4-element model, DOF 58 is the tip-node y position.
% Several output DOFs may be supplied, for example: outdof = [57,58,59].
outdof = 58;

settings.numberOfPeriods = 1000;
settings.stepsPerPeriod = 100;
settings.rho_s = 0.9;
settings.newtonTolerance = 1e-6;
settings.maximumNewtonIterations = 100;
% 'reference': straight undeformed beam with zero velocity (default here).
% 'pssm':      reconstructed pSSM periodic-orbit initial condition.
settings.initialConditionSource = 'reference';
settings.projectInitialQuaternionConstraints = true;
settings.suppressTimeStepOutput = true;

% The selected output histories are always saved.  Enabling this option
% additionally saves q, dq, ddq, a and lambda for every system DOF and can
% create a much larger MAT-file.
settings.saveFullSystemHistory = false;
settings.saveFigures = true;

resultsDirectory = fullfile(exampleDir,'full-time-history-results');
directoryOverride = getenv('ROTBEAM_FULL_HISTORY_OUTPUT_DIR');
if ~isempty(directoryOverride)
    resultsDirectory = directoryOverride;
end

% Optional short run used only for implementation checks.
if strcmp(getenv('ROTBEAM_FULL_HISTORY_SMOKE'),'1')
    settings.numberOfPeriods = 2;
    settings.stepsPerPeriod = 20;
end

%% Load and select one pSSM FRC point
loaded = load(frcFile,'FRC');
assert(isfield(loaded,'FRC') && isstruct(loaded.FRC), ...
    'The selected MAT-file does not contain an FRC structure.');
FRC = loaded.FRC;

requiredFields = {'Z0_frc','frequency','parameters','stable','dof_sys'};
assert(all(isfield(FRC,requiredFields)), ...
    'The FRC structure does not contain all required fields.');
assert(size(FRC.Z0_frc,1)==2*FRC.dof_sys, ...
    'FRC.Z0_frc must contain [q;dq] in each column.');

if isempty(selectedFrcIndex)
    stableIndices = find(FRC.stable);
    assert(~isempty(stableIndices),'No stable points are stored in FRC.');
    [~,nearestLocalIndex] = min(abs( ...
        FRC.frequency(stableIndices)-targetFrequency));
    frcIndex = stableIndices(nearestLocalIndex);
else
    validateattributes(selectedFrcIndex,{'numeric'}, ...
        {'scalar','integer','>=',1,'<=',numel(FRC.frequency)});
    frcIndex = selectedFrcIndex;
    assert(FRC.stable(frcIndex), ...
        'The selected FRC point is not marked as stable.');
end

param = genParam();
dof_sys = FRC.dof_sys;
expectedDofs = 7*(2*param.elmnum+1);
assert(dof_sys==expectedDofs, ...
    'FRC.dof_sys is inconsistent with the rotating-beam model.');
validateattributes(outdof,{'numeric'}, ...
    {'vector','integer','>=',1,'<=',dof_sys});
outdof = reshape(outdof,1,[]);

omega = FRC.frequency(frcIndex);
rotationSpeed = FRC.parameters(1,frcIndex);
epsilon = FRC.parameters(2,frcIndex);
forcingPhase = 0;
if isfield(FRC,'forcingPhase')
    forcingPhase = FRC.forcingPhase(frcIndex);
end

period = 2*pi/omega;
dt = period/settings.stepsPerPeriod;
numberOfSteps = settings.numberOfPeriods*settings.stepsPerPeriod;
time = (0:numberOfSteps)*dt;

pssmQInitial = FRC.Z0_frc(1:dof_sys,frcIndex);
pssmDqInitial = FRC.Z0_frc(dof_sys+1:2*dof_sys,frcIndex);
initialConditionSource = validatestring( ...
    settings.initialConditionSource,{'reference','pssm'});
switch initialConditionSource
    case 'reference'
        [qInitial,dqInitial,ddqInitial] = reference_state(param);
    case 'pssm'
        qInitial = pssmQInitial;
        dqInitial = pssmDqInitial;
        ddqInitial = zeros(dof_sys,1);
        if settings.projectInitialQuaternionConstraints
            [qInitial,dqInitial] = project_initial_constraints( ...
                qInitial,dqInitial,param.elmnum);
        end
end
aInitial = zeros(dof_sys,1);

fprintf(['Selected FRC column %d: Omega=%.10g rad/s, ' ...
    'omegaR=%.10g rad/s, epsilon=%.6g.\n'], ...
    frcIndex,omega,rotationSpeed,epsilon);
fprintf('Initial condition source: %s.\n',initialConditionSource);
fprintf(['Integrating %d periods with %d steps per period ' ...
    '(%d time steps).\n'],settings.numberOfPeriods, ...
    settings.stepsPerPeriod,numberOfSteps);

%% Construct and integrate the full-order model
% Assembly is used once here only to construct the same Rayleigh damping
% matrix as the pSSM example.  Time stepping uses Solver.G_alpha and the
% GEF element-level assembly.
[Damp,dampingInfo] = reference_rayleigh_damping(param);
[nodes,elements,constraints,rotationFrame] = build_beam_model( ...
    param,time,rotationSpeed,omega,epsilon,forcingPhase, ...
    qInitial,dqInitial,ddqInitial);

solver = Solver.G_alpha(elements,nodes,constraints,dt,dof_sys);
solver.rotationFrame = rotationFrame;
solver.rho_s = settings.rho_s;
solver.tol = settings.newtonTolerance;
solver.max_it = settings.maximumNewtonIterations;
solver.Damp = Damp;
solver.initialAuxiliaryAcceleration = aInitial;

integrationTimer = tic;
if settings.suppressTimeStepOutput
    solverLog = evalc('solver.Integrate(time);'); %#ok<NASGU>
    clear solverLog
else
    solver.Integrate(time);
end
integrationRuntime = toc(integrationTimer);
fprintf('Time integration completed in %.2f s.\n',integrationRuntime);

%% Retain the complete time history
solution = solver.Solution;
historyTime = reshape(solution.t,1,[]);
cycleNumber = historyTime/period;
periodEndIndices = 1:settings.stepsPerPeriod:numel(historyTime);

result = struct();
result.failed = false;
result.frcIndex = frcIndex;
if isfield(FRC,'labels')
    result.frcLabel = FRC.labels(frcIndex);
end
result.frequency = omega;
result.rotationSpeed = rotationSpeed;
result.rotationAcceleration = 0;
result.epsilon = epsilon;
result.forceAmplitude = epsilon*param.forcAmp;
result.forcingPhase = forcingPhase;
result.period = period;
result.timeStep = dt;
result.stepsPerPeriod = settings.stepsPerPeriod;
result.numberOfPeriods = settings.numberOfPeriods;
result.integrationRuntime = integrationRuntime;
result.outputDof = outdof;
result.initialConditionSource = initialConditionSource;
result.pssmInitialCondition = FRC.Z0_frc(:,frcIndex);
result.integrationInitialCondition = [qInitial;dqInitial];
result.dampingInfo = dampingInfo;
result.outputHistory = struct( ...
    't',historyTime, ...
    'cycle',cycleNumber, ...
    'q',solution.q_sys(outdof,:), ...
    'dq',solution.dq_sys(outdof,:), ...
    'ddq',solution.ddq_sys(outdof,:), ...
    'periodEndIndices',periodEndIndices);
result.finalRestartState = struct( ...
    'q',solution.q_sys(:,end), ...
    'dq',solution.dq_sys(:,end), ...
    'ddq',solution.ddq_sys(:,end), ...
    'a',solution.a(:,end));

if settings.saveFullSystemHistory
    result.fullSystemHistory = struct( ...
        't',historyTime, ...
        'q',solution.q_sys, ...
        'dq',solution.dq_sys, ...
        'ddq',solution.ddq_sys, ...
        'a',solution.a, ...
        'lambda',solution.lm);
end

if ~isfolder(resultsDirectory)
    mkdir(resultsDirectory);
end
frequencyTag = strrep(sprintf('%.6f',omega),'.','p');
fileStem = sprintf('RotBeam_full_history_%s_idx%03d_Omega%s', ...
    initialConditionSource,frcIndex,frequencyTag);
resultFile = fullfile(resultsDirectory,[fileStem,'.mat']);
save(resultFile,'result','settings','-v7.3');
fprintf('Complete output histories saved to:\n%s\n',resultFile);

%% Plot complete displacement, velocity and phase-history evolution
[timeFigure,phaseFigure] = plot_full_histories(result);
if settings.saveFigures
    exportgraphics(timeFigure, ...
        fullfile(resultsDirectory,[fileStem,'_time_history.png']), ...
        'Resolution',200);
    exportgraphics(phaseFigure, ...
        fullfile(resultsDirectory,[fileStem,'_phase_history.png']), ...
        'Resolution',200);
end

function [nodes,elements,constraints,rotationFrame] = build_beam_model( ...
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

function [qProjected,dqProjected] = project_initial_constraints(q,dq,elmnum)
qProjected = q;
dqProjected = dq;
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
[nodes,elements,constraints,~] = build_beam_model( ...
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
frequencies = frequencies(isfinite(frequencies) & frequencies>0);
assert(numel(frequencies)>=2, ...
    'At least two elastic modes are required for Rayleigh damping.');

omega1 = frequencies(1);
omega2 = frequencies(2);
alphaPerRatio = ...
    2*omega1*omega2*(omega2-omega1)/(omega2^2-omega1^2);
betaPerRatio = 2*(omega2-omega1)/(omega2^2-omega1^2);
Damp = param.dampRatio*(alphaPerRatio*M+betaPerRatio*K);
info = struct( ...
    'fixedDampingRatio',param.dampRatio, ...
    'alpha',param.dampRatio*alphaPerRatio, ...
    'beta',param.dampRatio*betaPerRatio, ...
    'referenceFrequencies',frequencies);
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

function [timeFigure,phaseFigure] = plot_full_histories(result)
t = result.outputHistory.t;
cycle = result.outputHistory.cycle;
q = result.outputHistory.q;
dq = result.outputHistory.dq;
outdofs = result.outputDof;
nOutputs = numel(outdofs);

timeFigure = figure('Name','Rotating-beam complete time histories', ...
    'Color','w');
timeLayout = tiledlayout(timeFigure,nOutputs,2, ...
    'TileSpacing','compact','Padding','compact');
title(timeLayout,sprintf(['Complete response histories, ' ...
    '\\Omega = %.6g rad/s'],result.frequency));
for k = 1:nOutputs
    nexttile(timeLayout)
    plot(t,q(k,:),'LineWidth',0.9)
    grid on
    ylabel(sprintf('q_{%d}',outdofs(k)))
    if k==nOutputs
        xlabel('Time [s]')
    end

    nexttile(timeLayout)
    plot(t,dq(k,:),'LineWidth',0.9)
    grid on
    ylabel(sprintf('dq_{%d}',outdofs(k)))
    if k==nOutputs
        xlabel('Time [s]')
    end
end

phaseFigure = figure('Name','Rotating-beam phase-history evolution', ...
    'Color','w');
phaseLayout = tiledlayout(phaseFigure,nOutputs,1, ...
    'TileSpacing','compact','Padding','compact');
title(phaseLayout,sprintf(['Phase histories colored by period number, ' ...
    '\\Omega = %.6g rad/s'],result.frequency));
for k = 1:nOutputs
    ax = nexttile(phaseLayout);
    surface(ax,[q(k,:);q(k,:)],[dq(k,:);dq(k,:)], ...
        zeros(2,numel(cycle)),[cycle;cycle], ...
        'FaceColor','none','EdgeColor','interp','LineWidth',1.0);
    hold(ax,'on')
    plot(ax,q(k,1),dq(k,1),'ko','MarkerFaceColor','w')
    plot(ax,q(k,end),dq(k,end),'kp','MarkerFaceColor','k')
    hold(ax,'off')
    view(ax,2)
    grid(ax,'on')
    xlabel(ax,sprintf('q_{%d}',outdofs(k)))
    ylabel(ax,sprintf('dq_{%d}',outdofs(k)))
    colorScale = colorbar(ax);
    colorScale.Label.String = 'Period number, t/T';
end
end
