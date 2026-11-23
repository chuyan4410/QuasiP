%% NREL 5MW blade modal analysis from original OpenFAST BeamDyn data
% Data source:
% OpenFAST/r-test/glue-codes/openfast/5MW_Baseline/
%   NRELOffshrBsline5MW_BeamDyn.dat
%   NRELOffshrBsline5MW_BeamDyn_Blade.dat
%
% Axis mapping:
%   current x-axis  <- OpenFAST BeamDyn z-axis (blade span)
%   current y-axis  <- OpenFAST BeamDyn y-axis
%   current z-axis  <- OpenFAST BeamDyn x-axis
% Therefore, the 6-by-6 OpenFAST sectional matrices are mapped as
%   CA_mB = K([3 2 1], [3 2 1])
%   CI_mB = K([6 5 4], [6 5 4])
%   Jrho  = M([6 5 4], [6 5 4])
% and the blade root is placed at the origin.
clear
clc

%% Paths
thisFile = mfilename('fullpath');
exampleDir = fileparts(thisFile);
repoRoot = fileparts(fileparts(exampleDir));
addpath(genpath(fullfile(repoRoot, 'src')));
addpath(genpath(fullfile(repoRoot, 'ext')));

openfastDir = fullfile(exampleDir, 'OpenFAST_5MW_Baseline');
[beamDynFile, beamDynBladeFile] = ensureOpenFASTBeamDynFiles(openfastDir);

%% Analysis options
numModes = 12;
enableRotFrame = true;
rotSpeedRPM = 0;  % Set to 12.1 for rated-speed rotating-frame linearization.
azimuth0 = 0;
OmegaHub = rotSpeedRPM * 2 * pi / 60;
gravity = [0; 0; 0];
reduce = 1;

%% Read original OpenFAST BeamDyn key points and sectional matrices
[kp_zr, initialTwist] = readBeamDynGeometry(beamDynFile);
[eta, K_openfast, M_openfast] = readBeamDynBladeSectionData(beamDynBladeFile);

if numel(kp_zr) ~= numel(eta)
    error('OpenFAST input mismatch: BeamDyn geometry and blade section counts differ.');
end

radius = kp_zr(:);
twist = initialTwist(:);
nodenum = numel(radius);
elmnum = (nodenum - 1) / 2;
dof_sys = 7 * nodenum;

if abs(elmnum - round(elmnum)) > eps
    error('The current L3 element mesh expects an odd number of OpenFAST stations.');
end

%% Load containers required by the node class
timenum = 1;
AMList = {{zeros(3,1), zeros(3,3)}};
gList = {{gravity}};
fList = {{zeros(6,1)}};
RM = [0; 0; 0];
dummyMaterial = Material.IsotropicLE(0, 0, 0);

%% Build OOP quaternion beam model from OpenFAST sectional resultants
myNodes = cell(nodenum, 1);
q0 = zeros(7,1);
dq0 = zeros(7,1);
ddq0 = zeros(7,1);

for n = 1:nodenum
    [rhoA, Jrho, CA_mB, CI_mB] = mapOpenFASTSectionToCurrentAxes( ...
        K_openfast(:,:,n), M_openfast(:,:,n));
    mySection = CrossSection.DirectRod3D(rhoA, Jrho, CA_mB, CI_mB);

    dofID = transpose(7 * n - 6 : 7 * n);
    q0(:) = 0;
    q0(1) = radius(n);
    q0(4:7) = quatFromAxisAngle([1; 0; 0], deg2rad(twist(n)));

    myNodes{n} = Node.HM3DCurvRodNode_quat(dofID, dummyMaterial, mySection, ...
        RM, q0, dq0, ddq0, timenum, AMList, gList, fList, 1);
end

myElements = cell(elmnum, 1);
for e = 1:elmnum
    nodes = {myNodes{2*e-1}; myNodes{2*e}; myNodes{2*e+1}};
    myElements{e} = Element.HM3DCurvRodL3_quat(nodes, reduce);
end

rotFrame = buildOperatingPointRotationFrame(enableRotFrame, azimuth0, OmegaHub);
for e = 1:elmnum
    myElements{e}.setRotationFrame(rotFrame);
end

%% Cantilever root constraint plus quaternion unit constraints
[q0_sys, dq0_sys, ddq0_sys] = Base.sysInitCond(myNodes);
dof_constraint = 6;
constrainedDOF = {{1, [1, 2, 3, 5, 6, 7]}};
dispSet = [0; 0; 0; 0; 0; 0];
CqSet = zeros(dof_constraint, dof_sys);
CqSet(1:3, 1:3) = eye(3);
CqSet(4:6, 5:7) = eye(3);
KbSet = zeros(dof_sys, dof_sys);
dispList = {{dispSet, CqSet, KbSet}};
constrainedType = {{'displacement boundary'}; {[]}};
myDispConstraint = Constraint.Constraint(constrainedDOF, ...
    constrainedType, dispList, dof_constraint);

forcList = {{zeros(6,1)}};
dof_applied = {{nodenum, 2}};
myForcConstraint = Constraint.Constraint(dof_applied, ...
    {{'force boundary'}}, forcList, dof_constraint);

myConstraints.Dconstraints = myDispConstraint;
myConstraints.Fconstraints = myForcConstraint;

%% Assemble and solve the constrained generalized eigenvalue problem
myAssembly = Assembly.Assembly(myElements, myNodes, myConstraints, ...
    q0_sys, dq0_sys, ddq0_sys);
myAssembly.genSysMatVec(1);

M = full(myAssembly.M);
Kt = full(myAssembly.Kt);
Cq = full(myAssembly.Cq);

U = null(Cq, 'r');
M_hat = U.' * M * U;
Kt_hat = U.' * Kt * U;

[V_hat, D_hat] = eig(Kt_hat, M_hat);
lambda = real(diag(D_hat));
valid = isfinite(lambda) & lambda > 0;
[omega, sortID] = sort(sqrt(lambda(valid)));
V_hat = V_hat(:, valid);
V_hat = V_hat(:, sortID);

freqHz = omega / (2 * pi);
numModes = min(numModes, numel(freqHz));
freqHz = freqHz(1:numModes);
omega = omega(1:numModes);
modeShape = U * V_hat(:, 1:numModes);

modalResults = table((1:numModes).', omega, freqHz, ...
    'VariableNames', {'Mode', 'CircularFrequency_rad_s', 'Frequency_Hz'});

disp('NREL 5MW blade modal frequencies from OpenFAST 49-station BeamDyn data:')
disp(modalResults)

%% Plot modal displacement components in current axes
figure('Name', 'NREL 5MW OpenFAST BeamDyn 49-station blade modes');
tiledlayout(3, 1)

nexttile
plot(eta, modeShape(3:7:end, [1;3;5]), 'LineWidth', 1.1)
xlabel('r/R')
ylabel('Out-of-plane component')
grid on
title('Flapwise-like modal components')

nexttile
plot(eta, modeShape(2:7:end, [2;4]), 'LineWidth', 1.1)
xlabel('r/R')
ylabel('In-plane component')
grid on
title('Edgewise-like modal components')

nexttile
plot(eta, modeShape(5:7:end, 6), 'LineWidth', 1.1)
xlabel('r/R')
ylabel('torsional component')
grid on
title('Torsional-like modal components')

%% Local helpers
function [beamDynFile, beamDynBladeFile] = ensureOpenFASTBeamDynFiles(openfastDir)
    if ~exist(openfastDir, 'dir')
        mkdir(openfastDir);
    end

    baseUrl = ['https://raw.githubusercontent.com/OpenFAST/r-test/', ...
        'refs/heads/main/glue-codes/openfast/5MW_Baseline/'];
    beamDynFile = fullfile(openfastDir, 'NRELOffshrBsline5MW_BeamDyn.dat');
    beamDynBladeFile = fullfile(openfastDir, 'NRELOffshrBsline5MW_BeamDyn_Blade.dat');

    downloadIfMissing(beamDynFile, [baseUrl, 'NRELOffshrBsline5MW_BeamDyn.dat']);
    downloadIfMissing(beamDynBladeFile, [baseUrl, 'NRELOffshrBsline5MW_BeamDyn_Blade.dat']);
end

function downloadIfMissing(fileName, url)
    if exist(fileName, 'file')
        return
    end

    try
        websave(fileName, url);
    catch ME
        error('OpenFAST:DownloadFailed', ...
            ['Could not download %s.\nPlace the OpenFAST file at:\n%s\n', ...
             'Original error: %s'], url, fileName, ME.message);
    end
end

function [kp_zr, initialTwist] = readBeamDynGeometry(fileName)
    lines = string(readlines(fileName));
    header = find(contains(lines, 'kp_xr'), 1, 'first');
    if isempty(header)
        error('OpenFAST:InvalidBeamDynFile', 'Could not find BeamDyn key-point table.');
    end

    n = sscanf(lines(header - 2), '%d');
    if isempty(n)
        n = 49;
    end

    dataStart = header + 2;
    kp_zr = zeros(n, 1);
    initialTwist = zeros(n, 1);
    for i = 1:n
        vals = sscanf(lines(dataStart + i - 1), '%f');
        kp_zr(i) = vals(3);
        initialTwist(i) = vals(4);
    end
end

function [eta, K, M] = readBeamDynBladeSectionData(fileName)
    lines = string(readlines(fileName));
    stationLine = find(contains(lines, 'station_total'), 1, 'first');
    if isempty(stationLine)
        error('OpenFAST:InvalidBladeFile', 'Could not find station_total.');
    end
    stationTotal = sscanf(lines(stationLine), '%d');

    dataStart = find(strtrim(lines) == '0.000000', 1, 'first');
    if isempty(dataStart)
        error('OpenFAST:InvalidBladeFile', 'Could not find first BeamDyn station.');
    end

    eta = zeros(stationTotal, 1);
    K = zeros(6, 6, stationTotal);
    M = zeros(6, 6, stationTotal);
    iLine = dataStart;
    for station = 1:stationTotal
        while strlength(strtrim(lines(iLine))) == 0
            iLine = iLine + 1;
        end
        eta(station) = str2double(strtrim(lines(iLine)));
        iLine = iLine + 1;

        while strlength(strtrim(lines(iLine))) == 0
            iLine = iLine + 1;
        end
        for row = 1:6
            K(row,:,station) = sscanf(lines(iLine), '%f').';
            iLine = iLine + 1;
        end

        while strlength(strtrim(lines(iLine))) == 0
            iLine = iLine + 1;
        end
        for row = 1:6
            M(row,:,station) = sscanf(lines(iLine), '%f').';
            iLine = iLine + 1;
        end
    end
end

function [rhoA, Jrho, CA_mB, CI_mB] = mapOpenFASTSectionToCurrentAxes(K, M)
    translationalMap = [3, 2, 1];
    rotationalMap = [6, 5, 4];

    rhoA = mean(diag(M(translationalMap, translationalMap)));
    CA_mB = K(translationalMap, translationalMap);
    CI_mB = K(rotationalMap, rotationalMap);
    Jrho = M(rotationalMap, rotationalMap);
end

function rotFrame = buildOperatingPointRotationFrame(enabled, azimuth, omega)
    if ~enabled
        rotFrame = struct('enabled', false);
        return
    end

    er = [cos(azimuth / 2); 0; 0; sin(azimuth / 2)];
    der = 0.5 * omega * [-sin(azimuth / 2); 0; 0; cos(azimuth / 2)];
    dder = -0.25 * omega^2 * [cos(azimuth / 2); 0; 0; sin(azimuth / 2)];

    rotFrame = struct('enabled', true, ...
                      'e', er, ...
                      'de', der, ...
                      'dde', dder, ...
                      'Lambda', Base.genA(er), ...
                      'Omega', 2 * Base.genG(er) * der, ...
                      'dotOmega', 2 * Base.genG(er) * dder);
end

function q = quatFromAxisAngle(axis, angle)
    axis = axis(:) / norm(axis);
    q = [cos(angle / 2); axis * sin(angle / 2)];
end
