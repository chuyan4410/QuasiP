%% EXAMPLE: cant1
clear
clc

tic
start = tic;
%% basical parameters
load('cant1_ic.mat')
idx = 4;
Omega0 = FRC.om(idx);
T0 = 2 * pi / Omega0;
disp(['Omega0 = ',num2str(Omega0), 'rad/s'])
disp(['Omega0 = ',num2str(T0), 's'])

Nsteps = 100;
nT = 1000;
timestep = T0/Nsteps;
totalTime = timestep*Nsteps*nT;
t = 0:timestep:totalTime;
timenum = length(t);
%------------------------------material properties-------------------------
% parameter_set1
% E   = 6.895e10;% parameters from RenHui,JVA18,Spin-up Maneuver;Li ZhanWei[D]
% nu  = 0.3;
% G   = E/(1+nu)/2;
% len = 8;
% wid = 1.986e-2;
% hgt = 3.675e-2;
% rho = 2766.67;
% Area = wid * hgt;
% I2 = hgt * wid^3 / 12;
% I3 = wid * hgt^3 / 12;
% I1 = I2 + I3;
% shearFactor = [5/6; 5/6];

% parameter_set2
E   = 70e9;
nu  = 0.3;
G   = E/(1+nu)/2;
len = 4;
wid = 5e-2;
hgt = 4e-2;
rho = 2700;
Area = wid * hgt;
I2 = hgt * wid^3 / 12;
I3 = wid * hgt^3 / 12;
I1 = I2 + I3;
dampRatio = 1e-6;
shearFactor = [5/6; 5/6];

RM = [0;0;0];
AMList = cell(timenum,1);
gList  = cell(timenum,1);
fList  = cell(timenum,1);
KI = 1;
for i = 1:timenum
    AMList{i} = {zeros(3,1), zeros(3,3)};
    gList{i}  = {zeros(3,1)};
    fList{i}  = {zeros(6,1)};
end

myMaterial = Material.IsotropicLE(E, nu, rho);
mySection  = CrossSection.Rod3D(Area, I1, I2, I3, shearFactor);
%------------------------------discretization------------------------------
elmnum  = 4;
nodenum = 2 * elmnum + 1;
dof_sys = 7 * nodenum;
elmlen  = len / elmnum;
myNodes = cell(nodenum,1);

zInit = FRC.Z0_frc(1:2*dof_sys,idx);

q0   = zeros(7,1);
dq0  = zeros(7,1);
ddq0 = zeros(7,1);
q1_sys   = zeros(dof_sys,1);
dq1_sys  = zeros(dof_sys,1);
ddq1_sys = zeros(dof_sys,1);
for n = 1:nodenum
    dofID = transpose(7 * n - 6: 7 * n);
    q0(1) = elmlen / 2 * (n - 1);
    q0(4) = 1;
    myNodes{n} = Node.HM3DCurvRodNode_quat(dofID, myMaterial, mySection, ...
                      RM, q0, dq0, ddq0, timenum, AMList, gList, fList, KI);
    q1_sys(dofID)  = q0;
end

myElements = cell(elmnum,1);
reduce = 1;
for e = 1:elmnum
    nodes = {myNodes{2*e-1}; myNodes{2*e}; myNodes{2*e+1}};
    myElements{e} = Element.HM3DCurvRodL3_quat(nodes, reduce);
end

%% define the constraints
dof_constraint = 6;
constrainedDOF = {{1, [1, 2, 3, 5, 6, 7]}};
dispSet = [0; 0; 0; 0; 0; 0];
CqSet = zeros(dof_constraint, dof_sys);
CqSet(1:3,1:3) = eye(3);
CqSet(4:6,5:7) = eye(3);
KbSet = zeros(dof_sys, dof_sys);
dispList = cell(timenum,1);
for i = 1:timenum
    dispList{i} = {dispSet, CqSet, KbSet};
end
constrainedType = {{'displacement boundary'}; {[]}};
myDispConstraint = Constraint.Constraint(constrainedDOF, ...
                   constrainedType, dispList, dof_constraint);

forcList = cell(timenum,1);
forcAmp = 1e1;
Omega = Omega0;
for i = 1:timenum
    forcList{i} = {[0; forcAmp * cos(Omega*t(i)); 0; 0; 0; 0]};
end
dof_applied = {{nodenum, 2}};
myForcConstraint = Constraint.Constraint(dof_applied, ...
    {{'force boundary'}}, forcList, dof_constraint);

myConstraints.Dconstraints = myDispConstraint;
myConstraints.Fconstraints = myForcConstraint; 

%% solver
mySolver = Solver.G_alpha(myElements, myNodes, myConstraints, timestep, dof_sys);
Damp = gainRefDamping(nodenum, elmnum, myConstraints, dampRatio, ...
    elmlen, myMaterial, mySection, RM, timenum, AMList, gList, fList, KI);
set(mySolver, 'rho_s', 0, 'tol', 1e-5, 'Damp', Damp)
set(mySolver, 'Zic', zInit)
mySolver.Integrate(t);

runtime = toc(start);
toc

%% plot
% 3D plot
myPlot3 = PostProcess.plotConfiguration(3, mySolver);
Xlim = [-1.1*len, 1.1*len];
Ylim = [-1.1*len, 1.1*len];
Zlim = [-2, 2];
set(myPlot3, 'Xlim', Xlim, 'Ylim', Ylim, 'Zlim', Zlim)
set(myPlot3, 'xTicks', [-8,-4,0,4,8], 'yTicks', [-8,-4,0,4,8], 'zTicks', [-1,0,1])
set(myPlot3, 'rigidRotation', 1)
% myPlot3.plotBeamDyn();
myPlot3.plotBeamStatic(t(500));

% 2D plot
myPlot2 = PostProcess.plotConfiguration(2, mySolver);
Xlim = [-1.1*len, 1.1*len];
Ylim = [-1.1*len, 1.1*len];
set(myPlot2, 'Xlim', Xlim, 'Ylim', Ylim)
set(myPlot2, 'xTicks', [-8,-4,0,4,8], 'yTicks', [-8,-4,0,4,8])
set(myPlot2, 'rigidRotation', 1)
% myPlot2.plotBeamDyn();
myPlot2.plotBeamStatic(t(500));

% phase portrait
figure;
plot(mySolver.Solution.q_sys(end-5,end-Nsteps:end), ...
    mySolver.Solution.dq_sys(end-5,end-Nsteps:end))


%% generating damping matrix using reference configuration
function Damp = gainRefDamping(nodenum, elmnum, constraints, dampRatio, ...
    elmlen, myMaterial, mySection, RM, timenum, AMList, gList, fList, KI)
myNodes = cell(nodenum,1);

q0   = zeros(7,1);
dq0  = zeros(7,1);
ddq0 = zeros(7,1);
for n = 1:nodenum
    dofID = transpose(7 * n - 6: 7 * n);
    q0(1) = elmlen / 2 * (n - 1);
    q0(4) = 1;
    myNodes{n} = Node.HM3DCurvRodNode_quat(dofID, myMaterial, mySection, ...
                      RM, q0, dq0, ddq0, timenum, AMList, gList, fList, KI);
end

myElements = cell(elmnum,1);
reduce = 1;
for e = 1:elmnum
    nodes = {myNodes{2*e-1}; myNodes{2*e}; myNodes{2*e+1}};
    myElements{e} = Element.HM3DCurvRodL3_quat(nodes, reduce);
end
Damp = Base.genDamp(myElements, myNodes, constraints, dampRatio);
end