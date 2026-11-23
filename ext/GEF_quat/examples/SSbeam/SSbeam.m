%% EXAMPLE: pure bending of cantilever beam
clear
clc

tic
start = tic;
%% basical parameters
timestep = 0.01;
totalTime = 20;
t = 0:timestep:totalTime;
timenum = length(t);
%------------------------------material properties-------------------------
E   = 1e4;% parameters from Zupan et al. Computational Mechanics, 2013
nu  = 0;
G   = 1e4;
len = 10;
rho = 1;
Area = 1;
I1 = 1e-1;
I2 = 1e-1;
I3 = 1e-1;
shearFactor = [1; 1];

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
elmnum  = 10;
nodenum = 2 * elmnum + 1;
dof_sys = 7 * nodenum;
elmlen  = len / elmnum;
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

%% define the constraints
dof_constraint = 5;
constrainedDOF = {{1,       [1, 2, 3]};
                  {nodenum, [2, 3]}};
dispSet = [0; 0; 0; 0; 0];
CqSet = zeros(dof_constraint, dof_sys);
CqSet(1:3,1:3) = eye(3);
CqSet(4:5,end-5:end-4) = eye(2);
KbSet = zeros(dof_sys, dof_sys);
dispList = cell(timenum,1);
for i = 1:timenum
    dispList{i} = {dispSet, CqSet, KbSet};
end
constrainedType = {{'displacement boundary'}; {[]}};
myDispConstraint = Constraint.Constraint(constrainedDOF, ...
                   constrainedType, dispList, dof_constraint);

forcList = cell(timenum,1);
forcAmp = -1;
for i = 1:timenum
    forcList{i} = {[0; forcAmp; 0; 0; 0; 0]};
end
dof_applied = {{elmnum + 1, 2}};
myForcConstraint = Constraint.Constraint(dof_applied, ...
    {{'force boundary'}}, forcList, dof_constraint);

myConstraints.Dconstraints = myDispConstraint;
myConstraints.Fconstraints = myForcConstraint; 

%% solver
mySolver = Solver.G_alpha(myElements, myNodes, myConstraints, timestep, dof_sys);
set(mySolver, 'dynamic', false, 'rho_s', 0)
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
set(myPlot3, 'saveGif', 1, 'gifFileName', 'cantBeam3D.gif')
myPlot3.plotBeamDyn();
myPlot3.plotBeamStatic(101 * timestep);

% 2D plot
myPlot2 = PostProcess.plotConfiguration(2, mySolver);
Xlim = [-1.1*len, 1.1*len];
Ylim = [-1.1*len, 1.1*len];
set(myPlot2, 'Xlim', Xlim, 'Ylim', Ylim)
set(myPlot2, 'xTicks', [-8,-4,0,4,8], 'yTicks', [-8,-4,0,4,8])
set(myPlot2, 'pickDims', [1;2])
set(myPlot2, 'rigidRotation', 1)
set(myPlot2, 'saveGif', 1, 'gifFileName', 'cantBeam2D.gif')
myPlot2.plotBeamDyn();
myPlot2.plotBeamStatic(101 * timestep);