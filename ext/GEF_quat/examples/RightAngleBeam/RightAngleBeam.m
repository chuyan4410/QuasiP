%% EXAMPLE: Right angle beam
clear
clc

tic
start = tic;
%% basical parameters
timestep = 0.01;
totalTime = 30;
t = 0:timestep:totalTime;
timenum = length(t);
%------------------------------material properties-------------------------
E    = 1e6;
nu   = 0;
G    = E/(1+nu)/2;
len  = 10;
rho  = 1;
Area = 1;
I2   = 1e-3;
I3   = 1e-3;
I1   = I2 + I3;
shearFactor = [2; 2];

%-----------------------------magnetic input-------------------------------
RM = [0; 0; 0];
AMList = cell(timenum,1);
gList  = cell(timenum,1);
fList  = cell(timenum,1);
KI = 1e4;
for i = 1:timenum
    AMList{i} = {zeros(3,1), zeros(3,3)};
    gList{i}  = {zeros(3,1)};
    fList{i}  = {zeros(6,1)};
end

myMaterial = Material.IsotropicLE(E, nu, rho);
mySection  = CrossSection.Rod3D(Area, I1, I2, I3, shearFactor);
%------------------------------discretization------------------------------
elmnum  = 5;
nodenum = 2 * elmnum + 1;
dof_sys = 7 * nodenum * 2;
elmlen  = len / elmnum;
myNodes = cell(2*nodenum,1);

dq0  = zeros(7,1);
ddq0 = zeros(7,1);
% segment 1
for n = 1:nodenum
    q0    = zeros(7,1);
    dofID = transpose(7 * n - 6: 7 * n);
    q0(2) = elmlen / 2 * (n - 1);
    q0(4) = cos(pi / 4);
    q0(7) = sin(pi / 4);
    myNodes{n} = Node.HM3DCurvRodNode_quat(dofID, myMaterial, mySection, ...
                      RM, q0, dq0, ddq0, timenum, AMList, gList, fList, KI);
end
% segment 2
for n = 1:nodenum
    q0    = zeros(7,1);
    dofID = transpose(7 * (n + nodenum) - 6: 7 * (n + nodenum));
    q0(1) = elmlen / 2 * (n - 1);
    q0(2) = len;
    q0(4) = 1;
    myNodes{n+nodenum} = Node.HM3DCurvRodNode_quat(dofID, myMaterial, mySection, ...
                              RM, q0, dq0, ddq0, timenum, AMList, gList, fList, KI);
end

myElements = cell(2*elmnum,1);
reduce = 1;
% segment 1
for e = 1:elmnum
    nodes = {myNodes{2*e-1}; myNodes{2*e}; myNodes{2*e+1}};
    myElements{e} = Element.HM3DCurvRodL3_quat(nodes, reduce);
end
% segment 2
for e = elmnum+1:2*elmnum
    nodes = {myNodes{2*e}; myNodes{2*e+1}; myNodes{2*e+2}};
    myElements{e} = Element.HM3DCurvRodL3_quat(nodes, reduce);
end

%% define the constraints
dof_constraint = 12;
constrainedDOF = {{1, [1, 2, 3, 5, 6, 7]}};
dispSet = [0; 0; 0; 0; 0; sin(pi/4)];
% here, the Cq only corresponding to constrainedDOF
CqSet = zeros(6, dof_sys);
CqSet(1:3,1:3) = eye(3);
CqSet(4:6,5:7) = eye(3);
KbSet = zeros(dof_sys, dof_sys);
dispList = cell(timenum,1);
for i = 1:timenum
    dispList{i} = {dispSet, CqSet, KbSet};
end
constrainedType = {{'displacement boundary'}; 
                   {'gradient discontinuity',{nodenum, nodenum+1,[cos(pi/4);0;0;-sin(pi/4)]}}};
% {nodenum, nodenum+1,[cos(pi/4);0;0;-sin(pi/4)]}
% ids of node_A, node_B, the rotation quaternion from node_A to node_B
myDispConstraint = Constraint.Constraint(constrainedDOF, ...
                   constrainedType, dispList, dof_constraint);

forcList = cell(timenum,1);
forcAmp = 50;
for i = 1:timenum
    if t(i) <= 1
        forcList{i} = {[0; 0; forcAmp*t(i); 0; 0; 0]};
    elseif t(i) > 1 && t(i) <= 2
        forcList{i} = {[0; 0; 2*forcAmp-forcAmp * t(i); 0; 0; 0]};
    else
        forcList{i} = {[0; 0; 0; 0; 0; 0]};
    end
end
dof_applied = {{nodenum, 3}};
myForcConstraint = Constraint.Constraint(dof_applied, ...
    {{'force boundary'}}, forcList, dof_constraint);

myConstraints.Dconstraints = myDispConstraint;
myConstraints.Fconstraints = myForcConstraint; 

%% solver
mySolver = Solver.G_alpha(myElements, myNodes, myConstraints, timestep, dof_sys);
set(mySolver, 'rho_s', 0)
mySolver.Integrate(t);

runtime = toc(start);
toc

%% plot
% 3D plot
myPlot3 = PostProcess.plotConfiguration(3, mySolver);
Xlim = [-1*len, 2*len];
Ylim = [-1*len, 2*len];
Zlim = [-len, len];
set(myPlot3,'Xlim', Xlim, 'Ylim', Ylim, 'Zlim', Zlim)
set(myPlot3, 'xTicks', [-10,-5,0,5,10,15,20], 'yTicks', [-10,-5,0,5,10,15,20], 'zTicks', [-10,-5,0,5,10])
set(myPlot3, 'saveGif', 1, 'gifFileName', 'rightAngle3D.gif')
myPlot3.plotBeamDyn();
myPlot3.plotBeamStatic(14.5);

% 2D plot
myPlot2 = PostProcess.plotConfiguration(2, mySolver);
Xlim = [-1*len, 2*len];
Ylim = [-1*len, 2*len];
set(myPlot2, 'Xlim', Xlim, 'Ylim', Ylim)
set(myPlot2, 'xTicks', [-10,-5,0,5,10,15,20], 'yTicks', [-10,-5,0,5,10,15,20])
set(myPlot2, 'pickDims', [1;2])
set(myPlot2, 'saveGif', 1, 'gifFileName', 'rightAngle2D.gif')
myPlot2.plotBeamDyn();
myPlot2.plotBeamStatic(14.5);