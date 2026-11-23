%% EXAMPLE: rotationary beam
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
E   = 6.895e10;% parameters from RenHui,JVA18,Spin-up Maneuver;Li ZhanWei[D]
nu  = 0.3;
G   = E/(1+nu)/2;
len = 8;
wid = 1.986e-2;
hgt = 3.675e-2;
rho = 2766.67;
Area = wid * hgt;
I2 = hgt * wid^3 / 12;
I3 = wid * hgt^3 / 12;
I1 = I2 + I3;
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
dof_constraint = 6;
constrainedDOF = {{1, [1, 2, 3, 5, 6, 7]}};
CqSet = zeros(dof_constraint, dof_sys);
CqSet(1:3,1:3) = eye(3);
CqSet(4:6,5:7) = eye(3);
KbSet = zeros(dof_sys, dof_sys);

omegas = 4;
Ts = 15;
phiSet = zeros(timenum,1);

dispList = cell(timenum,1);
for i = 1:timenum
    if t(i) < Ts
        phiSet(i) = omegas / Ts * (t(i)^2 / 2 + ( Ts / 2 / pi)^ 2 * (cos(2 * pi * t(i) / Ts) - 1));
    else
        phiSet(i) = omegas * (t(i) - Ts/2);
    end
    dispSet = [0; 0; 0; 0; 0; sin(phiSet(i) / 2)];
    dispList{i} = {dispSet, CqSet, KbSet};
end
constrainedType = {{'displacement boundary'}; {[]}};
myDispConstraint = Constraint.Constraint(constrainedDOF, ...
                   constrainedType, dispList, dof_constraint);

forcList = cell(timenum,1);
forcAmp = 0;
for i = 1:timenum
    forcList{i} = {[0; forcAmp * sin(t(i)); 0; 0; 0; 0]};
end
dof_applied = {{nodenum, 2}};
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
Xlim = [-1.1*len, 1.1*len];
Ylim = [-1.1*len, 1.1*len];
Zlim = [-2, 2];
set(myPlot3, 'Xlim', Xlim, 'Ylim', Ylim, 'Zlim', Zlim)
set(myPlot3, 'xTicks', [-8,-4,0,4,8], 'yTicks', [-8,-4,0,4,8], 'zTicks', [-1,0,1])
set(myPlot3, 'rigidRotation', 1)
set(myPlot3, 'saveGif', 1, 'gifFileName', 'beam3D_animation.gif')
% myPlot3.plotBeamDyn();
myPlot3.plotBeamStatic(14.5);

% 2D plot
myPlot2 = PostProcess.plotConfiguration(2, mySolver);
Xlim = [-1.1*len, 1.1*len];
Ylim = [-1.1*len, 1.1*len];
set(myPlot2, 'Xlim', Xlim, 'Ylim', Ylim)
set(myPlot2, 'xTicks', [-8,-4,0,4,8], 'yTicks', [-8,-4,0,4,8])
set(myPlot2, 'rigidRotation', 1)
set(myPlot2, 'saveGif', 1, 'gifFileName', 'beam2D_animation.gif')
% myPlot2.plotBeamDyn();
myPlot2.plotBeamStatic(14.5);