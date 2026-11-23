%% EXAMPLE: Right angle beam
% clear
% clc

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

q1_sys   = [];
dq1_sys  = [];
ddq1_sys = [];
for n = 1:numel(myNodes)
    q1_sys   = [q1_sys;   myNodes{n}.q(:,1)];
    dq1_sys  = [dq1_sys;  myNodes{n}.dq(:,1)];
    ddq1_sys = [ddq1_sys; myNodes{n}.ddq(:,1)];
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

%% Assembly
myAssembly = Assembly.Assembly(myElements, myNodes, myConstraints, ...
             q1_sys, dq1_sys, ddq1_sys);
myAssembly.genSysMatVec(2);
zti_sys = [q1_sys; dq1_sys];
M0 = myAssembly.genSysM0();
Mnl= myAssembly.genSysMnl(zti_sys);
Fi = myAssembly.genSysFi(zti_sys);
Fe1= myAssembly.genSysFe1(zti_sys);
Fe2= myAssembly.genSysFe2(zti_sys);
Fe3= myAssembly.genSysFe3(zti_sys);
Fe5= myAssembly.genSysFe5(zti_sys);
Fex= myAssembly.genSysFext(zti_sys);
C  = myAssembly.genSysC(zti_sys);
Cq = myAssembly.genSysCq(zti_sys);

