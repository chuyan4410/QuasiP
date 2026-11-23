%% EXAMPLE: cantilever beam 
clear all
clc

tic
start = tic;
%% basical parameters
timestep = 0.01;
totalTime = 10;
t = 0:timestep:totalTime;
timenum = length(t);
%------------------------------material properties-------------------------
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
    q1_sys(dofID) = q0;
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
forcAmp = 10;
for i = 1:timenum
    forcList{i} = {[0; forcAmp; 0; 0; 0; 0]};
end
dof_applied = {{nodenum, 2}};
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
Ct = myAssembly.genSysCt(zti_sys);
Kt = myAssembly.genSysKt(zti_sys);
Fi = myAssembly.genSysFi(zti_sys);
Fe = myAssembly.genSysFe(zti_sys);
Fe1= myAssembly.genSysFe1(zti_sys);
Fe2= myAssembly.genSysFe2(zti_sys);
Fe3= myAssembly.genSysFe3(zti_sys);
Fe5= myAssembly.genSysFe5(zti_sys);
Fex= myAssembly.genSysFext(zti_sys);
C  = myAssembly.genSysC(zti_sys);
Cq = myAssembly.genSysCq(zti_sys);
Kb = myAssembly.genSysKb(zti_sys,ones(27,1));
dMdq = myAssembly.genSysdMdq(zti_sys);
Cqlmda2order = myAssembly.genSysCqlmda2order([zti_sys;ones(15,1)]);