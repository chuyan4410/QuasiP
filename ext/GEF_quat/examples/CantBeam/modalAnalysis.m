%% EXAMPLE: cantilever beam 
clear
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

% parameter_set2
E   = 70e9;
nu  = 0.3;
G   = E/(1+nu)/2;
len = 1;
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

q0_sys   = [];
dq0_sys  = [];
ddq0_sys = [];
for n = 1:numel(myNodes)
    q0_sys   = [q0_sys;   myNodes{n}.q(:,1)];
    dq0_sys  = [dq0_sys;  myNodes{n}.dq(:,1)];
    ddq0_sys = [ddq0_sys; myNodes{n}.ddq(:,1)];
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

myAssembly = Assembly.Assembly(myElements, myNodes, myConstraints, ...
             q0_sys, dq0_sys, ddq0_sys);
myAssembly.genSysMatVec(2);

%% modal analysis
% second order form
M  = myAssembly.M;
M0 = myAssembly.genSysM0();
Kt = myAssembly.Kt;
Cq = myAssembly.Cq;
[U1,la] = eigs(Cq.' * Cq, dof_sys);
idx = diag(la) == 0;
U = U1(:,idx);
tmp = Cq*U;
M_hat  = U.' * M  * U;
Kt_hat = U.' * Kt * U;
E_max = 10;
[V2_hat,D2_hat] = eigs(Kt_hat, M_hat, E_max, 'smallestabs');
D2 = sqrt(diag(D2_hat));
V2 = U * V2_hat;

% first order form
dampRatio = 1e-3;
alpha = (2 * D2(1) * D2(2) * (D2(2) * dampRatio - D2(1) * dampRatio)) / (D2(2)^2-D2(1)^2);
beta = 2 * (D2(2) * dampRatio - D2(1) * dampRatio) / (D2(2)^2 - D2(1)^2);
Damp = alpha * M + beta * Kt;
A = [-Kt     sparse(dof_sys,dof_sys)     -Cq.'
     sparse(dof_sys,dof_sys)  speye(dof_sys)   sparse(dof_sys,size(Cq,1))
     Cq          sparse(size(Cq,1),dof_sys)  sparse(size(Cq,1),size(Cq,1))];
B = [Damp   M0   sparse(dof_sys,size(Cq,1))
     speye(dof_sys)   sparse(dof_sys,dof_sys)   sparse(dof_sys,size(Cq,1))
     sparse(size(Cq,1),dof_sys)  sparse(size(Cq,1),dof_sys) sparse(size(Cq,1),size(Cq,1))];
[V1_hat,D1_hat] = eigs(A, B, E_max, 'smallestabs');
D1 = diag(D1_hat);

figure
plot(imag(V1_hat(3:7:dof_sys,1)))

% rescaling of first order form
k1 = 1e6;
k2 = 1e-3;% t = k2 * \tau
dampRatio = 1e-3;
alpha = (2 * D2(1) * D2(2) * (D2(2) * dampRatio - D2(1) * dampRatio)) / (D2(2)^2-D2(1)^2);
beta = 2 * (D2(2) * dampRatio - D2(1) * dampRatio) / (D2(2)^2 - D2(1)^2);
Damp = alpha * M + beta * Kt;
A_r = [-Kt     sparse(dof_sys,dof_sys)     -Cq.'*k1
      sparse(dof_sys,dof_sys)  speye(dof_sys)   sparse(dof_sys,size(Cq,1))
      Cq*k1          sparse(size(Cq,1),dof_sys)  sparse(size(Cq,1),size(Cq,1))];
B_r = [Damp/k2   M0/(k2^2)   sparse(dof_sys,size(Cq,1))
     speye(dof_sys)   sparse(dof_sys,dof_sys)   sparse(dof_sys,size(Cq,1))
     sparse(size(Cq,1),dof_sys)  sparse(size(Cq,1),dof_sys) sparse(size(Cq,1),size(Cq,1))];
[V1_hat_r,D1_hat_r] = eigs(A_r, B_r, E_max, 'smallestabs');
D1_r = diag(D1_hat_r);

figure
plot(imag(V1_hat_r(3:7:dof_sys,1)))

%% plot modes
Y = zeros(nodenum, E_max);
X = q0_sys(1:7:end) / len;
% bending modes
Y(:,1) = V2(3:7:end,1);
Y(:,2) = V2(2:7:end,2);
Y(:,3) = V2(3:7:end,3);
Y(:,4) = V2(2:7:end,4);
Y(:,5) = V2(3:7:end,5);
Y(:,6) = V2(2:7:end,6);

figure
chooseMode = 3;
plot(X, Y(:,chooseMode))