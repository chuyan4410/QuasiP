function test_frot_pdim_separation
%TEST_FROT_PDIM_SEPARATION Verify the dedicated PDIM Frot coefficients.

[myAssembly, z] = build_single_element_assembly();
axisR = [1; 2; 3];
axisR = axisR / norm(axisR);
mu0 = [1.2; -0.15];
delta = [0.3; 0.07];
mu = mu0 + delta;

% The time-domain path receives an equivalent prescribed quaternion
% state at zero frame angle. It remains independent of the PDIM routines.
er = [1; 0; 0; 0];
der = [0; 0.5 * mu(1) * axisR];
dder = [-0.25 * mu(1)^2; 0.5 * mu(2) * axisR];
myAssembly.rotationFrame = struct( ...
    'enabled', true, ...
    'e', repmat(er,1,2), ...
    'de', repmat(der,1,2), ...
    'dde', repmat(dder,1,2));
FrotFull = myAssembly.genSysFrot(z);

FrotExpanded = sparse(myAssembly.dof_sys,1);
for stateOrder = [1,3]
    FrotExpanded = FrotExpanded ...
        + myAssembly.genSysFrotPdimPart(z,axisR,mu0,stateOrder,[0,0]) ...
        + delta(1) * myAssembly.genSysFrotPdimPart( ...
            z,axisR,mu0,stateOrder,[1,0]) ...
        + delta(2) * myAssembly.genSysFrotPdimPart( ...
            z,axisR,mu0,stateOrder,[0,1]) ...
        + delta(1)^2 * myAssembly.genSysFrotPdimPart( ...
            z,axisR,mu0,stateOrder,[2,0]);
end
assert_relative_error(FrotExpanded,FrotFull,5e-12, ...
    'Delta-expanded coefficients do not match the complete Frot.');

% Each coefficient must remain homogeneous in the state variables.
stateScale = 1.37;
for stateOrder = [1,3]
    indices = {[0,0],[1,0],[0,1],[2,0]};
    for k = 1:numel(indices)
        part = myAssembly.genSysFrotPdimPart( ...
            z,axisR,mu0,stateOrder,indices{k});
        scaledPart = myAssembly.genSysFrotPdimPart( ...
            stateScale*z,axisR,mu0,stateOrder,indices{k});
        assert_relative_error(scaledPart,stateScale^stateOrder*part, ...
            2e-11,'A rotating-frame coefficient has the wrong state degree.');
    end
end

% The state polarization must recover the original homogeneous polynomial
% on diagonal inputs. Linear terms bypass the polarization helper.
for stateOrder = [1,3]
    inputx = repmat({z},1,stateOrder);
    multilinearValue = get_frot_pdim_multilinear_form( ...
        myAssembly,axisR,mu0,[1,0],inputx);
    polynomialValue = myAssembly.genSysFrotPdimPart( ...
        z,axisR,mu0,stateOrder,[1,0]);
    assert_relative_error(multilinearValue,polynomialValue,2e-11, ...
        'The multilinear diagonal value does not recover the polynomial.');
end

% Polarization must produce a symmetric multilinear map for distinct state
% directions. Adjacent transpositions generate every permutation, so
% checking each adjacent swap verifies the full permutation symmetry.
rng(731,'twister');
symmetryCases = { ...
    struct('stateOrder',3,'indices',{{[0,0],[1,0],[0,1],[2,0]}})};
for caseID = 1:numel(symmetryCases)
    stateOrder = symmetryCases{caseID}.stateOrder;
    inputx = cell(1,stateOrder);
    for inputID = 1:stateOrder
        inputx{inputID} = randn(size(z));
    end

    indices = symmetryCases{caseID}.indices;
    for indexID = 1:numel(indices)
        referenceValue = get_frot_pdim_multilinear_form( ...
            myAssembly,axisR,mu0,indices{indexID},inputx);
        for swapID = 1:stateOrder-1
            swappedInputs = inputx;
            swappedInputs([swapID,swapID+1]) = ...
                swappedInputs([swapID+1,swapID]);
            swappedValue = get_frot_pdim_multilinear_form( ...
                myAssembly,axisR,mu0,indices{indexID},swappedInputs);
            assert_relative_error(swappedValue,referenceValue,2e-10, ...
                'The polarized rotating-frame force is not symmetric.');
        end
    end
end

% There are no omega*dotOmega or dotOmega^2 structural terms.
for stateOrder = [1,3]
    crossPart = myAssembly.genSysFrotPdimPart( ...
        z,axisR,mu0,stateOrder,[1,1]);
    accelerationSquaredPart = myAssembly.genSysFrotPdimPart( ...
        z,axisR,mu0,stateOrder,[0,2]);
    assert(nnz(crossPart) == 0 && nnz(accelerationSquaredPart) == 0, ...
        'Unexpected mixed or dotOmega-squared parameter contribution.');
end

% Non-unit nodal quaternions and nonzero velocities exercise derivatives
% away from the unit-quaternion identities used by the old expression.
n = myAssembly.dof_sys;
myAssembly.genSysFrot(z);
elm = myAssembly.elements{1};
Krot = elm.Krot;
Crot = elm.Crot;
Kfd = zeros(n);
Cfd = zeros(n);
h = 1e-5;
for j = 1:n
    dz = zeros(2*n,1);
    dz(j) = h;
    Kfd(:,j) = (myAssembly.genSysFrot(z+dz) ...
        - myAssembly.genSysFrot(z-dz))/(2*h);
    dz = zeros(2*n,1);
    dz(n+j) = h;
    Cfd(:,j) = (myAssembly.genSysFrot(z+dz) ...
        - myAssembly.genSysFrot(z-dz))/(2*h);
end
kError = norm(Krot-Kfd,'fro')/max(1,norm(Krot,'fro'));
cError = norm(Crot-Cfd,'fro')/max(1,norm(Crot,'fro'));
assert(kError < 2e-8 && cError < 2e-8, ...
    'The rotating-frame analytic tangents disagree with central differences.');
skewError = norm(Crot+Crot.','fro')/max(1,norm(Crot,'fro'));
assert(skewError < 2e-11,'The rotating velocity tangent is not skew-symmetric.');
velocityForce = myAssembly.genSysFrot(z) ...
    - myAssembly.genSysFrot([z(1:n);zeros(n,1)]);
v = z(n+1:2*n);
powerError = abs(v.'*velocityForce)/max(1,norm(v)*norm(velocityForce));
assert(powerError < 2e-11,'The velocity-dependent rotating force does work.');

% The time-integration entry point must use the new analytic rotation
% tangents too. This fixture has q=q0, so both integration measures agree.
myAssembly.genSysFrot(z);
elm.get_K_P();
assert_relative_error(elm.Crot,Crot,2e-11, ...
    'The time integrator does not use the consistent rotation velocity tangent.');
assert_relative_error(elm.Krot,Krot,2e-11, ...
    'The time integrator does not use the consistent rotation position tangent.');

% Fifth-order rotation requests are unsupported, not silently zeroed.
rejected = false;
try
    myAssembly.genSysFrotPdimPart(z,axisR,mu0,5,[0,0]);
catch err
    rejected = strcmp(err.identifier,'RotFramePDIM:UnsupportedStateOrder');
end
assert(rejected,'A fifth-order rotation request was not rejected.');
fprintf('Krot FD: %.3e; Crot FD: %.3e; skew: %.3e; power: %.3e\n', ...
    kError,cError,skewError,powerError);
fprintf('test_frot_pdim_separation passed.\n');
end

function [myAssembly,z] = build_single_element_assembly()
material = Material.IsotropicLE(70e9,0.3,2700);
section = CrossSection.Rod3D(2.5e-3,1.8e-6,0.9e-6,0.9e-6,[0.85;0.85]);
timenum = 2;
zeroAM = repmat({zeros(3,1)},timenum,1);
zeroBodyForce = repmat({zeros(3,1)},timenum,1);
zeroDistributedForce = repmat({zeros(6,1)},timenum,1);

positions = [0, 0.5, 1.0];
quaternions = [ ...
     0.97,  0.93,  0.91; ...
     0.10, -0.07,  0.12; ...
    -0.06,  0.11,  0.04; ...
     0.08,  0.05, -0.09];
myNodes = cell(3,1);
q = zeros(21,1);
dq = zeros(21,1);
for nodeID = 1:3
    dofID = (7*nodeID-6:7*nodeID).';
    qNode = [positions(nodeID); 0.03*nodeID; -0.02*nodeID; ...
        quaternions(:,nodeID)];
    dqNode = (0.01:0.01:0.07).' * (0.5 + 0.2*nodeID);
    myNodes{nodeID} = Node.HM3DCurvRodNode_quat( ...
        dofID,material,section,zeros(3,1),qNode,zeros(7,1), ...
        zeros(7,1),timenum,zeroAM,zeroBodyForce, ...
        zeroDistributedForce,1);
    q(dofID) = qNode;
    dq(dofID) = dqNode;
end

myElement = Element.HM3DCurvRodL3_quat(myNodes,1);
constraints.Dconstraints = struct();
constraints.Fconstraints = struct();
myAssembly = Assembly.Assembly({myElement},myNodes,constraints, ...
    q,zeros(21,1),zeros(21,1));
z = [q;dq];
end

function assert_relative_error(actual,expected,tolerance,message)
scale = max([1,norm(actual),norm(expected)]);
assert(norm(actual-expected) <= tolerance*scale,message);
end
