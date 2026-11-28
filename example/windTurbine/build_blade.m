function m = build_blade(direction,amplitude,rpm)
% Build the NREL 5 MW structural model at its rotating equilibrium.
arguments
    direction (1,1) string {mustBeMember(direction,["flap","edge","torsion"])} = "flap"
    amplitude (1,1) double {mustBePositive,mustBeFinite} = 60
    rpm (1,1) double {mustBeReal,mustBeFinite} = 12.1
end
blade=nrel5mwBladeSectionData();
nn=numel(blade.radius); n=7*nn; nc=nn+6;
speed=rpm*2*pi/60; cone=2.5; zeta=.0048;
q=zeros(n,1); nodes=cell(nn,1); elements=cell((nn-1)/2,1);
material=Material.IsotropicLE(0,0,0);
for j=1:nn
    id=7*j-6:7*j; twist=blade.twist(j);
    q(id(1:3))=blade.radius(j)*[cosd(cone);0;sind(cone)];
    q(id(4:7))=Base.quatMultiply([cosd(cone/2);0;-sind(cone/2);0], ...
        [cosd(twist/2);sind(twist/2);0;0]);
    inertia=[blade.flpIner(j)+blade.edgIner(j),blade.flpIner(j),blade.edgIner(j)];
    section=CrossSection.DirectRod3D(blade.mass(j),diag(inertia), ...
        diag([blade.EA(j),blade.GA(j),blade.GA(j)]), ...
        diag([blade.GJ(j),blade.flpStff(j),blade.edgStff(j)]));
    nodes{j}=Node.HM3DCurvRodNode_quat(id.',material,section,zeros(3,1), ...
        q(id),zeros(7,1),zeros(7,1),2,repmat({{zeros(3,1),zeros(3)}},2,1), ...
        repmat({{zeros(3,1)}},2,1),repmat({{zeros(6,1)}},2,1),1);
end
F=zeros(n,1); loadTangent=sparse(n,n);
switch direction
    case "flap"
        component=3; mode=1;
    case "edge"
        component=2; mode=2;
    case "torsion"
        component=1; mode=6;
end
axisLoad=zeros(3,1); axisLoad(component)=amplitude;
% For a material-axis couple, virtual work is 2*(G(e)*delta_e)'*moment.
% G is linear in e, so the generalized follower load has a constant tangent.
dG=Base.genElmdG_de(); torqueMap=zeros(4);
for j=1:4, torqueMap(:,j)=2*dG(:,:,j).'*axisLoad; end
for k=1:numel(elements)
    el=Element.HM3DCurvRodL3_quat(nodes(2*k-1:2*k+1),1);
    elements{k}=el; [w,xi]=Base.gaussQuadrature(el.numGauss);
    for g=1:numel(w)
        [S,N,Ds,Dn]=Base.genShapeF(xi(g)); el.updateSectionProperties(xi(g));
        ds=norm(Ds*el.q0_elm); weight=w(g)*ds;
        gp=struct('S',S,'N',N,'Ds',Ds/ds,'Dn',Dn/ds,'weight',weight, ...
            'mass',el.rhoA,'J',el.Jrho,'CA',el.CA_mB,'CI',el.CI_mB);
        e=N*el.q0_elm; A=Base.genA(e); r=gp.Ds*el.q0_elm;
        gp.gamma=A.'*r; gp.kappa=2*Base.genG(e)*gp.Dn*el.q0_elm;
        gp.c2=Base.genc(A,r)^2; m.gp(k,g)=gp;
        id=el.dofID;
        if direction=="torsion"
            loadTangent(id,id)=loadTangent(id,id)+N.'*torqueMap*N*weight;
        else
            F(id)=F(id)+S.'*axisLoad*weight;
        end
    end
end
if direction~="torsion"
    assert(abs(sum(F(component:7:end))/amplitude-diff(blade.radius([1 end])))<1e-8);
end
root=[1:3,5:7]; I=speye(n); R=I(root,:);
constraints.Dconstraints=Constraint.Constraint({{1,root}},{{'displacement boundary'};{[]}}, ...
    repmat({{q(root),R,sparse(n,n)}},2,1),6);
constraints.Fconstraints=Constraint.Constraint({{nn,3}},{{'force boundary'}}, ...
    repmat({{zeros(6,1)}},2,1),6);
asm=Assembly.Assembly(elements,nodes,constraints,q,zeros(n,1),zeros(n,1));
asm.rotationFrame=struct('enabled',true,'e',repmat([1;0;0;0],1,2), ...
    'de',repmat([0;0;0;speed/2],1,2),'dde',repmat([-speed^2/4;0;0;0],1,2));
m.reference=q; lm=zeros(nc,1); zero=zeros(n,1);
for it=1:30
    asm.genSysForceTangents([q;zero]);
    H=asm.genSysCq([q;zero]); c=asm.genSysC([q;zero]);
    K=asm.Kt+asm.genSysKb([q;zero],lm);
    residual=[asm.Fe+asm.Frot+H.'*lm;c];
    if norm(residual(1:n),inf)<1e-4 && norm(c,inf)<1e-12, break; end
    delta=[K,H.';H,sparse(nc,nc)]\(-residual);
    q=q+delta(1:n); lm=lm+delta(n+1:end);
end
assert(norm(residual(1:n),inf)<1e-4 && norm(c,inf)<1e-12,'Equilibrium did not converge.');
Z=null(full(H)); M=asm.M; G=asm.Ct; Mr=Z.'*M*Z; Kr=Z.'*K*Z;
% Classical modal damping in the equilibrium tangent space. Check the
% actual pole ratios with the rotating-frame gyroscopic matrix retained.
[V,d]=eig((Kr+Kr.')/2,(Mr+Mr.')/2,'vector'); assert(all(d>0));
V=V./sqrt(diag(V.'*Mr*V)).'; weights=M*Z*V;
D=weights*diag(2*zeta*sqrt(d))*weights.';
Gr=Z.'*G*Z; Dr=Z.'*D*Z;
poles=polyeig(Kr,Gr+Dr,Mr); poles=poles(imag(poles)>0);
m.ratios=-real(poles)./abs(poles);
[shapes,p]=polyeig(Kr,Gr,Mr); pos=find(imag(p)>0);
[~,ord]=sort(imag(p(pos))); pos=pos(ord);
if direction=="torsion"
    tip=n-3:n; m.unit='deg';
else
    tip=n-7+component; m.unit='m';
    shape=Z*shapes(:,pos(mode));
    assert(abs(shape(tip))>abs(shape(n-7+5-component)), ...
        'Mode %d is not %s dominated at the tip.',mode,direction);
end
m.direction=char(direction); m.mode=mode;
m.runPrefix=sprintf('%s_%srpm_A%s',direction,strrep(num2str(rpm,15),'.','p'), ...
    strrep(num2str(amplitude,17),'.','p'));
m.omega=imag(p(pos(mode))); m.tip=tip; m.q=q; m.lm=lm; m.R=R;
m.speed=speed; m.loadTangent=loadTangent; m.F=F+loadTangent*q;
m.K=K; m.M=M; m.G=G; m.D=sparse(D); m.H=H;
qScale=repmat([5;5;5;.1;.1;.1;.1],nn,1); seedOmega=.95*m.omega;
m.scale=[qScale;seedOmega*qScale;seedOmega^2*qScale;1e6*ones(nc,1)];
m.n=n; m.nc=nc; m.equilibriumResidual=norm(residual(1:n),inf);
m.elasticFrequencies=sort(imag(p(pos)));
fprintf('Equilibrium %.3g; %s mode %d: %.9g rad/s; damping range %.6g to %.6g\n', ...
    m.equilibriumResidual,direction,mode,m.omega,min(m.ratios),max(m.ratios));
% Independently check the batched force against the native GEF assembly.
rngState=rng; rng(71); restore=onCleanup(@()rng(rngState));
qt=q+.002*randn(n,1); v=.04*randn(n,1); a=.1*randn(n,1);
asm.genSysForceTangents([qt;v]);
native=asm.M*a+asm.Fe-asm.Fi+asm.Frot;
batched=blade_force(m,qt,v,a);
m.forceCheck=norm(native-batched)/max(1,norm(native));
assert(m.forceCheck<1e-10,'Batched and native GEF force differ by %.3g',m.forceCheck);
end

