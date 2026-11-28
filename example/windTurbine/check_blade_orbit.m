function report=check_blade_orbit(order,location,steps,varargin)
% Check an HB orbit by direct constrained Newmark integration for one period.
% From the example directory: check_blade_orbit(5,'peak',256).
if nargin<1, order=5; end
if nargin<2, location='peak'; end
if nargin<3, steps=256; end
m=build_blade(varargin{:}); n=m.n; nc=m.nc;
runName=sprintf('%s_H%d',m.runPrefix,order);
if endsWith(m.runPrefix,'_A60')
    [~,found]=coco_fname(runName,'bd.mat');
    if ~found, runName=strrep(runName,'_A60_H','_H'); end
end
bd=coco_bd_read(runName);
freq=coco_bd_col(bd,'k_OM1'); amps=blade_branch_amplitude(bd,m);
if strcmp(location,'peak'), [~,i]=max(amps);
elseif strcmp(location,'start'), i=1;
else, error('Select start or peak.'); end
C=coco_bd_col(bd([1,i+1],:),'C'); C=reshape(C,2*order+1,[]);
theta=(0:steps)*2*pi/steps;
Y=m.scale.*(C(1,:).'+C(2:2:end,:).'*cos((1:order).'*theta) ...
    +C(3:2:end,:).'*sin((1:order).'*theta));
q=m.q+Y(1:n,1); v=Y(n+(1:n),1); a=Y(2*n+(1:n),1); lm=m.lm+Y(3*n+1:end,1);
quat=reshape(1:n,7,[]); quat=quat(4:7,:);
E=reshape(m.q(quat(:))+Y(quat(:),:),4,[],steps+1);
unitError=max(abs(sum(E.^2,1)-1),[],'all');
tip=zeros(1,steps+1); tip(1)=blade_tip_response(m,Y(m.tip,1));
dt=2*pi/(steps*freq(i)); beta=dt^2/4; gamma=dt/2;
for k=1:steps
    qp=q+dt*v+beta*a; vp=v+gamma*a;
    for iteration=1:12
        q=qp+beta*a; v=vp+gamma*a;
        [force,J]=blade_force(m,q,v,a);
        e=reshape(q(quat(:)),4,[]);
        H=[m.R;sparse(repelem(1:nc-6,4),quat(:),2*e(:),nc-6,n)];
        constraint=[m.R*(q-m.reference);sum(e.^2,1).'-1];
        phase=cos(theta(k+1));
        force=force+m.D*v+H.'*lm-(m.F+m.loadTangent*(q-m.q))*phase;
        if norm(force,inf)<2e-4 && norm(constraint,inf)<1e-12, break; end
        kb=sparse(quat(:),quat(:),repelem(2*lm(7:end),4),n,n);
        tangent=beta*(J(:,1:n)+kb-m.loadTangent*phase)+gamma*(J(:,n+(1:n))+m.D)+J(:,2*n+(1:n));
        du=[tangent,full(H.');full(H),zeros(nc)]\(-[force;constraint/beta]);
        a=a+du(1:n); lm=lm+du(n+1:end);
    end
    assert(iteration<12,'Newmark failed at step %d.',k);
    tip(k+1)=blade_tip_response(m,q(m.tip)-m.q(m.tip));
end
reference=blade_tip_response(m,Y(m.tip,:)); amplitude=max(abs(reference));
report=struct('omega',freq(i),'steps',steps,'HBamplitude',amplitude, ...
    'timeAmplitude',max(abs(tip)),'waveErrorOverPeak',max(abs(tip-reference))/amplitude, ...
    ['tipClosure_' m.unit],abs(tip(end)-tip(1)),'HBunitConstraintError',unitError);
disp(report);
end
