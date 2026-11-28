function result=run_blade_frc(order,varargin)
% Run COCO startup and QuasiP install, then run this example in MATLAB.
% run_blade_frc(order,direction,amplitude,speed), e.g. run_blade_frc(5,'edge',20,9.68).
% speed is the rotor speed in rpm (default 12.1).
if nargin<1, order=5; end
validateattributes(order,{'numeric'},{'scalar','integer','positive'});
folder=fileparts(mfilename('fullpath')); previous=pwd;
cleanup=onCleanup(@()cd(previous)); cd(folder);
m=build_blade(varargin{:}); n=m.n; nc=m.nc; N=3*n+nc;
bounds=m.omega*[.95,1.05]; om=bounds(1);
% Scale Fourier coordinates for comparable continuation tangent components.
qScale=m.scale(1:n);
m.rows=[1./(om*qScale);1./(om^2*qScale);1e-5*ones(n,1);ones(nc,1)];
B=spdiags([ones(2*n,1);zeros(n+nc,1)],0,N,N);
C=sparse(N,N); C(1:n,n+(1:n))=-speye(n); C(n+(1:n),2*n+(1:n))=-speye(n);
C(2*n+(1:n),n+(1:n))=m.D;
left=spdiags(m.rows,0,N,N); right=spdiags(m.scale,0,N,N);
B=left*B*right; C=left*C*right;
% Linear harmonic response around the centrifugally deformed equilibrium.
dynamic=[m.K-om^2*m.M+1i*om*(m.G+m.D),m.H.';m.H,sparse(nc,nc)];
% Balance the mixed force and constraint units before solving the seed.
row=max(abs(dynamic),[],2); scaled=spdiags(1./row,0,n+nc,n+nc)*dynamic;
col=max(abs(scaled),[],1).';
sol=(scaled*spdiags(1./col,0,n+nc,n+nc)\([m.F;zeros(nc,1)]./row))./col;
q1=sol(1:n);
seed.k_wave=1; seed.Xc=([q1;1i*om*q1;-om^2*q1;sol(n+1:end)]./m.scale).';
samples=2^nextpow2(6*order+1); % fifth-degree force and retained harmonics
prob=coco_prob();
phi=struct('param','MP','method','mVCF','isvar',true,'sparse_jacobian',true, ...
    'type',{{'HB'}},'k',{{order}},'P',{{1}},'S',{{samples}},'S_wave',{{1024}});
prob=coco_set(prob,'QuasiP','ode_opts',...
    struct('n',N,'d',1,'e',1,'isfwithu',false), ...
    'amp_opts',struct('isamp',true,'outdof',m.tip,'Nsampling',1024), ...
    'sta_opts',struct('issta',false),'PHI_INITIAL',phi);
prob=coco_set(prob,'cont','PtMX',[0,500],'h0',.08,'h_min',1e-5,'h_max',.6, ...
    'NSV',20,'TOL',1e-8,'MaxRes',100);
prob=coco_set(prob,'corr','ItMX',25,'TOL',1e-8);
prob=ode_isol2QuasiP(prob,'',@(x,p,t)blade_rhs(x,t,m), ...
    @(x,p,t)jacobian(x,t,m),[],@(p)operators(B,C),[],seed,om);
runName=sprintf('%s_H%d',m.runPrefix,order);
bd=coco(prob,runName,[],1,{'k_OM1'},bounds);
freq=coco_bd_col(bd,'k_OM1'); amp=blade_branch_amplitude(bd,m);
result=table(freq(:),amp(:),'VariableNames',{'omega_rad_s',...
    ['tip_' m.direction '_' m.unit]});
assert(max(result.omega_rad_s)>=bounds(2)-2e-6,...
    'Continuation did not reach the upper frequency.');
writetable(result,fullfile(folder,[runName '.csv']));
fig=figure('Visible','off'); plot(freq,amp,'LineWidth',1.5);
xlabel('Excitation frequency [rad/s]','FontSize',13); 
ylabel(sprintf('Tip %s amplitude [%s]',m.direction,m.unit),'FontSize',13);
set(gca, 'FontSize', 12);
box off;
grid on; 
exportgraphics(fig,fullfile(folder,[runName '.png']),'Resolution',180);
set(fig,'Visible','on');
savefig(fig,fullfile(folder,[runName '.fig']));
close(fig);
end

function J=jacobian(x,t,m)
[~,J]=blade_rhs(x,t,m);
end

function [flags,depends,A,B,C,Ap,Bp,Cp]=operators(B,C)
flags=[0,1,1]; depends=[0,0,0]; A=[]; Ap=[]; Bp=[]; Cp=[];
end

