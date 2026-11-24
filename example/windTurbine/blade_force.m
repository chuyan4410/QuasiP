function [f,J] = blade_force(m,q,v,a)
% Reference-length GEF element force, batched over phase samples. Complex
% step differentiation is local to each element (63 input components).
n=m.n; ns=size(q,2); f=zeros(n,ns);
if nargout>1, J=zeros(n,3*n,ns); end
for k=1:size(m.gp,1)
    id=(k-1)*14+(1:21);
    X=[q(id,:);v(id,:);a(id,:)];
    if nargout>1
        % One baseline and 63 perturbed copies per sample.
        X=reshape(reshape(X,63,1,ns)+1i*1e-20*[zeros(63,1),eye(63)],63,[]);
    end
    fe=zeros(21,size(X,2));
    for g=1:size(m.gp,2)
        fe=fe+gauss_force(m.gp(k,g),m.speed,X);
    end
    if nargout>1
        block=reshape(fe,21,64,ns);
        cols=[id,id+n,id+2*n];
        J(id,cols,:)=J(id,cols,:)+imag(block(:,2:end,:))/1e-20;
        fe=reshape(real(block(:,1,:)),21,ns);
    end
    f(id,:)=f(id,:)+fe;
end
end

function f=gauss_force(g,omega,X)
persistent dG d2A
if isempty(dG), dG=Base.genElmdG_de(); d2A=Base.genElmd2A_de2(); end
q=X(1:21,:); v=X(22:42,:); a=X(43:63,:);
e=g.N*q; u=g.N*v; b=g.N*a; ep=g.Dn*q; r=g.Ds*q;
G=reshape(reshape(dG,12,4)*e,3,4,[]);
Gu=reshape(reshape(dG,12,4)*u,3,4,[]);
A=zeros(3,3,size(e,2)); Adot=A; dA=cell(4,1);
for j=1:4
    dA{j}=reshape(reshape(d2A(:,:,:,j),9,4)*e,3,3,[]);
    A=A+.5*dA{j}.*reshape(e(j,:),1,1,[]);
    Adot=Adot+dA{j}.*reshape(u(j,:),1,1,[]);
end
spin=[0;0;omega]; spinHat=Base.genSkew(spin);
axial=g.CA*(trans(A,r)-g.gamma)*g.c2;
bend=g.CI*(2*mult(G,ep)-g.kappa)*g.c2;
w=2*mult(G,u)+trans(A,spin); momentum=g.J*w;
dw=2*(mult(Gu,u)+mult(G,b))+trans(Adot,spin);
rot=2*(trans(Gu,momentum)+trans(G,g.J*dw));
for j=1:4
    angularGradient=2*dG(:,:,j)*u+trans(dA{j},spin);
    rot(j,:)=rot(j,:)-sum(angularGradient.*momentum,1) ...
        +sum(trans(dA{j},r).*axial,1)+2*sum((dG(:,:,j)*ep).*bend,1);
end
translation=g.mass*(g.S*a+2*spinHat*(g.S*v)+spinHat^2*(g.S*q));
f=g.weight*(g.S.'*translation+g.Ds.'*mult(A,axial) ...
    +g.N.'*rot+2*g.Dn.'*trans(G,bend));
end

function y=mult(A,x)
y=reshape(pagemtimes(A,reshape(x,size(x,1),1,[])),size(A,1),[]);
end

function y=trans(A,x)
y=reshape(pagemtimes(A,'transpose',reshape(x,size(x,1),1,[]),'none'),size(A,2),[]);
end
