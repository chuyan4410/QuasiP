function [f,J]=blade_rhs(x,theta,m)
n=m.n; nc=m.nc; ns=size(x,2); X=m.scale.*x;
q=m.q+X(1:n,:); v=X(n+(1:n),:); a=X(2*n+(1:n),:);
lm=m.lm+X(3*n+1:end,:);
if nargout>1
    [force,K]=blade_force(m,q,v,a);
    J=zeros(3*n+nc,3*n+nc,ns);
else
    force=blade_force(m,q,v,a);
end
f=zeros(3*n+nc,ns); quat=reshape(1:n,7,[]); quat=quat(4:7,:);
for j=1:ns
    e=reshape(q(quat(:),j),4,[]);
    H=[m.R;sparse(repelem(1:nc-6,4),quat(:),2*e(:),nc-6,n)];
    f(2*n+(1:n),j)=(m.F+m.loadTangent*X(1:n,j))*cos(theta(j))-force(:,j)-H.'*lm(:,j);
    f(3*n+1:end,j)=-[m.R*(q(:,j)-m.reference);sum(e.^2,1).'-1];
    if nargout>1
        kb=sparse(quat(:),quat(:),repelem(2*lm(7:end,j),4),n,n);
        K(:,1:n,j)=K(:,1:n,j)+kb-m.loadTangent*cos(theta(j));
        J(2*n+(1:n),:,j)=-[K(:,:,j),full(H.')];
        J(3*n+1:end,1:n,j)=-H;
    end
end
f=m.rows.*f;
if nargout>1, J=J.*m.rows.*m.scale.'; end
end
