q_sys = mySolver.Solution.q_sys;
coord_XY  = zeros(3,timenum);
deflexion = zeros(1,timenum);
for i=1:timenum
    e = q_sys(4:7,i);
    A = Base.genA(e);
%     theta=phi(3,i);
%     R=[cos(theta)  -sin(theta)  0
%        sin(theta)   cos(theta)  0
%        0            0           1];
    coord_XY(:,i) = A'*q_sys(end-6:end-4,i);
%     disQua=(coord_XY(1,i)-q_sys(end-5,i))^2+(coord_XY(2,i)-q_sys(end-4,i))^2;
% %     if 
% %         deflexion(i)=sqrt(disQua);
% %     end
%     deflexion(i)=-sqrt(disQua);
    
end