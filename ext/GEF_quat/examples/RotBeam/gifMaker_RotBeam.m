%用于梁
close;
aviobj = VideoWriter('test.avi','Uncompressed AVI');%创建avi文件
% open(aviobj)%打开创建的avi文件
fig=figure;
q_sys = mySolver.Solution.q_sys;
position1=zeros(2*nodenum,timenum);
position2=zeros(2*nodenum,timenum);
coords=zeros(3*nodenum,timenum);
for i=1:timenum
    for j=1:elmnum
        position1(4*(j-1)+1,i)=q_sys(9*(j-1)+1,i);
        position1(4*(j-1)+2,i)=q_sys(9*(j-1)+2,i);
        position1(4*(j-1)+3,i)=q_sys(9*(j-1)+7,i);
        position1(4*(j-1)+4,i)=q_sys(9*(j-1)+8,i);
    end
    position1(end-1,i)=q_sys(dof_sys-5,i);
    position1(end,i)=q_sys(dof_sys-4,i);
    c=q_sys(4:7,i);
    A=Base.genA(c);
    for k=1:nodenum
        coords(3*k-2:3*k,i)=A*[elmlen*(k-1)/2;0;0];
        position2(2*(k-1)+1,i)=coords(3*k-2,i);
        position2(2*(k-1)+2,i)=coords(3*k-1,i);
    end
end

coord_XY=zeros(3,timenum);
for i=1:timenum
    c=q_sys(4:7,i);
    A=Base.genA(c);
    coord_XY(:,i)=A'*q_sys(end-6:end-4,i);
end

%% make gif
for i=1:timenum
    p = q_sys(:,i);
    if (mod(i,10) == 1)
        X1=position1(1:2:end,i);
        Y1=position1(2:2:end,i);
        X2=position2(1:2:end,i);
        Y2=position2(2:2:end,i);
        
        plot(X1,Y1,'bo-','LineWidth',2,'MarkerSize',3)
        hold on
        plot(X2,Y2,'r*-','LineWidth',2,'MarkerSize',3)
        hold off
        set(gcf,'Color',[1,1,1]); % modify background color
        set(gca,'xtick',x1,'fontname','Times New Roman','fontsize',12);
        set(gca,'ytick',y1,'fontname','Times New Roman','fontsize',12);
        xlabel('X(m)','fontname','Times New Roman','fontsize',12)
        ylabel('Y(m)','fontname','Times New Roman','fontsize',12)
        xlim([-len,len])
        ylim([-len,len])
        axis square;
        pause(0.01)
        %axis([-1 1 -1 1]);
        F = getframe(fig);
%         writeVideo(aviobj,F);
        im = frame2im(F);
        [I,map] = rgb2ind(im,256);
        %写入 GIF89a 格式文件   
        if i == 1
           imwrite(I,map,'test.gif','GIF', 'Loopcount',inf,'DelayTime',0.01);
        else
           imwrite(I,map,'test.gif','GIF','WriteMode','append','DelayTime',0.01);
        end   
    end
end