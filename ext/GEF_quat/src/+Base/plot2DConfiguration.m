%% data pre-procession
q_sys = mySolver.Solution.q_sys;

position1 = zeros(2*nodenum,timenum);
position2 = zeros(2*nodenum,timenum);
coords = zeros(3*nodenum,timenum);
for i = 1:timenum
    for j = 1:nodenum
        position1(2*(j-1)+1,i) = q_sys(7*(j-1)+1,i);
        position1(2*(j-1)+2,i) = q_sys(7*(j-1)+2,i);
    end
    e = q_sys(4:7,i);
    A = Base.genA(e);
    for k = 1:nodenum
        coords(3*k-2:3*k,i) = A*[elmlen*(k-1)/2;0;0];
        position2(2*(k-1)+1,i) = coords(3*k-2,i);
        position2(2*(k-1)+2,i) = coords(3*k-1,i);
    end
end
x1 = linspace(-len,len,5);
y1 = linspace(-len,len,5);

%% plot dynamical figure

fig = figure;
for i = 1:timenum
    X1 = position1(1:2:end,i);
    Y1 = position1(2:2:end,i);
    X2 = position2(1:2:end,i);
    Y2 = position2(2:2:end,i);
    
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
end

%% state at time: i
figure
i  = 101;%101:circle\
X1 = position1(1:2:end,i);
Y1 = position1(2:2:end,i);
plot(X1,Y1,'bo-','LineWidth',2,'MarkerSize',3)
set(gcf,'Color',[1,1,1]); % modify background color
set(gca,'xtick',x1,'fontname','Times New Roman','fontsize',12);
set(gca,'ytick',y1,'fontname','Times New Roman','fontsize',12);
xlabel('X(m)','fontname','Times New Roman','fontsize',12)
ylabel('Y(m)','fontname','Times New Roman','fontsize',12)
xlim([-len,len])
ylim([-len,len])
axis square;