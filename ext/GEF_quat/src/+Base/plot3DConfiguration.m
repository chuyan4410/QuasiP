%% data pre-procession (3D)

q_sys = mySolver.Solution.q_sys;
t = mySolver.Solution.t;  

nodeNUmber = 2 * nodenum;
Xlim = [-1*len, 2.5*len];
Ylim = [0, 2*len];
Zlim = [-len, len];

position1 = zeros(3*nodeNUmber, timenum); % current positions
position2 = zeros(3*nodeNUmber, timenum); % initial positions

coords = zeros(3*nodeNUmber, timenum);

for i = 1:timenum
    %% node coordinates
    for j = 1:nodeNUmber
        position1(3*(j-1)+1,i) = q_sys(7*(j-1)+1,i); % X
        position1(3*(j-1)+2,i) = q_sys(7*(j-1)+2,i); % Y
        position1(3*(j-1)+3,i) = q_sys(7*(j-1)+3,i); % Z

        position2(3*(j-1)+1,i) = q_sys(7*(j-1)+1,1);
        position2(3*(j-1)+2,i) = q_sys(7*(j-1)+2,1);
        position2(3*(j-1)+3,i) = q_sys(7*(j-1)+3,1);
    end
end


%% plot 3D dynamical figure

fig = figure;
ax = gca;
% timeText = text(ax, ...
%     0.02, 0.95, '', ...
%     'Units','normalized', ...
%     'FontSize',12, ...
%     'FontName','Times New Roman', ...
%     'FontWeight','bold');


for i = 1:timenum
    X1 = position1(1:3:end,i);
    Y1 = position1(2:3:end,i);
    Z1 = position1(3:3:end,i);

    X2 = position2(1:3:end,i);
    Y2 = position2(2:3:end,i);
    Z2 = position2(3:3:end,i);

    plot3(X1, Y1, Z1, 'bo-', ...
          'LineWidth', 2, 'MarkerSize', 4);
    hold on
    plot3(X2, Y2, Z2, 'r*-', ...
          'LineWidth', 2, 'MarkerSize', 4);
    hold off

    set(gcf,'Color',[1,1,1]);
    axis equal
    grid on
    xlabel('X (m)','FontName','Times New Roman')
    ylabel('Y (m)','FontName','Times New Roman')
    zlabel('Z (m)','FontName','Times New Roman')
    
    timeText = text(gca, ...
        0.02, 0.95, '', ...
        'Units','normalized', ...
        'FontSize',12, ...
        'FontName','Times New Roman');
    set(timeText, 'String', sprintf('t = %.4f s', t(i)));

    xlim(Xlim)
    ylim(Ylim)
    zlim(Zlim)


    view(3)        % key
    drawnow
    pause(0.01)
end

%% state at time: i (3D)

figure
i = 101;

X1 = position1(1:3:end,i);
Y1 = position1(2:3:end,i);
Z1 = position1(3:3:end,i);

plot3(X1, Y1, Z1, 'bo-', ...
      'LineWidth', 2, 'MarkerSize', 4);

set(gcf,'Color',[1,1,1]);
axis equal
grid on
xlabel('X (m)')
ylabel('Y (m)')
zlabel('Z (m)')
xlim(Xlim)
ylim(Ylim)
zlim(Zlim)
view(3)
