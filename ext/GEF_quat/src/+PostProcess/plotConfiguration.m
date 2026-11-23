classdef plotConfiguration < matlab.mixin.SetGetExactNames
    %PLOTCONFIGURATION : plot the beam configuration
    
    properties
        Dim                 % 2D or 3D plot
        q_sys               % systematical generalized coordinates
        t                   % time
        timenum             % the number of time steps
        nodeNumber          % node number of beam
        XYZ                 % the position of last node
        Xlim = []           % The range of values on the X-axis
        Ylim = []           % The range of values on the Y-axis
        Zlim = []           % The range of values on the Z-axis
        positiont           % current positions
        position0           % initial positions
        rigidRotation = 0   % plot the beam's rigid rotation
        pickDims = [1;2;3]  % picked 1: X-direction; 2: Y-direction; 3: Z-direction;
        xTicks = []
        yTicks = []
        zTicks = []
        xyzLabel = {'X(m)', 'Y(m)', 'Z(m)'};

        saveGif = 0         % Switch to enable/disable GIF saving (0: off, 1: on)
        gifFileName  = 'beam_animation.gif' % Name of the saved GIF file
        gifDelayTime = 0.02;% Delay time between frames in seconds (smaller = faster)
    end
    
    methods
        function obj = plotConfiguration(Dim, Solver)
            obj.Dim = Dim;
            obj.q_sys = Solver.Solution.q_sys;
            obj.t = Solver.Solution.t;
            obj.timenum = length(obj.t);
            obj.nodeNumber = numel(Solver.nodes);
            obj.XYZ = Solver.nodes{end}.q(1:3,end);
        end

        function getPositin(obj)
            obj.positiont = zeros(obj.Dim * obj.nodeNumber, obj.timenum);
            obj.position0 = zeros(obj.Dim * obj.nodeNumber, obj.timenum);
            coords = zeros(3 * obj.nodeNumber, obj.timenum);
            for i = 1:obj.timenum
                if obj.Dim == 2
                    for j = 1:obj.nodeNumber
                        obj.positiont(obj.Dim*(j-1)+1,i) = obj.q_sys(7*(j-1)+obj.pickDims(1),i);
                        obj.positiont(obj.Dim*(j-1)+2,i) = obj.q_sys(7*(j-1)+obj.pickDims(2),i);

                        obj.position0(obj.Dim*(j-1)+1,i) = obj.q_sys(7*(j-1)+obj.pickDims(1),1);
                        obj.position0(obj.Dim*(j-1)+2,i) = obj.q_sys(7*(j-1)+obj.pickDims(2),1);
                    end

                    if obj.rigidRotation == 1
                        
                        e = obj.q_sys(4:7,i);
                        A = Base.genA(e);
                        for k = 1:obj.nodeNumber
                            coords(3*k-2:3*k,i) = A * obj.q_sys(7*(k-1)+1:7*(k-1)+3,1);
                            obj.position0(obj.Dim*(k-1)+1,i) = coords(3*(k-1)+obj.pickDims(1),i);
                            obj.position0(obj.Dim*(k-1)+2,i) = coords(3*(k-1)+obj.pickDims(2),i);
                        end
                    end
                elseif obj.Dim == 3
                    for j = 1:obj.nodeNumber
                        obj.positiont(obj.Dim*(j-1)+1,i) = obj.q_sys(7*(j-1)+obj.pickDims(1),i);
                        obj.positiont(obj.Dim*(j-1)+2,i) = obj.q_sys(7*(j-1)+obj.pickDims(2),i); 
                        obj.positiont(obj.Dim*(j-1)+3,i) = obj.q_sys(7*(j-1)+obj.pickDims(3),i); 

                        obj.position0(obj.Dim*(j-1)+1,i) = obj.q_sys(7*(j-1)+obj.pickDims(1),1);
                        obj.position0(obj.Dim*(j-1)+2,i) = obj.q_sys(7*(j-1)+obj.pickDims(2),1);
                        obj.position0(obj.Dim*(j-1)+3,i) = obj.q_sys(7*(j-1)+obj.pickDims(3),1);
                    end

                    if obj.rigidRotation == 1
                        e = obj.q_sys(4:7,i);
                        A = Base.genA(e);
                        for k = 1:obj.nodeNumber
                            coords(3*k-2:3*k,i) = A * obj.q_sys(7*(k-1)+1:7*(k-1)+3,1);
                            obj.position0(obj.Dim*(k-1)+1,i) = coords(3*(k-1)+obj.pickDims(1),i);
                            obj.position0(obj.Dim*(k-1)+2,i) = coords(3*(k-1)+obj.pickDims(2),i);
                            obj.position0(obj.Dim*(k-1)+3,i) = coords(3*(k-1)+obj.pickDims(3),i);
                        end
                    end
                end            
            end
        end
        
        function set.Xlim(obj, val)
            obj.Xlim = val;
        end

        function set.Ylim(obj, val)
            obj.Ylim = val;
        end

        function set.Zlim(obj, val)
            obj.Zlim = val;
        end

        function set.rigidRotation(obj, val)
            obj.rigidRotation = val;
        end

        function set.xTicks(obj, val)
            obj.xTicks = val;
        end

        function set.yTicks(obj, val)
            obj.yTicks = val;
        end

        function set.zTicks(obj, val)
            obj.zTicks = val;
        end

        function set.pickDims(obj, val)
            assert(isvector(val), 'pickDims should be a 2x1 or 3x1 vector')
            obj.pickDims = val;
        end

        function set.xyzLabel(obj, stringCell)
            assert(iscell(stringCell), 'xyz label input should be cell array')
            obj.xyzLabel = stringCell;
        end

        function val = get.Xlim(obj)
            if isempty(obj.Xlim)
                val = [-1.2 * abs(obj.XYZ(1)), 1.2 * abs(obj.XYZ(1))];
            else
                val = obj.Xlim;
            end
        end

        function val = get.Ylim(obj)
            if isempty(obj.Ylim)
                val = [-1.2 * abs(obj.XYZ(2)), 1.2 * abs(obj.XYZ(2))];
            else
                val = obj.Ylim;
            end
        end

        function val = get.Zlim(obj)
            if isempty(obj.Zlim)
                val = [-1.2 * abs(obj.XYZ(3)), 1.2 * abs(obj.XYZ(3))];
            else
                val = obj.Zlim;
            end
        end

        function val = get.xTicks(obj)
            if isempty(obj.xTicks)
                val = linspace(obj.Xlim(1), obj.Xlim(2), 5);
            else
                val = obj.xTicks;
            end
        end

        function val = get.yTicks(obj)
            if isempty(obj.yTicks)
                val = linspace(obj.Ylim(1), obj.Ylim(2), 5);
            else
                val = obj.yTicks;
            end
        end

        function val = get.zTicks(obj)
            if isempty(obj.zTicks)
                val = linspace(obj.Zlim(1), obj.Zlim(2), 5);
            else
                val = obj.zTicks;
            end
        end

        function set.saveGif(obj, val)
            assert(val == 0 || val == 1, 'saveGif must be 0 or 1');
            obj.saveGif = val;
        end

        function plotBeamDyn(obj)
            obj.getPositin();
            if obj.Dim == 2
                figure;
                for i = 1:obj.timenum
                    X1 = obj.positiont(1:2:end,i);
                    Y1 = obj.positiont(2:2:end,i);
                    X2 = obj.position0(1:2:end,i);
                    Y2 = obj.position0(2:2:end,i);
    
                    plot(X1,Y1,'bo-','LineWidth',2,'MarkerSize',3)
                    hold on
                    plot(X2,Y2,'r*-','LineWidth',2,'MarkerSize',3)
                    hold off
                    set(gcf,'Color',[1,1,1]); % modify background color
                    set(gca,'xtick',obj.xTicks,'fontname','Times New Roman','fontsize',12);
                    set(gca,'ytick',obj.yTicks,'fontname','Times New Roman','fontsize',12);
                    xlabel(obj.xyzLabel{obj.pickDims(1)},'fontname','Times New Roman','fontsize',12)
                    ylabel(obj.xyzLabel{obj.pickDims(2)},'fontname','Times New Roman','fontsize',12)

                    timeText = text(gca, 0.02, 0.95, '', ...
                        'Units','normalized', ...
                        'FontSize',12, ...
                        'FontName','Times New Roman');
                    set(timeText, 'String', sprintf('t = %.4f s', obj.t(i)));

                    xlim(obj.Xlim)
                    ylim(obj.Ylim)
                    axis square;

                    drawnow; % Force update before capturing frame
                    
                    % --- GIF Saving Block ---
                    if obj.saveGif == 1
                        frame = getframe(gcf);
                        im = frame2im(frame);
                        [imind, cm] = rgb2ind(im, 256);
                        
                        % Write to the GIF file
                        if i == 1
                            imwrite(imind, cm, obj.gifFileName, 'gif', 'Loopcount', inf, 'DelayTime', obj.gifDelayTime);
                        else
                            imwrite(imind, cm, obj.gifFileName, 'gif', 'WriteMode', 'append', 'DelayTime', obj.gifDelayTime);
                        end
                    end
                    % ------------------------

                    pause(0.01)
                end
            elseif obj.Dim == 3
                figure;
                for i = 1:obj.timenum
                    X1 = obj.positiont(1:3:end,i);
                    Y1 = obj.positiont(2:3:end,i);
                    Z1 = obj.positiont(3:3:end,i);

                    X2 = obj.position0(1:3:end,i);
                    Y2 = obj.position0(2:3:end,i);
                    Z2 = obj.position0(3:3:end,i);

                    plot3(X1, Y1, Z1, 'bo-', ...
                        'LineWidth', 2, 'MarkerSize', 4);
                    hold on
                    plot3(X2, Y2, Z2, 'r*-', ...
                        'LineWidth', 2, 'MarkerSize', 4);
                    hold off

                    set(gcf,'Color',[1,1,1]);
                    set(gca,'xtick',obj.xTicks,'fontname','Times New Roman','fontsize',12);
                    set(gca,'ytick',obj.yTicks,'fontname','Times New Roman','fontsize',12);
                    set(gca,'ztick',obj.zTicks,'fontname','Times New Roman','fontsize',12);

                    axis equal
                    grid on
                    xlabel(obj.xyzLabel{obj.pickDims(1)},'FontName','Times New Roman')
                    ylabel(obj.xyzLabel{obj.pickDims(2)},'FontName','Times New Roman')
                    zlabel(obj.xyzLabel{obj.pickDims(3)},'FontName','Times New Roman')
    
                    timeText = text(gca, 0.02, 0.95, '', ...
                        'Units','normalized', ...
                        'FontSize',12, ...
                        'FontName','Times New Roman');
                    set(timeText, 'String', sprintf('t = %.4f s', obj.t(i)));

                    xlim(obj.Xlim)
                    ylim(obj.Ylim)
                    zlim(obj.Zlim)

                    view(3)        % key

                    drawnow; % Force update before capturing frame
                    
                    % --- GIF Saving Block ---
                    if obj.saveGif == 1
                        frame = getframe(gcf);
                        im = frame2im(frame);
                        [imind, cm] = rgb2ind(im, 256);
                        
                        % Write to the GIF file
                        if i == 1
                            imwrite(imind, cm, obj.gifFileName, 'gif', 'Loopcount', inf, 'DelayTime', obj.gifDelayTime);
                        else
                            imwrite(imind, cm, obj.gifFileName, 'gif', 'WriteMode', 'append', 'DelayTime', obj.gifDelayTime);
                        end
                    end
                    % ------------------------
                    pause(0.01)
                end
            else
                error('only support 2D or 3D plot')
            end
        end

        function plotBeamStatic(obj, ti)
            obj.getPositin();
            idx = find(obj.t == ti);
            if obj.Dim == 2
                figure;
                X1 = obj.positiont(1:2:end,idx);
                Y1 = obj.positiont(2:2:end,idx);
                X2 = obj.position0(1:2:end,idx);
                Y2 = obj.position0(2:2:end,idx);
    
                plot(X1,Y1,'bo-','LineWidth',2,'MarkerSize',3)
                hold on
                plot(X2,Y2,'r*-','LineWidth',2,'MarkerSize',3)
                hold off
                set(gcf,'Color',[1,1,1]); % modify background color
                set(gca,'xtick',obj.xTicks,'fontname','Times New Roman','fontsize',12);
                set(gca,'ytick',obj.yTicks,'fontname','Times New Roman','fontsize',12);
                xlabel(obj.xyzLabel{obj.pickDims(1)},'fontname','Times New Roman','fontsize',12)
                ylabel(obj.xyzLabel{obj.pickDims(2)},'fontname','Times New Roman','fontsize',12)

                timeText = text(gca, 0.02, 0.95, '', ...
                        'Units','normalized', ...
                        'FontSize',12, ...
                        'FontName','Times New Roman');
                set(timeText, 'String', sprintf('t = %.4f s', obj.t(idx)));

                xlim(obj.Xlim)
                ylim(obj.Ylim)
                axis square;
                pause(0.01)
            elseif obj.Dim == 3
                figure;
                X1 = obj.positiont(1:3:end,idx);
                Y1 = obj.positiont(2:3:end,idx);
                Z1 = obj.positiont(3:3:end,idx);

                X2 = obj.position0(1:3:end,idx);
                Y2 = obj.position0(2:3:end,idx);
                Z2 = obj.position0(3:3:end,idx);

                plot3(X1, Y1, Z1, 'bo-', ...
                        'LineWidth', 2, 'MarkerSize', 4);
                hold on
                plot3(X2, Y2, Z2, 'r*-', ...
                        'LineWidth', 2, 'MarkerSize', 4);
                hold off

                set(gcf,'Color',[1,1,1]);
                set(gca,'xtick',obj.xTicks,'fontname','Times New Roman','fontsize',12);
                set(gca,'ytick',obj.yTicks,'fontname','Times New Roman','fontsize',12);
                set(gca,'ztick',obj.zTicks,'fontname','Times New Roman','fontsize',12);

                axis equal
                grid on
                xlabel(obj.xyzLabel{obj.pickDims(1)},'FontName','Times New Roman')
                ylabel(obj.xyzLabel{obj.pickDims(2)},'FontName','Times New Roman')
                zlabel(obj.xyzLabel{obj.pickDims(3)},'FontName','Times New Roman')
    
                timeText = text(gca, 0.02, 0.95, '', ...
                        'Units','normalized', ...
                        'FontSize',12, ...
                        'FontName','Times New Roman');
                set(timeText, 'String', sprintf('t = %.4f s', obj.t(idx)));

                xlim(obj.Xlim)
                ylim(obj.Ylim)
                zlim(obj.Zlim)

                view(3)        % key
                drawnow
                pause(0.01)
            else
                error('only support 2D or 3D plot')
            end
        end
    end
end

