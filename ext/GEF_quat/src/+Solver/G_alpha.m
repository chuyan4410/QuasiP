classdef G_alpha < matlab.mixin.SetGetExactNames
    % class for generalized-alpha scheme for time integration
    
    properties
        elements            % element class array
        nodes               % node class array
        constraints         % constraints class
        Dconstraints        % displacement constraints class
        Fconstraints        % force constraints class

        %------Generalized-alpha parameters--------------------------------
        rho_s = 0           % spectral radius: [0, 1]. rho_s = 1: Newmark
        alpha_f
        alpha_m
        gamma
        beta
        gamma_hat
        beta_hat

        %------integration quantities--------------------------------------
        dynamic = 1         % dynamical problem: 1; static problem: 0
        dof_sys             % number of DOF for mechanical system
        nodenum             % number of nodes
        dof_constraint      % number of DOF for constraint (no quaternion)
        t                   % sequence of time steps
        timenum             % number of time steps
        dt                  % time step size
        Solution            % solution structure data
        max_it      = 100   % maximum number of Newton iteration 
        tol         = 1e-6  % tolerance of Newton iteration
        itNum       = []    % storage for numbers of Newton iteration 

        DL                  % left pre-processing matrix
        DR                  % right pre-processing matrix
        M                   % systematical mass matrix
        Ct                  % systematical gyro matrix 
        Kt                  % systematical tangent stiffness matrix
        Fi                  % systematical force vector corresponding to inertia 
        Fg                  % systematical body force vector
        Frot                % extra inertia force caused by rotating frame
        Fe                  % systematical elastic force vector
        Jac                 % Jacobi matrix in DAEs
        Res                 % Residual of governing equations

        qti_sys             % current state position vector
        dqti_sys            % current state velocity vector
        ddqti_sys           % current state acceleration vector
        lmti_sys            % current Lagrange multiplier
        ati_sys             % current auxiliary acceleration vector
        Zic = []            % initial condition
        initialAuxiliaryAcceleration = [] % restart value for Solution.a(:,1)

        dampRatio = 0       % ratio of Rayleigh damping
        Damp = []           % Rayleigh damping matrix

        rotationFrame = struct('enabled', false) % prescribed rotating frame data
        rotFrameCurrent = struct('enabled', false)
        aeroForceModel = []
        aeroForceUpdateMode = 'timeStep'
        aeroForceCache = struct('stepID', [], 'forceSys', [], 'forceOut', [])
        aeroForceHistory = {}
    end
    
    methods
        function obj = G_alpha(elements, nodes, constraints, dt, dof_sys)
            obj.elements       = elements;
            obj.nodes          = nodes;
            obj.constraints    = constraints;
            obj.Dconstraints   = constraints.Dconstraints;
            obj.Fconstraints   = constraints.Fconstraints;

            obj.dof_constraint = obj.Dconstraints.dof_constraint;
            obj.dt             = dt;
            obj.dof_sys        = dof_sys;
            obj.nodenum        = numel(nodes);
        end

        function set.rho_s(obj, rho_s)
            obj.rho_s = rho_s;
        end

        function set.dampRatio(obj, dampRatio)
            obj.dampRatio = dampRatio;
        end

        function set.Damp(obj, Damp)
            obj.Damp = Damp;
        end

        function set.rotationFrame(obj, rotationFrame)
            obj.rotationFrame = rotationFrame;
        end

        function val = get.alpha_f(obj)
            val = obj.rho_s / (obj.rho_s + 1);
        end

        function val = get.alpha_m(obj)
            val = (2 * obj.rho_s - 1) / (obj.rho_s + 1);
        end

        function val = get.gamma(obj)
            val = 1/2 - obj.alpha_m + obj.alpha_f;
        end

        function val = get.beta(obj)
            val = 1/4 * (1 - obj.alpha_m + obj.alpha_f)^2;
        end

        function val = get.gamma_hat(obj)
            val = obj.gamma / obj.dt / obj.beta;
        end

        function val = get.beta_hat(obj)
            val = (1 - obj.alpha_m) / obj.dt^2 / obj.beta / (1 - obj.alpha_f);
        end
        
        function set.max_it(obj, max_it)
            obj.max_it = max_it;
        end
        
        function set.tol(obj, tol)
            obj.tol = tol;
        end

        function set.dynamic(obj, dynamic)
            obj.dynamic = dynamic;
        end

        function Integrate(obj, t)
            [q0_sys, dq0_sys, ddq0_sys] = Base.sysInitCond(obj.nodes);
            obj.timenum    = length(t);
            obj.Solution.t = t;
            obj.Solution.q_sys   = zeros(obj.dof_sys, obj.timenum);
            obj.Solution.dq_sys  = zeros(obj.dof_sys, obj.timenum);
            obj.Solution.ddq_sys = zeros(obj.dof_sys, obj.timenum);
            obj.Solution.lm      = zeros(obj.dof_constraint + obj.nodenum, obj.timenum);
            obj.Solution.a       = zeros(obj.dof_sys, obj.timenum);
            obj.Solution.q_sys(:,1)   = q0_sys;
            obj.Solution.dq_sys(:,1)  = dq0_sys;
            obj.Solution.ddq_sys(:,1) = ddq0_sys;
            if ~isempty(obj.initialAuxiliaryAcceleration)
                if numel(obj.initialAuxiliaryAcceleration) ~= obj.dof_sys
                    error('Solver:InvalidInitialAuxiliaryAcceleration', ...
                        'initialAuxiliaryAcceleration must have dof_sys entries.');
                end
                obj.Solution.a(:,1) = obj.initialAuxiliaryAcceleration(:);
            end
            
            startStep = 2;
            if ~isempty(obj.Zic)
                startStep = 3;
                obj.Solution.q_sys(:,2)   = obj.Zic(1:obj.dof_sys);
                obj.Solution.dq_sys(:,2)  = obj.Zic(obj.dof_sys+1 : 2*obj.dof_sys);
            end
            
            obj.preProcess();
            if isempty(obj.Damp)
                obj.Damp = Base.genDamp(obj.elements, obj.nodes, ...
                                        obj.constraints, obj.dampRatio);
            end

            for i = startStep:obj.timenum
                obj.updateRotationFrame(i);
                % prediction
                obj.qti_sys  = obj.Solution.q_sys(:,i-1) + obj.dt * obj.Solution.dq_sys(:,i-1) + ...
                           obj.dt^2 * (0.5 - obj.beta) * obj.Solution.a(:,i-1);
                obj.dqti_sys = obj.Solution.dq_sys(:,i-1) + obj.dt * (1 - obj.gamma) * obj.Solution.a(:,i-1);
                obj.lmti_sys = zeros(obj.dof_constraint + obj.nodenum,1);
                obj.ati_sys  = (obj.alpha_f * obj.Solution.ddq_sys(:,i-1) ...
                              - obj.alpha_m * obj.Solution.a(:,i - 1)) / (1 - obj.alpha_m);
                obj.qti_sys  = obj.qti_sys + obj.dt^2 * obj.beta * obj.ati_sys;
                obj.dqti_sys = obj.dqti_sys + obj.dt * obj.gamma * obj.ati_sys;
                obj.ddqti_sys= zeros(obj.dof_sys,1);

                obj.nodeUpdate(i);
                obj.prepareAeroForceModel(i);
                % Newton iteration
                k = 0;
                while k <= obj.max_it
                    %% assmble and compute the element matrices and vectors
                    obj.assemble();

                    [obj.Jac, obj.Res] = obj.Dconstraints.change_K_P(obj, i);
                    [obj.Jac, obj.Res] = obj.Fconstraints.change_K_P(obj, i);
                    obj.applyAeroForceModel(i);

                    JacPre   = obj.DL * obj.Jac * obj.DR;
                    rPre     = obj.DL * obj.Res;
                    DeltaPre =-JacPre \ rPre;
                    Delta    = obj.DR * DeltaPre;
                    delta_q  = Delta(1:obj.dof_sys);
                    delta_lm = Delta(obj.dof_sys+1:end);
                    obj.qti_sys   = obj.qti_sys   + delta_q;
                    obj.dqti_sys  = obj.dqti_sys  + delta_q * obj.gamma_hat;
                    obj.ddqti_sys = obj.ddqti_sys + delta_q * obj.beta_hat;
                    obj.lmti_sys  = obj.lmti_sys  + delta_lm;

                    dif = norm(obj.Res);
                    if dif < obj.tol
                        obj.itNum = [obj.itNum; k];
                        break
                    elseif k == obj.max_it
                        error('Newton iteration does not converge')
                    end
                    %% update node quantities in iteration steps
                    obj.nodeUpdate(i);
                    k = k + 1;
                end
                obj.Solution.a(:,i)       = obj.ati_sys + (1 - obj.alpha_f)/(1 - obj.alpha_m) * obj.ddqti_sys;
                obj.Solution.q_sys(:,i)   = obj.qti_sys;
                obj.Solution.dq_sys(:,i)  = obj.dqti_sys;
                obj.Solution.ddq_sys(:,i) = obj.ddqti_sys;
                obj.Solution.lm(:,i)      = obj.lmti_sys;

                disp(['TimeStep: ',num2str(i),', ItNum: ',num2str(k)])
                tNow = datetime('now','Format','dd:HH:mm:ss.SSS');
                disp([char(tNow), '  Time:', num2str(t(i)), 's is finished']);
            end
        end

        function nodeUpdate(obj, i)
            for n = 1:numel(obj.nodes)
                id = obj.nodes{n}.dofID;
                qti_node   = obj.qti_sys(id);
                dqti_node  = obj.dqti_sys(id);
                ddqti_node = obj.ddqti_sys(id);
                obj.nodes{n}.updateNode(i, qti_node, dqti_node, ddqti_node);
            end
        end

        function assemble(obj)
            n = obj.dof_sys;
            obj.M  = zeros(n,n);
            obj.Ct = zeros(n,n);
            obj.Kt = zeros(n,n);
            obj.Fi = zeros(n,1);
            obj.Fg = zeros(n,1);
            obj.Frot = zeros(n,1);
            obj.Fe = zeros(n,1);

            for e = 1:numel(obj.elements)
                elm = obj.elements{e};
                if obj.isRotationFrameEnabled()
                    elm.setRotationFrame(obj.rotFrameCurrent);
                else
                    elm.setRotationFrame(struct('enabled', false));
                end

                elm.get_K_P();

                id = elm.dofID;

                obj.M(id,id)  = obj.M(id,id)  + elm.M;
                obj.Ct(id,id) = obj.Ct(id,id) + elm.Ct;
                obj.Kt(id,id) = obj.Kt(id,id) + elm.Kt;
                obj.Fi(id)    = obj.Fi(id)    + elm.Fi;
                obj.Fg(id)    = obj.Fg(id)    + elm.Fg;
                obj.Frot(id)  = obj.Frot(id)  + elm.Frot;
                obj.Fe(id)    = obj.Fe(id)    + elm.Fe;
            end
        end

        function prepareAeroForceModel(obj, i)
            if isempty(obj.aeroForceModel) || strcmpi(obj.aeroForceUpdateMode, 'newton')
                return
            end
            if ~strcmpi(obj.aeroForceUpdateMode, 'timeStep')
                error('Solver:InvalidAeroForceUpdateMode', ...
                    'aeroForceUpdateMode must be ''timeStep'' or ''newton''.');
            end

            [forceSys, forceOut] = obj.aeroForceModel.evaluate(obj, i);
            obj.validateAeroForce(forceSys);
            obj.aeroForceCache = struct('stepID', i, ...
                'forceSys', forceSys(:), 'forceOut', forceOut);
            obj.storeAeroForceHistory(i, forceOut);
        end

        function applyAeroForceModel(obj, i)
            if isempty(obj.aeroForceModel)
                return
            end

            if strcmpi(obj.aeroForceUpdateMode, 'newton')
                [forceSys, forceOut] = obj.aeroForceModel.evaluate(obj, i);
                obj.validateAeroForce(forceSys);
                obj.storeAeroForceHistory(i, forceOut);
            elseif strcmpi(obj.aeroForceUpdateMode, 'timeStep')
                if isempty(obj.aeroForceCache.stepID) || obj.aeroForceCache.stepID ~= i
                    obj.prepareAeroForceModel(i);
                end
                forceSys = obj.aeroForceCache.forceSys;
            else
                error('Solver:InvalidAeroForceUpdateMode', ...
                    'aeroForceUpdateMode must be ''timeStep'' or ''newton''.');
            end

            obj.Res(1:obj.dof_sys) = obj.Res(1:obj.dof_sys) - forceSys(:);
        end

        function validateAeroForce(obj, forceSys)
            if numel(forceSys) ~= obj.dof_sys
                error('Solver:InvalidAeroForce', ...
                    'The aero force model must return a vector with dof_sys entries.');
            end
        end

        function storeAeroForceHistory(obj, i, forceOut)
            if numel(obj.aeroForceHistory) < obj.timenum
                obj.aeroForceHistory{obj.timenum} = [];
            end
            obj.aeroForceHistory{i} = forceOut;
        end

        function tf = isRotationFrameEnabled(obj)
            tf = isstruct(obj.rotationFrame) && isfield(obj.rotationFrame, 'enabled') ...
                && obj.rotationFrame.enabled;
        end

        function updateRotationFrame(obj, i)
            if ~obj.isRotationFrameEnabled()
                obj.rotFrameCurrent = struct('enabled', false);
                return
            end

            er   = obj.getRotationFrameVector('e', i);
            der  = obj.getRotationFrameVector('de', i);
            dder = obj.getRotationFrameVector('dde', i);

            obj.rotFrameCurrent.enabled  = true;
            obj.rotFrameCurrent.e        = er;
            obj.rotFrameCurrent.de       = der;
            obj.rotFrameCurrent.dde      = dder;
            obj.rotFrameCurrent.Lambda   = Base.genA(er);
            obj.rotFrameCurrent.Omega    = 2 * Base.genG(er) * der;
            obj.rotFrameCurrent.dotOmega = 2 * Base.genG(er) * dder;
        end

        function v = getRotationFrameVector(obj, fieldName, i)
            if ~isfield(obj.rotationFrame, fieldName)
                if strcmp(fieldName, 'e')
                    v = [1; 0; 0; 0];
                elseif isfield(obj.rotationFrame, 'e')
                    v = obj.estimateRotationFrameDerivative(fieldName, i);
                else
                    v = zeros(4,1);
                end
                return
            end

            data = obj.rotationFrame.(fieldName);
            if iscell(data)
                v = data{i};
            elseif ismatrix(data)
                if size(data, 1) == 4
                    v = data(:,i);
                elseif size(data, 2) == 4
                    v = data(i,:).';
                else
                    error('RotationFrame:InvalidInput', ...
                          'rotationFrame.%s must have one dimension equal to 4.', fieldName);
                end
            else
                error('RotationFrame:InvalidInput', ...
                      'rotationFrame.%s must be a cell array or numeric matrix.', fieldName);
            end
            v = v(:);
        end

        function v = estimateRotationFrameDerivative(obj, fieldName, i)
            if strcmp(fieldName, 'de')
                if i == 1
                    v = (obj.getRotationFrameVector('e', 2) - obj.getRotationFrameVector('e', 1)) / obj.dt;
                elseif i == obj.timenum
                    v = (obj.getRotationFrameVector('e', i) - obj.getRotationFrameVector('e', i-1)) / obj.dt;
                else
                    v = (obj.getRotationFrameVector('e', i+1) - obj.getRotationFrameVector('e', i-1)) / (2 * obj.dt);
                end
            elseif strcmp(fieldName, 'dde')
                if i == 1 || i == obj.timenum
                    v = zeros(4,1);
                else
                    v = (obj.getRotationFrameVector('e', i+1) - 2 * obj.getRotationFrameVector('e', i) ...
                       + obj.getRotationFrameVector('e', i-1)) / obj.dt^2;
                end
            else
                v = zeros(4,1);
            end
        end

        function preProcess(obj)
            % dof_constraint: the number of constrained DOF, without number of
            % nodes (intrinsic constraint for quaternion)
            if obj.dynamic == true
                obj.DL = zeros(obj.dof_sys + obj.dof_constraint + obj.nodenum);
                obj.DR = zeros(obj.dof_sys + obj.dof_constraint + obj.nodenum);
                obj.DL(1:obj.dof_sys,1:obj.dof_sys) = obj.beta * obj.dt^2 * eye(obj.dof_sys);
                obj.DL(obj.dof_sys+1:end,obj.dof_sys+1:end) = eye(obj.dof_constraint + obj.nodenum);
                obj.DR(1:obj.dof_sys,1:obj.dof_sys) = eye(obj.dof_sys);
                obj.DR(obj.dof_sys+1:end,obj.dof_sys+1:end) = (obj.beta * obj.dt^2) \ eye(obj.dof_constraint + obj.nodenum);
            else
                obj.DL = eye(obj.dof_sys + obj.dof_constraint + obj.nodenum);
                obj.DR = eye(obj.dof_sys + obj.dof_constraint + obj.nodenum);
            end
        end
    end
end

