classdef Assembly < matlab.mixin.SetGetExactNames
    %ASSEMBLY : assemble the systematical matrices and vectors
    % Use this class to extract the systematical matrices and vectors
    % at specific time 
    
    properties
        dof_sys             % number of DOFs for mechanical system
        dof_constraint      % number of constrained DOF
        qti_sys             % current position vector
        dqti_sys            % current velocity vector
        ddqti_sys           % current acceleration vector
        nonlinOrder         % order of nonlinearity
        rotationFrame = struct('enabled', false)
        rotFrameCurrent = struct('enabled', false)
        rotationFrameStep = 2

        elements            % element classes
        nodes               % node classes
        Dconstraints        % displacement constraint classes
        Fconstraints        % force constraint classes

        %------systematical matrices and vectors--------------------------------
        M                   % systematical mass matrix
        M0                  % constant part of systematical mass matrix
        Mnl                 % nonlinear part of systematical mass matrix
        Ct                  % systematical gyro matrix 
        Kt                  % systematical tangent stiffness matrix
        Fi                  % systematical force vector corresponding to inertia 
        Fg                  % systematical body force vector
        Frot                % extra inertia force caused by rotating frame
        Fe                  % systematical elastic force vector
        Fe1                 % 1st-order part of systematical elastic force vector
        Fe2                 % 2nd-order part of systematical elastic force vector
        Fe3                 % 3rd-order part of systematical elastic force vector
        Fe5                 % 5th-order part of systematical elastic force vector
        Fext                % systematical external force vector 
        
        %------about the constraints--------------------------------
        C                   % the constraint equations
        Cq                  % the constraint Jacobian
        Kb                  % the added stiffness matrix of constraint
    end
    
    methods
        function obj = Assembly(elements, nodes, constraints, qti_sys, dqti_sys, ddqti_sys)
            obj.elements     = elements;
            obj.nodes        = nodes;
            obj.Dconstraints = constraints.Dconstraints;
            obj.Fconstraints = constraints.Fconstraints;
            obj.dof_sys      = length(qti_sys);
            obj.qti_sys      = qti_sys;
            obj.dqti_sys     = dqti_sys;
            obj.ddqti_sys    = ddqti_sys;
        end

        function val = get.dof_constraint(obj)
            val = obj.Dconstraints.dof_constraint;
        end

        function set.rotationFrame(obj, rotationFrame)
            obj.rotationFrame = rotationFrame;
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
            if isfield(obj.rotationFrame, 'dt')
                dt = obj.rotationFrame.dt;
            elseif isfield(obj.rotationFrame, 't')
                t = obj.rotationFrame.t;
                if i > 1 && i < numel(t)
                    dt = (t(i+1) - t(i-1)) / 2;
                elseif i > 1
                    dt = t(i) - t(i-1);
                else
                    dt = t(2) - t(1);
                end
            else
                error('RotationFrame:MissingDerivative', ...
                      ['rotationFrame.%s is missing. Provide de/dde explicitly, ' ...
                       'or provide rotationFrame.dt or rotationFrame.t for finite differences.'], fieldName);
            end

            numSteps = obj.getRotationFrameStepCount();
            if strcmp(fieldName, 'de')
                if i == 1
                    v = (obj.getRotationFrameVector('e', 2) - obj.getRotationFrameVector('e', 1)) / dt;
                elseif i == numSteps
                    v = (obj.getRotationFrameVector('e', i) - obj.getRotationFrameVector('e', i-1)) / dt;
                else
                    v = (obj.getRotationFrameVector('e', i+1) - obj.getRotationFrameVector('e', i-1)) / (2 * dt);
                end
            elseif strcmp(fieldName, 'dde')
                if i == 1 || i == numSteps
                    v = zeros(4,1);
                else
                    v = (obj.getRotationFrameVector('e', i+1) - 2 * obj.getRotationFrameVector('e', i) ...
                       + obj.getRotationFrameVector('e', i-1)) / dt^2;
                end
            else
                v = zeros(4,1);
            end
        end

        function nSteps = getRotationFrameStepCount(obj)
            data = obj.rotationFrame.e;
            if iscell(data)
                nSteps = numel(data);
            elseif size(data,1) == 4
                nSteps = size(data,2);
            elseif size(data,2) == 4
                nSteps = size(data,1);
            else
                error('RotationFrame:InvalidInput', ...
                      'rotationFrame.e must have one dimension equal to 4.');
            end
        end

        function updateNodesFromState(obj, stepID, q_sys, dq_sys, ddq_sys)
            for n = 1:numel(obj.nodes)
                id = obj.nodes{n}.dofID;
                obj.nodes{n}.updateNode(stepID, q_sys(id), dq_sys(id), ddq_sys(id));
            end
        end

        function setElementRotationFrame(obj, elm)
            if obj.isRotationFrameEnabled()
                elm.setRotationFrame(obj.rotFrameCurrent);
            else
                elm.setRotationFrame(struct('enabled', false));
            end
        end
        
        function genSysMatVec(obj, i)
            % i: the i-th step of time
            obj.updateRotationFrame(i);
            n = obj.dof_sys;
            obj.M  = sparse(n,n);
            obj.Ct = sparse(n,n);
            obj.Kt = sparse(n,n);
            obj.Fi = sparse(n,1);
            obj.Fg = sparse(n,1);
            obj.Frot = sparse(n,1);
            obj.Fe = sparse(n,1);
            obj.Fext = sparse(n,1);

            obj.updateNodesFromState(i, obj.qti_sys, obj.dqti_sys, obj.ddqti_sys);
            
            for e = 1:numel(obj.elements)
                elm = obj.elements{e};
                obj.setElementRotationFrame(elm);

                elm.get_K_P();
                if obj.isRotationFrameEnabled()
                    CrotApprox = elm.Crot;
                    KrotApprox = elm.Krot;
                    elm.get_consistentRotFrameTangents();
                    elm.Ct = elm.Ct - CrotApprox + elm.Crot;
                    elm.Kt = elm.Kt - KrotApprox + elm.Krot;
                end
                elm.Kt = elm.Kt + elm.genClassicInertiaKt();

                id = elm.dofID;

                obj.M(id,id)  = obj.M(id,id)  + elm.M;
                obj.Ct(id,id) = obj.Ct(id,id) + elm.Ct;
                obj.Kt(id,id) = obj.Kt(id,id) + elm.Kt;
                obj.Fi(id)    = obj.Fi(id)    + elm.Fi;
                obj.Fg(id)    = obj.Fg(id)    + elm.Fg;
                obj.Frot(id)  = obj.Frot(id)  + elm.Frot;
                obj.Fe(id)    = obj.Fe(id)    + elm.Fe;
            end

            if strcmpi(obj.Dconstraints.kind{1}{1}, 'displacement boundary')
                load = obj.Dconstraints.loadList{i};
                obj.C  = sparse(obj.dof_constraint + numel(obj.nodes), 1);
                obj.Cq = sparse(obj.dof_constraint + numel(obj.nodes), obj.dof_sys);
                obj.Kb = sparse(obj.dof_sys, obj.dof_sys);
                kk = 1;
                for ii = 1:numel(obj.Dconstraints.dof)
                    nodeID  = obj.Dconstraints.dof{ii}{1};
                    dofList = obj.Dconstraints.dof{ii}{2};
                    dofID   = obj.nodes{nodeID}.dofID;
                    for jj = dofList
                        obj.C(kk) = obj.qti_sys(dofID(jj)) - load{1}(kk);
                        kk = kk + 1; 
                    end
                end
                obj.Cq(1:kk-1,:) = load{2};

                % special operation for gradient discontinuity
                if strcmpi(obj.Dconstraints.kind{2}{1}, 'gradient discontinuity')
                    nodeA = obj.Dconstraints.kind{2}{2}{1};
                    nodeB = obj.Dconstraints.kind{2}{2}{2};
                    e_rel = obj.Dconstraints.kind{2}{2}{3};
                    R_rel = Base.genQuatRightMultiMatrix(e_rel);
                    dofID_A = obj.nodes{nodeA}.dofID;
                    dofID_B = obj.nodes{nodeB}.dofID;
                    obj.C(kk:kk+2) = obj.qti_sys(dofID_B(1:3)) - obj.qti_sys(dofID_A(1:3));
                    e_A = obj.qti_sys(dofID_A(4:7));
                    e_B = obj.qti_sys(dofID_B(4:7));
                    tmp = Base.quatMultiply(e_A, e_rel);
                    e_constraint = e_B - tmp;
                    obj.C(kk+3:kk+5) = e_constraint(2:4);
                    obj.Cq(kk:kk+2,dofID_A(1:3))   =-eye(3);
                    obj.Cq(kk:kk+2,dofID_B(1:3))   = eye(3);
                    obj.Cq(kk+3:kk+5,dofID_A(4:7)) =-R_rel(2:4,1:4);
                    obj.Cq(kk+3:kk+5,dofID_B(5:7)) = eye(3);
                end

                for n = 1:numel(obj.nodes)
                    obj.C(obj.dof_constraint+n) = sum(obj.qti_sys(7 * n - 3 : 7 * n).^2) - 1;
                    obj.Cq(obj.dof_constraint+n,obj.nodes{n}.dofID(4:7)) = 2 * obj.qti_sys(obj.nodes{n}.dofID(4:7));
                    % obj.Kb(7*n-3:7*n, 7*n-3:7*n) = 2 * obj.lmti(n+obj.dof_constraint) * eye(4);
                end
                obj.Kb = obj.Kb + load{3};
            end

            if strcmpi(obj.Fconstraints.kind{1}{1}, 'force boundary')
                load = obj.Fconstraints.loadList{i};
                for ii = 1:numel(obj.Fconstraints.dof)
                    nodeID   = obj.Fconstraints.dof{ii}{1};
                    dofID    = obj.nodes{nodeID}.dofID;
                    loadVec = load{ii};
                    if obj.isRotationFrameEnabled()
                        loadVec = [obj.rotFrameCurrent.Lambda.' * loadVec(1:3); ...
                                   obj.rotFrameCurrent.Lambda.' * loadVec(4:6)];
                    end
                    quatLoad = 2 * obj.nodes{nodeID}.H{i}.' * loadVec(4:6);

                    obj.Fext(dofID(1:3)) = obj.Fext(dofID(1:3)) + loadVec(1:3);
                    obj.Fext(dofID(4:7)) = obj.Fext(dofID(4:7)) + quatLoad;
                end
            end
        end

        function genSysForceTangents(obj, zti_sys)
            %GENSYSFORCETANGENTS Assemble force values and exact tangents once.
            %
            % This path is intended for shooting/variational calculations.
            % At a fixed [q;dq], it evaluates M, Fe, Fi, Fg, Frot and the
            % complete Kt/Ct in one element/Gauss-point traversal. The
            % regular time integrator continues to use genSysMatVec/get_K_P
            % with its existing tangent convention.
            obj.updateRotationFrame(obj.rotationFrameStep);
            n = obj.dof_sys;
            assert(numel(zti_sys)==2*n, ...
                'Assembly:InvalidStateSize', ...
                'zti_sys must contain [q;dq] with 2*dof_sys entries.');

            obj.qti_sys = zti_sys(1:n);
            obj.dqti_sys = zti_sys(n+1:2*n);
            obj.ddqti_sys = zeros(n,1);
            obj.M = sparse(n,n);
            obj.Ct = sparse(n,n);
            obj.Kt = sparse(n,n);
            obj.Fi = sparse(n,1);
            obj.Fg = sparse(n,1);
            obj.Frot = sparse(n,1);
            obj.Fe = sparse(n,1);

            obj.updateNodesFromState(obj.rotationFrameStep, ...
                obj.qti_sys,obj.dqti_sys,obj.ddqti_sys);

            for e = 1:numel(obj.elements)
                elm = obj.elements{e};
                obj.setElementRotationFrame(elm);
                elm.get_K_P(true);
                elm.Kt = elm.Kt+elm.genClassicInertiaKt();

                id = elm.dofID;
                obj.M(id,id) = obj.M(id,id)+elm.M;
                obj.Ct(id,id) = obj.Ct(id,id)+elm.Ct;
                obj.Kt(id,id) = obj.Kt(id,id)+elm.Kt;
                obj.Fi(id) = obj.Fi(id)+elm.Fi;
                obj.Fg(id) = obj.Fg(id)+elm.Fg;
                obj.Frot(id) = obj.Frot(id)+elm.Frot;
                obj.Fe(id) = obj.Fe(id)+elm.Fe;
            end
        end

        function M0 = genSysM0(obj)
            n = obj.dof_sys;
            obj.M0 = sparse(n,n);

            for n = 1:numel(obj.nodes)
                id = obj.nodes{n}.dofID;
                qti_node   = obj.qti_sys(id);
                dqti_node  = obj.dqti_sys(id);
                ddqti_node = obj.ddqti_sys(id);
                obj.nodes{n}.updateNode(2, qti_node, dqti_node, ddqti_node);
            end

            for e = 1:numel(obj.elements)
                elm = obj.elements{e};

                elm.get_M0();

                id = elm.dofID;

                obj.M0(id,id) = obj.M0(id,id) + elm.M0;
            end
            M0 = sparse(obj.M0);
        end

        function Mnl = genSysMnl(obj, zti_sys)
            n = obj.dof_sys;
            obj.Mnl = sparse(n,n);
            obj.nonlinOrder.Mnl = 2;

            xti_sys   = zti_sys(1:n);
            dxti_sys  = zti_sys(n+1:2*n);
            ddxti_sys = zeros(n,1);

            for n = 1:numel(obj.nodes)
                id = obj.nodes{n}.dofID;
                qti_node   = xti_sys(id);
                dqti_node  = dxti_sys(id);
                ddqti_node = ddxti_sys(id);
                obj.nodes{n}.updateNode(2, qti_node, dqti_node, ddqti_node);
            end

            for e = 1:numel(obj.elements)
                elm = obj.elements{e};

                elm.get_Mnl();

                id = elm.dofID;

                obj.Mnl(id,id) = obj.Mnl(id,id) + elm.Mnl;
            end
            Mnl = sparse(obj.Mnl);
        end

        function dM_dq = genSysdMdq(obj, zti_sys)
            NumGlobalDof = obj.dof_sys;
            xti_sys   = zti_sys(1:NumGlobalDof);
            dxti_sys  = zti_sys(NumGlobalDof+1:2*NumGlobalDof);
            ddxti_sys = zeros(NumGlobalDof,1);

            NumElem     = numel(obj.elements);
            DofsPerElem = obj.elements{1}.numDOF;

            EntriesPerElem = DofsPerElem^3;
            TotalEntries   = NumElem * EntriesPerElem;

            I_list = zeros(TotalEntries, 1); % global index i
            J_list = zeros(TotalEntries, 1); % gloabl index j
            K_list = zeros(TotalEntries, 1); % gloabl index k (derivative w.r.t. q_k)
            V_list = zeros(TotalEntries, 1); % corresponding value

            [local_j, local_i, local_k] = meshgrid(1:DofsPerElem, 1:DofsPerElem, 1:DofsPerElem);
            local_i_vec = local_i(:);
            local_j_vec = local_j(:);
            local_k_vec = local_k(:);

            for n = 1:numel(obj.nodes)
                id = obj.nodes{n}.dofID;
                qti_node   = xti_sys(id);
                dqti_node  = dxti_sys(id);
                ddqti_node = ddxti_sys(id);
                obj.nodes{n}.updateNode(2, qti_node, dqti_node, ddqti_node);
            end

            ptr = 0;
            for e = 1:numel(obj.elements)
                elm = obj.elements{e};

                elm.get_dMdq();

                id = elm.dofID;
                global_I = id(local_i_vec).'; 
                global_J = id(local_j_vec).';
                global_K = id(local_k_vec).';

                range = ptr + (1:EntriesPerElem);
                I_list(range) = global_I;
                J_list(range) = global_J;
                K_list(range) = global_K;
                V_list(range) = elm.dM_dq(:);
                
                ptr = ptr + EntriesPerElem;
            end
            subs  = [I_list, J_list, K_list];
            dM_dq = sptensor(subs, V_list, [NumGlobalDof, NumGlobalDof, NumGlobalDof]);
        end

        function Ct = genSysCt(obj, zti_sys)
            obj.updateRotationFrame(obj.rotationFrameStep);
            n = obj.dof_sys;
            obj.Ct = sparse(n,n);

            xti_sys   = zti_sys(1:n);
            dxti_sys  = zti_sys(n+1:2*n);
            ddxti_sys = zeros(n,1);

            obj.updateNodesFromState(obj.rotationFrameStep, xti_sys, dxti_sys, ddxti_sys);

            for e = 1:numel(obj.elements)
                elm = obj.elements{e};
                obj.setElementRotationFrame(elm);

                elm.get_Ct();

                id = elm.dofID;

                obj.Ct(id,id) = obj.Ct(id,id) + elm.Ct;
            end
            Ct = sparse(obj.Ct);
        end

        function Kt = genSysKt(obj, zti_sys)
            obj.updateRotationFrame(obj.rotationFrameStep);
            n = obj.dof_sys;
            obj.Kt = sparse(n,n);

            xti_sys   = zti_sys(1:n);
            dxti_sys  = zti_sys(n+1:2*n);
            ddxti_sys = zeros(n,1);

            obj.updateNodesFromState(obj.rotationFrameStep, xti_sys, dxti_sys, ddxti_sys);

            for e = 1:numel(obj.elements)
                elm = obj.elements{e};
                obj.setElementRotationFrame(elm);

                elm.get_Kt();
                elm.Kt = elm.Kt + elm.genClassicInertiaKt();

                id = elm.dofID;

                obj.Kt(id,id) = obj.Kt(id,id) + elm.Kt;
            end
            Kt = sparse(obj.Kt);
        end

        function Fi = genSysFi(obj, zti_sys)
            n = obj.dof_sys;
            obj.Fi = sparse(n,1);
            obj.nonlinOrder.Fi = 3;

            xti_sys   = zti_sys(1:n);
            dxti_sys  = zti_sys(n+1:2*n);
            ddxti_sys = zeros(n,1);

            for n = 1:numel(obj.nodes)
                id = obj.nodes{n}.dofID;
                qti_node   = xti_sys(id);
                dqti_node  = dxti_sys(id);
                ddqti_node = ddxti_sys(id);
                obj.nodes{n}.updateNode(2, qti_node, dqti_node, ddqti_node);
            end
            
            for e = 1:numel(obj.elements)
                elm = obj.elements{e};

                elm.get_Fi();

                id = elm.dofID;

                obj.Fi(id) = obj.Fi(id) + elm.Fi;
            end
            Fi = sparse(obj.Fi);
        end

        function Fg = genSysFg(obj, zti_sys)
            obj.updateRotationFrame(obj.rotationFrameStep);
            n = obj.dof_sys;
            obj.Fg = sparse(n,1);

            xti_sys   = zti_sys(1:n);
            dxti_sys  = zti_sys(n+1:2*n);
            ddxti_sys = zeros(n,1);

            obj.updateNodesFromState(obj.rotationFrameStep, xti_sys, dxti_sys, ddxti_sys);

            for e = 1:numel(obj.elements)
                elm = obj.elements{e};
                obj.setElementRotationFrame(elm);

                elm.get_K_P();

                id = elm.dofID;
                obj.Fg(id) = obj.Fg(id) + elm.Fg;
            end
            Fg = sparse(obj.Fg);
        end

        function Frot = genSysFrot(obj, zti_sys)
            obj.updateRotationFrame(obj.rotationFrameStep);
            n = obj.dof_sys;
            obj.Frot = sparse(n,1);

            xti_sys   = zti_sys(1:n);
            dxti_sys  = zti_sys(n+1:2*n);
            ddxti_sys = zeros(n,1);

            obj.updateNodesFromState(obj.rotationFrameStep, xti_sys, dxti_sys, ddxti_sys);

            for e = 1:numel(obj.elements)
                elm = obj.elements{e};
                obj.setElementRotationFrame(elm);

                if obj.isRotationFrameEnabled()
                    elm.get_consistentRotFrameTangents();
                else
                    elm.Frot = zeros(elm.numDOF,1);
                end

                id = elm.dofID;
                obj.Frot(id) = obj.Frot(id) + elm.Frot;
            end
            Frot = sparse(obj.Frot);
        end

        function FrotPdimPart = genSysFrotPdimPart(obj, zti_sys, axisR, ...
                mu0, stateOrder, parameterMultiIndex)
            %GENSYSFROTPDIMPART Assemble one delta-expanded Frot coefficient.
            %
            % mu = [omega; dotOmega] = mu0 + delta. parameterMultiIndex
            % selects [ellOmega, ellDotOmega], while stateOrder selects the
            % homogeneous degree in z = [q; dq]. This path is used only by
            % semi-intrusive PDIM calculations and never reconstructs the
            % complete time-domain Frot.
            validateattributes(axisR, {'numeric'}, ...
                {'vector','numel',3,'real','finite'}, mfilename, 'axisR');
            axisR = axisR(:);
            if abs(norm(axisR) - 1) > 1e-12
                error('RotFramePDIM:AxisNotUnit', ...
                    'axisR must be a unit vector.');
            end

            validateattributes(mu0, {'numeric'}, ...
                {'vector','numel',2,'real','finite'}, mfilename, 'mu0');
            mu0 = mu0(:);
            validateattributes(stateOrder, {'numeric'}, ...
                {'scalar','integer','positive'}, mfilename, 'stateOrder');
            validateattributes(parameterMultiIndex, {'numeric'}, ...
                {'vector','numel',2,'integer','nonnegative'}, ...
                mfilename, 'parameterMultiIndex');
            parameterMultiIndex = parameterMultiIndex(:).';
            if sum(parameterMultiIndex) > 2
                error('RotFramePDIM:UnsupportedParameterOrder', ...
                    'The rotating-frame force is at most quadratic in the parameters.');
            end

            n = obj.dof_sys;
            if numel(zti_sys) < 2*n
                error('RotFramePDIM:InvalidStateSize', ...
                    'zti_sys must contain at least [q; dq].');
            end
            xti_sys = zti_sys(1:n);
            dxti_sys = zti_sys(n+1:2*n);
            ddxti_sys = zeros(n,1);
            obj.updateNodesFromState(obj.rotationFrameStep, ...
                xti_sys, dxti_sys, ddxti_sys);

            FrotPdimPart = zeros(n,1);
            for e = 1:numel(obj.elements)
                elm = obj.elements{e};
                FrotElm = elm.get_FrotPdimPart(axisR, mu0, ...
                    stateOrder, parameterMultiIndex);
                id = elm.dofID;
                FrotPdimPart(id) = FrotPdimPart(id) + FrotElm;
            end
            FrotPdimPart = sparse(FrotPdimPart);
        end

        function Fe = genSysFe(obj, zti_sys)
            n = obj.dof_sys;
            obj.Fe = sparse(n,1);
            
            xti_sys   = zti_sys(1:n);
            dxti_sys  = zti_sys(n+1:2*n);
            ddxti_sys = zeros(n,1);

            for n = 1:numel(obj.nodes)
                id = obj.nodes{n}.dofID;
                qti_node   = xti_sys(id);
                dqti_node  = dxti_sys(id);
                ddqti_node = ddxti_sys(id);
                obj.nodes{n}.updateNode(2, qti_node, dqti_node, ddqti_node);
            end
            
            for e = 1:numel(obj.elements)
                elm = obj.elements{e};

                elm.get_Fe();

                id = elm.dofID;

                obj.Fe(id) = obj.Fe(id) + elm.Fe;
            end
            Fe = sparse(obj.Fe);
        end

        function Fe1 = genSysFe1(obj, zti_sys)
            n = obj.dof_sys;
            obj.Fe1 = sparse(n,1);
            obj.nonlinOrder.Fe1 = 1;

            xti_sys   = zti_sys(1:n);
            dxti_sys  = zti_sys(n+1:2*n);
            ddxti_sys = zeros(n,1);

            for n = 1:numel(obj.nodes)
                id = obj.nodes{n}.dofID;
                qti_node   = xti_sys(id);
                dqti_node  = dxti_sys(id);
                ddqti_node = ddxti_sys(id);
                obj.nodes{n}.updateNode(2, qti_node, dqti_node, ddqti_node);
            end
            
            for e = 1:numel(obj.elements)
                elm = obj.elements{e};

                elm.get_Fe1();

                id = elm.dofID;

                obj.Fe1(id) = obj.Fe1(id) + elm.Fe1;
            end
            Fe1 = sparse(obj.Fe1);
        end

        function Fe2 = genSysFe2(obj, zti_sys)
            n = obj.dof_sys;
            obj.Fe2 = sparse(n,1);
            obj.nonlinOrder.Fe2 = 2;

            xti_sys   = zti_sys(1:n);
            dxti_sys  = zti_sys(n+1:2*n);
            ddxti_sys = zeros(n,1);

            for n = 1:numel(obj.nodes)
                id = obj.nodes{n}.dofID;
                qti_node   = xti_sys(id);
                dqti_node  = dxti_sys(id);
                ddqti_node = ddxti_sys(id);
                obj.nodes{n}.updateNode(2, qti_node, dqti_node, ddqti_node);
            end
            
            for e = 1:numel(obj.elements)
                elm = obj.elements{e};

                elm.get_Fe2();

                id = elm.dofID;

                obj.Fe2(id) = obj.Fe2(id) + elm.Fe2;
            end
            Fe2 = sparse(obj.Fe2);
        end

        function Fe3 = genSysFe3(obj, zti_sys)
            n = obj.dof_sys;
            obj.Fe3 = sparse(n,1);
            obj.nonlinOrder.Fe3 = 3;

            xti_sys   = zti_sys(1:n);
            dxti_sys  = zti_sys(n+1:2*n);
            ddxti_sys = zeros(n,1);

            for n = 1:numel(obj.nodes)
                id = obj.nodes{n}.dofID;
                qti_node   = xti_sys(id);
                dqti_node  = dxti_sys(id);
                ddqti_node = ddxti_sys(id);
                obj.nodes{n}.updateNode(2, qti_node, dqti_node, ddqti_node);
            end
            
            for e = 1:numel(obj.elements)
                elm = obj.elements{e};

                elm.get_Fe3();

                id = elm.dofID;

                obj.Fe3(id) = obj.Fe3(id) + elm.Fe3;
            end
            Fe3 = sparse(obj.Fe3);
        end

        function Fe5 = genSysFe5(obj, zti_sys)
            n = obj.dof_sys;
            obj.Fe5 = sparse(n,1);
            obj.nonlinOrder.Fe5 = 5;

            xti_sys   = zti_sys(1:n);
            dxti_sys  = zti_sys(n+1:2*n);
            ddxti_sys = zeros(n,1);

            for n = 1:numel(obj.nodes)
                id = obj.nodes{n}.dofID;
                qti_node   = xti_sys(id);
                dqti_node  = dxti_sys(id);
                ddqti_node = ddxti_sys(id);
                obj.nodes{n}.updateNode(2, qti_node, dqti_node, ddqti_node);
            end
            
            for e = 1:numel(obj.elements)
                elm = obj.elements{e};

                elm.get_Fe5();

                id = elm.dofID;

                obj.Fe5(id) = obj.Fe5(id) + elm.Fe5;
            end
            Fe5 = sparse(obj.Fe5);
        end

        function Fext = genSysFext(obj, zti_sys)
            obj.updateRotationFrame(obj.rotationFrameStep);
            n = obj.dof_sys;
            xti_sys   = zti_sys(1:n);
            dxti_sys  = zti_sys(n+1:2*n); 
            ddxti_sys = zeros(n,1);
            obj.Fext = sparse(n,1);
            obj.updateNodesFromState(obj.rotationFrameStep, xti_sys, dxti_sys, ddxti_sys);

            if strcmpi(obj.Fconstraints.kind{1}{1}, 'force boundary')
                load = obj.Fconstraints.loadList{obj.rotationFrameStep};
                for ii = 1:numel(obj.Fconstraints.dof)
                    nodeID   = obj.Fconstraints.dof{ii}{1};
                    dofID    = obj.nodes{nodeID}.dofID;
                    loadVec = load{ii};
                    if obj.isRotationFrameEnabled()
                        loadVec = [obj.rotFrameCurrent.Lambda.' * loadVec(1:3); ...
                                   obj.rotFrameCurrent.Lambda.' * loadVec(4:6)];
                    end
                    quatLoad = 2 * obj.nodes{nodeID}.H{obj.rotationFrameStep}.' * loadVec(4:6);

                    obj.Fext(dofID(1:3)) = obj.Fext(dofID(1:3)) + loadVec(1:3);
                    obj.Fext(dofID(4:7)) = obj.Fext(dofID(4:7)) + quatLoad;
                end
            end
            Fext = sparse(obj.Fext);
        end

        function C = genSysC(obj, zti_sys)
            n = obj.dof_sys;
            xti_sys   = zti_sys(1:n);
            dxti_sys  = zti_sys(n+1:2*n); 
            ddxti_sys = zeros(n,1);
            for n = 1:numel(obj.nodes)
                id = obj.nodes{n}.dofID;
                qti_node   = xti_sys(id);
                dqti_node  = dxti_sys(id);
                ddqti_node = ddxti_sys(id);
                obj.nodes{n}.updateNode(2, qti_node, dqti_node, ddqti_node);
            end

            if strcmpi(obj.Dconstraints.kind{1}{1}, 'displacement boundary')
                load = obj.Dconstraints.loadList{2};
                obj.C  = sparse(obj.dof_constraint + numel(obj.nodes), 1);
                kk = 1;
                for ii = 1:numel(obj.Dconstraints.dof)
                    nodeID  = obj.Dconstraints.dof{ii}{1};
                    dofList = obj.Dconstraints.dof{ii}{2};
                    dofID   = obj.nodes{nodeID}.dofID;
                    for jj = dofList
                        obj.C(kk) = xti_sys(dofID(jj)) - load{1}(kk);
                        kk = kk + 1; 
                    end
                end

                % special operation for gradient discontinuity
                if strcmpi(obj.Dconstraints.kind{2}{1}, 'gradient discontinuity')
                    nodeA = obj.Dconstraints.kind{2}{2}{1};
                    nodeB = obj.Dconstraints.kind{2}{2}{2};
                    e_rel = obj.Dconstraints.kind{2}{2}{3};
                    dofID_A = obj.nodes{nodeA}.dofID;
                    dofID_B = obj.nodes{nodeB}.dofID;
                    obj.C(kk:kk+2) = xti_sys(dofID_B(1:3)) - xti_sys(dofID_A(1:3));
                    e_A = xti_sys(dofID_A(4:7));
                    e_B = xti_sys(dofID_B(4:7));
                    tmp = Base.quatMultiply(e_A, e_rel);
                    e_constraint = e_B - tmp;
                    obj.C(kk+3:kk+5) = e_constraint(2:4);
                end

                for n = 1:numel(obj.nodes)
                    obj.C(obj.dof_constraint+n) = sum(xti_sys(7 * n - 3 : 7 * n).^2) - 1;
                end
            end
            C  = sparse(obj.C);
        end

        function Cq = genSysCq(obj, zti_sys)
            n = obj.dof_sys;
            xti_sys   = zti_sys(1:n);
            dxti_sys  = zti_sys(n+1:2*n); 
            ddxti_sys = zeros(n,1);
            for n = 1:numel(obj.nodes)
                id = obj.nodes{n}.dofID;
                qti_node   = xti_sys(id);
                dqti_node  = dxti_sys(id);
                ddqti_node = ddxti_sys(id);
                obj.nodes{n}.updateNode(2, qti_node, dqti_node, ddqti_node);
            end

            if strcmpi(obj.Dconstraints.kind{1}{1}, 'displacement boundary')
                load = obj.Dconstraints.loadList{2};
                obj.Cq = sparse(obj.dof_constraint + numel(obj.nodes), obj.dof_sys);
                kk = 1;
                for ii = 1:numel(obj.Dconstraints.dof)
                    dofList = obj.Dconstraints.dof{ii}{2};
                    for jj = dofList
                        kk = kk + 1; 
                    end
                end
                obj.Cq(1:kk-1,:) = load{2};

                % special operation for gradient discontinuity
                if strcmpi(obj.Dconstraints.kind{2}{1}, 'gradient discontinuity')
                    nodeA = obj.Dconstraints.kind{2}{2}{1};
                    nodeB = obj.Dconstraints.kind{2}{2}{2};
                    e_rel = obj.Dconstraints.kind{2}{2}{3};
                    R_rel = Base.genQuatRightMultiMatrix(e_rel);
                    dofID_A = obj.nodes{nodeA}.dofID;
                    dofID_B = obj.nodes{nodeB}.dofID;

                    obj.Cq(kk:kk+2,dofID_A(1:3))   =-eye(3);
                    obj.Cq(kk:kk+2,dofID_B(1:3))   = eye(3);
                    obj.Cq(kk+3:kk+5,dofID_A(4:7)) =-R_rel(2:4,1:4);
                    obj.Cq(kk+3:kk+5,dofID_B(5:7)) = eye(3);
                end

                for n = 1:numel(obj.nodes)
                    obj.Cq(obj.dof_constraint+n,obj.nodes{n}.dofID(4:7)) = 2 * xti_sys(obj.nodes{n}.dofID(4:7));
                end
            end
            Cq = sparse(obj.Cq);
        end

        function Kb = genSysKb(obj, zti_sys, lmti)
            n = obj.dof_sys;
            xti_sys   = zti_sys(1:n);
            dxti_sys  = zti_sys(n+1:2*n); 
            ddxti_sys = zeros(n,1);

            for n = 1:numel(obj.nodes)
                id = obj.nodes{n}.dofID;
                qti_node   = xti_sys(id);
                dqti_node  = dxti_sys(id);
                ddqti_node = ddxti_sys(id);
                obj.nodes{n}.updateNode(2, qti_node, dqti_node, ddqti_node);
            end

            if strcmpi(obj.Dconstraints.kind{1}{1}, 'displacement boundary')
                load = obj.Dconstraints.loadList{2};
                obj.Kb = sparse(obj.dof_sys, obj.dof_sys);
                
                for n = 1:numel(obj.nodes)
                    obj.Kb(7*n-3:7*n, 7*n-3:7*n) = 2 * lmti(n+obj.dof_constraint) * eye(4);
                end
                obj.Kb = obj.Kb + load{3};
            end
            Kb = sparse(obj.Kb);
        end

        function Cqlmda2order = genSysCqlmda2order(obj, zti_sys)
            n = obj.dof_sys;
            xti_sys   = zti_sys(1:n);
            dxti_sys  = zti_sys(n+1:2*n); 
            lmti_sys  = zti_sys(2*n+1:end);
            ddxti_sys = zeros(n,1);
            for n = 1:numel(obj.nodes)
                id = obj.nodes{n}.dofID;
                qti_node   = xti_sys(id);
                dqti_node  = dxti_sys(id);
                ddqti_node = ddxti_sys(id);
                obj.nodes{n}.updateNode(2, qti_node, dqti_node, ddqti_node);
            end

            if strcmpi(obj.Dconstraints.kind{1}{1}, 'displacement boundary')
                load = obj.Dconstraints.loadList{2};
                Cq_tmp = zeros(obj.dof_constraint + numel(obj.nodes), obj.dof_sys);
                kk = 1;
                for ii = 1:numel(obj.Dconstraints.dof)
                    dofList = obj.Dconstraints.dof{ii}{2};
                    for jj = dofList
                        kk = kk + 1; 
                    end
                end
                Cq_tmp(1:kk-1,:) = load{2};

                % special operation for gradient discontinuity
                if strcmpi(obj.Dconstraints.kind{2}{1}, 'gradient discontinuity')
                    nodeA = obj.Dconstraints.kind{2}{2}{1};
                    nodeB = obj.Dconstraints.kind{2}{2}{2};
                    e_rel = obj.Dconstraints.kind{2}{2}{3};
                    R_rel = Base.genQuatRightMultiMatrix(e_rel);
                    dofID_A = obj.nodes{nodeA}.dofID;
                    dofID_B = obj.nodes{nodeB}.dofID;

                    Cq_tmp(kk:kk+2,dofID_A(1:3))   =-eye(3);
                    Cq_tmp(kk:kk+2,dofID_B(1:3))   = eye(3);
                    Cq_tmp(kk+3:kk+5,dofID_A(4:7)) =-R_rel(2:4,1:4);
                    Cq_tmp(kk+3:kk+5,dofID_B(5:7)) = eye(3);
                end

                for n = 1:numel(obj.nodes)
                    Cq_tmp(obj.dof_constraint+n,obj.nodes{n}.dofID(4:7)) = 2 * xti_sys(obj.nodes{n}.dofID(4:7));
                end
            end
            Cq_tmp(1:obj.dof_constraint,:) = 0;
            tmp = Cq_tmp.' * lmti_sys;
            Cqlmda2order = sparse(tmp);
        end
    end
end
