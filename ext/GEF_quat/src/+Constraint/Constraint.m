classdef Constraint < matlab.mixin.SetGetExactNames
    %CONSTRAINT : define the constraint class
   
    
    properties
        kind            % kind of constriants: force or displacement boundary
        dof             % constrained DOF: first element: element number, 
                        % second element: dof of specific node, (cell array)
        loadList        % for dispConstraint: loadList contains three parts: C, Cq, Kb
                        % for forcConstraint: loadList contains force
        dof_constraint  % number of constrained DOF
        C               % the constraint equations
        Cq              % the constraint Jacobian
        Kb              % the added stiffness matrix of constraint
        
        Jac             % Jacobian matrix for DAEs
        Res             % residual vector for DAEs
    end
    
    methods
        function obj = Constraint(dof, kind, loadList, dof_constraint)
            obj.dof      = dof;
            obj.kind     = kind;
            obj.loadList = loadList;
            obj.dof_constraint = dof_constraint;
        end
        
        function [Jac, Res] = change_K_P(obj, solver, i)
            if strcmpi(obj.kind{1}{1}, 'displacement boundary')
                load = obj.loadList{i};
                M    = solver.M;
                Damp = solver.Damp;
                Kt   = solver.Kt;
                Fe   = solver.Fe;
                Fi   = solver.Fi;
                Fg   = solver.Fg;
                Frot = solver.Frot;
                qti_sys   = solver.qti_sys;
                dqti_sys  = solver.dqti_sys;
                ddqti_sys = solver.ddqti_sys;
                lmti      = solver.lmti_sys;
                dof_sys   = solver.dof_sys;
                nodenum   = solver.nodenum;
                beta_hat  = solver.beta_hat;
                gamma_hat = solver.gamma_hat;

                % generate the C, Cq, and Kb
                obj.C  = zeros(obj.dof_constraint + nodenum, 1);
                obj.Cq = zeros(obj.dof_constraint + nodenum, dof_sys);
                obj.Kb = zeros(dof_sys, dof_sys);
                kk = 1;
                for ii = 1:numel(obj.dof)
                    nodeID  = obj.dof{ii}{1};
                    dofList = obj.dof{ii}{2};
                    dofID   = solver.nodes{nodeID}.dofID;
                    for jj = dofList
                        obj.C(kk) = qti_sys(dofID(jj)) - load{1}(kk);
                        kk = kk + 1; 
                    end
                end
                obj.Cq(1:kk-1,:) = load{2};

                % special operation for gradient discontinuity
                if strcmpi(obj.kind{2}{1}, 'gradient discontinuity')
                    nodeA = obj.kind{2}{2}{1};
                    nodeB = obj.kind{2}{2}{2};
                    e_rel = obj.kind{2}{2}{3};
                    R_rel = Base.genQuatRightMultiMatrix(e_rel);
                    dofID_A = solver.nodes{nodeA}.dofID;
                    dofID_B = solver.nodes{nodeB}.dofID;
                    obj.C(kk:kk+2) = qti_sys(dofID_B(1:3)) - qti_sys(dofID_A(1:3));
                    e_A = qti_sys(dofID_A(4:7));
                    e_B = qti_sys(dofID_B(4:7));
                    tmp = Base.quatMultiply(e_A, e_rel);
                    e_constraint = e_B - tmp;
                    obj.C(kk+3:kk+5) = e_constraint(2:4);
                    obj.Cq(kk:kk+2,dofID_A(1:3))   =-eye(3);
                    obj.Cq(kk:kk+2,dofID_B(1:3))   = eye(3);
                    obj.Cq(kk+3:kk+5,dofID_A(4:7)) =-R_rel(2:4,1:4);
                    obj.Cq(kk+3:kk+5,dofID_B(5:7)) = eye(3);
                end

                for n = 1:nodenum
                    obj.C(obj.dof_constraint+n) = sum(qti_sys(7 * n - 3 : 7 * n).^2) - 1;
                    obj.Cq(obj.dof_constraint+n,solver.nodes{n}.dofID(4:7)) = 2 * qti_sys(solver.nodes{n}.dofID(4:7));
                    obj.Kb(7*n-3:7*n, 7*n-3:7*n) = 2 * lmti(n+obj.dof_constraint) * eye(4);
                end
                obj.Kb = obj.Kb + load{3};

                if solver.isRotationFrameEnabled()
                    tangent = M * beta_hat + (Damp + solver.Ct) * gamma_hat + Kt;
                else
                    tangent = M * beta_hat + Damp * gamma_hat + Kt;
                end

                Jac = [tangent   obj.Cq.'  % + obj.Kb
                       obj.Cq                                 zeros(size(obj.Cq,1))];
                Res = [M * ddqti_sys + Damp * dqti_sys + Fe - Fi + Frot - Fg + obj.Cq.' * lmti
                       obj.C];
                obj.Jac = Jac;
                obj.Res = Res;

            elseif strcmpi(obj.kind{1}{1}, 'force boundary')
                load = obj.loadList{i};
                Fext = zeros(solver.dof_sys + solver.dof_constraint + solver.nodenum,1);
                for ii = 1:numel(obj.dof)
                    nodeID   = obj.dof{ii}{1};
                    % dofList  = obj.dof{ii}{2};
                    dofID    = solver.nodes{nodeID}.dofID;
                    loadVec = load{ii};
                    if solver.isRotationFrameEnabled()
                        LambdaR = solver.rotFrameCurrent.Lambda;
                        loadVec = [LambdaR.' * loadVec(1:3); LambdaR.' * loadVec(4:6)];
                    end
                    quatLoad = 2 * solver.nodes{nodeID}.H{i}.' * loadVec(4:6);
                    
                    Fext(dofID(1:3)) = Fext(dofID(1:3)) + loadVec(1:3);
                    Fext(dofID(4:7)) = Fext(dofID(4:7)) + quatLoad;
                    
                end
                Jac = solver.Jac;
                Res = solver.Res - Fext;
                obj.Jac = Jac;
                obj.Res = Res;
            else
                error('Constraint class only support displacement or force boundary')
            end
        end
    end
end

