classdef HM3DCurvRodL3_quat < matlab.mixin.SetGetExactNames
    %HM3DCurvRodL3_quat : define the beam element with three nodes
    
    properties
        %------element/node information------------------------------------
        numNodes = 3    % number of nodes per element
        numDOFperN      % number of DOF per node
        numDOF          % number of DOF per element
        nodes           % node class 
        dofID           % dof ID list of the element
        
        %------material/cross-section properties---------------------------
        E               % Young's modulus
        Gs              % shear modulus 
        Gm = 4*pi*1e-7  % magnetic permeability
        density         % density
        Area            % area of cross section
        ks = [1; 1]     % shear correction factor, [2 x 1] vector
        I1              % polar moment of inertia
        I2              % moment of inertia at 2-direction
        I3              % moment of inertia at 3-direction
        RM = [0; 0; 0]  % residual magnetic flux density
        KI = 1          % correction factor for mass and inertia moment
        rhoA            % linear density
        CA_mB           % corresponding to stretch-shear under material description
        CI_mB           % corresponding to bending under material description
        Jrho            % rho * diag([I1, I2, I3])
        sectionInputMode = 'standard'

        E_preIntp               
        G_preIntp                  
        density_preIntp          
        A_preIntp                
        ks1_preIntp     % put at the first position of k
        ks2_preIntp     % put at the second position of k
        I1_preIntp  
        I2_preIntp               
        I3_preIntp                    
        RM1_preIntp
        RM2_preIntp
        RM3_preIntp
        KI_preIntp              
        rhoA_node
        Jrho_node
        CA_mB_node
        CI_mB_node

        %------gauss quadrature information--------------------------------
        numGauss        % number of gauss points

        %------state quantities--------------------------------------------
        q_elm           % elementary generalized coordinates
        q0_elm          % initial generalized coordinates of element
        dq_elm

        %------element matrices and vectors--------------------------------
        M               % element mass matrix
        M0              % constant part of element mass matrix
        Mnl             % nonlinear part of element mass matrix
        Ct              % element gyro matrix
        Kt              % element tangent stiffness matrix
        Fi              % element force vector corresponding to inertia
        Fg              % element body force vector
        Fe              % element elastic force vector
        Fe1             % 1st order part of element elastic force vector
        Fe2             % 2nd order part of element elastic force vector
        Fe3             % 3rd order part of element elastic force vector
        Fe5             % 5th order part of element elastic force vector
        dM_dq           % third tensor

        %------rotating frame quantities-----------------------------------
        rotationFrame = struct('enabled', false)
        Frot            % extra inertia force caused by the rotating frame
        Crot            % approximate damping-like tangent from rotating frame
        Krot            % approximate stiffness-like tangent from rotating frame
        Gr1             % translational gyroscopic matrix from frame rotation
        Gr2             % translational centrifugal matrix from frame rotation
        dGr1            % time derivative of Gr1
    end
    
    methods
        function obj = HM3DCurvRodL3_quat(nodes, reduce)
            obj.nodes      = nodes;
            obj.numDOFperN = nodes{1}.numDOF;
            obj.numDOF     = obj.numNodes * obj.numDOFperN;
            
            obj.dofID  = [];
            obj.q0_elm = [];

            %---extract material/cross-section properties from nodes---
            obj.E_preIntp       = [];              
            obj.G_preIntp       = [];             
            obj.density_preIntp = [];           
            obj.A_preIntp       = [];           
            obj.ks1_preIntp     = [];  
            obj.ks2_preIntp     = [];
            obj.I1_preIntp      = [];             
            obj.I2_preIntp      = [];           
            obj.I3_preIntp      = [];           
            obj.RM1_preIntp     = [];          
            obj.KI_preIntp      = [];
            obj.rhoA_node       = zeros(obj.numNodes,1);
            obj.Jrho_node       = zeros(3,3,obj.numNodes);
            obj.CA_mB_node      = zeros(3,3,obj.numNodes);
            obj.CI_mB_node      = zeros(3,3,obj.numNodes);

            sectionModes = strings(obj.numNodes,1);

            for n = 1:obj.numNodes
                obj.dofID  = [obj.dofID;  nodes{n}.dofID(:)];
                obj.q0_elm = [obj.q0_elm; nodes{n}.q(:,1)];

                if isa(nodes{n}.crossSection, 'CrossSection.DirectRod3D')
                    sectionModes(n) = "direct";
                    obj.rhoA_node(n) = nodes{n}.crossSection.rhoA;
                    obj.Jrho_node(:,:,n) = nodes{n}.crossSection.Jrho;
                    obj.CA_mB_node(:,:,n) = nodes{n}.crossSection.CA_mB;
                    obj.CI_mB_node(:,:,n) = nodes{n}.crossSection.CI_mB;
                else
                    sectionModes(n) = "standard";
                    obj.E_preIntp       = [obj.E_preIntp; nodes{n}.Material.E * ones(obj.numDOFperN,1)];              
                    obj.G_preIntp       = [obj.G_preIntp; nodes{n}.Material.G * ones(obj.numDOFperN,1)];             
                    obj.density_preIntp = [obj.density_preIntp; nodes{n}.Material.density * ones(obj.numDOFperN,1)];           
                    obj.A_preIntp       = [obj.A_preIntp; nodes{n}.crossSection.A * ones(obj.numDOFperN,1)];           
                    obj.ks1_preIntp     = [obj.ks1_preIntp; nodes{n}.crossSection.k(1) * ones(obj.numDOFperN,1)];  
                    obj.ks2_preIntp     = [obj.ks2_preIntp; nodes{n}.crossSection.k(2) * ones(obj.numDOFperN,1)]; 
                    obj.I1_preIntp      = [obj.I1_preIntp; nodes{n}.crossSection.I1 * ones(obj.numDOFperN,1)];                 
                    obj.I2_preIntp      = [obj.I2_preIntp; nodes{n}.crossSection.I2 * ones(obj.numDOFperN,1)];           
                    obj.I3_preIntp      = [obj.I3_preIntp; nodes{n}.crossSection.I3 * ones(obj.numDOFperN,1)]; 
                end
                obj.RM1_preIntp     = [obj.RM1_preIntp; nodes{n}.RM(1) * ones(obj.numDOFperN,1)];
                obj.RM2_preIntp     = [obj.RM2_preIntp; nodes{n}.RM(2) * ones(obj.numDOFperN,1)]; 
                obj.RM3_preIntp     = [obj.RM3_preIntp; nodes{n}.RM(3) * ones(obj.numDOFperN,1)]; 
                obj.KI_preIntp      = [obj.KI_preIntp; nodes{n}.KI * ones(obj.numDOFperN,1)];
            end

            if any(sectionModes ~= sectionModes(1))
                error('Element:MixedSectionInput', ...
                    'A single element cannot mix Rod3D and DirectRod3D nodes.');
            end
            obj.sectionInputMode = char(sectionModes(1));

            % choose gauss quadrature scheme
            if reduce == 1
                obj.numGauss = 2; % Reduced
            elseif reduce == 0
                obj.numGauss = 3; % Full (Standard)
            else
                obj.numGauss = 4; % Over integration
            end

            obj.updateSectionProperties(0);
        end
        
        function get_K_P(obj, consistentTangent)
            if nargin < 2
                consistentTangent = false;
            end
            nDim   = obj.numDOF;
            obj.M  = zeros(nDim,nDim);
            obj.Ct = zeros(nDim,nDim);
            obj.Kt = zeros(nDim,nDim);
            obj.Fi = zeros(nDim,1);
            obj.Fg = zeros(nDim,1);
            obj.Fe = zeros(nDim,1);
            obj.Frot = zeros(nDim,1);
            obj.Crot = zeros(nDim,nDim);
            obj.Krot = zeros(nDim,nDim);
            obj.Gr1  = zeros(nDim,nDim);
            obj.Gr2  = zeros(nDim,nDim);
            obj.dGr1 = zeros(nDim,nDim);

            obj.q_elm  = [];
            obj.dq_elm = [];

            for n = 1:obj.numNodes
                obj.q_elm  = [obj.q_elm;  obj.nodes{n}.qti ];
                obj.dq_elm = [obj.dq_elm; obj.nodes{n}.dqti];     
            end

            [weights, points] = Base.gaussQuadrature(obj.numGauss);

            for g = 1 : size(points,1)
                xi = points(g);
                obj.updateSectionProperties(xi);
                [S, N, dS_dxi, dN_dxi] = Base.genShapeF(xi);
                if consistentTangent
                    % The shooting force/tangent callbacks use the same
                    % reference-domain integration measure as get_Fe,
                    % get_Fi, get_Kt and get_Ct.
                    Jac = dS_dxi * obj.q0_elm;
                else
                    Jac = dS_dxi * obj.q_elm;
                end
                Jac_det = norm(Jac);

                S_prim = dS_dxi / Jac_det;
                N_prim = dN_dxi / Jac_det;
                e  = N * obj.q_elm;
                e0 = N * obj.q0_elm;
                r_prim  = S_prim * obj.q_elm;
                r0_prim = S_prim * obj.q0_elm;
                e_prim  = N_prim * obj.q_elm;
                e0_prim = N_prim * obj.q0_elm;
                G       = Base.genG(e);    
                lambd   = Base.genA(e);
                G0      = Base.genG(e0);
                lambd0  = Base.genA(e0);

                c = Base.genc(lambd0, r0_prim);

                gam  = lambd.'  * r_prim;
                gam0 = lambd0.' * r0_prim;
                kap  = 2 * G  * e_prim;
                kap0 = 2 * G0 * e0_prim; 
                M_int  = obj.rhoA * (S.' * S) + 4 * N.' * G.' * obj.Jrho * G * N;
                Fi_int = Base.genElmFi(obj.dq_elm, N, obj.Jrho, G, obj.numDOF);
                Ct_int = Base.genElmCt(obj.dq_elm, N, obj.Jrho, G, obj.numDOF);
                Fe_int = Base.genElmFe(c, e, lambd, G, S_prim, r_prim, N, e_prim, ...
                            N_prim, obj.CA_mB, gam, gam0, obj.CI_mB, kap, kap0);
                if consistentTangent
                    Kt_int = Base.genElmExactKt(c, e, G, r_prim, N, e_prim, ...
                                lambd, S_prim, N_prim, obj.CA_mB, obj.CI_mB, gam, gam0,...
                                kap, kap0, obj.numDOF);
                else
                    Kt_int = Base.genElmKt(c, e, G, r_prim, N, e_prim, ...
                                lambd, S_prim, N_prim, obj.CA_mB, obj.CI_mB, gam, gam0,...
                                kap, kap0, obj.numDOF);
                end

                obj.M  = obj.M  + M_int  * weights(g) * Jac_det;
                obj.Ct = obj.Ct + Ct_int * weights(g) * Jac_det;
                obj.Kt = obj.Kt + Kt_int * weights(g) * Jac_det;
                obj.Fi = obj.Fi + Fi_int * weights(g) * Jac_det;
                obj.Fg = obj.Fg + obj.genBodyForce(S) * weights(g) * Jac_det;
                obj.Fe = obj.Fe + Fe_int * weights(g) * Jac_det;

                if obj.isRotationFrameEnabled()
                    [Frot_int, Crot_int, Krot_int, Gr1_int, Gr2_int, dGr1_int] = ...
                        obj.genRotFrameTerms(S, N, e, weights(g) * Jac_det, ...
                        true);
                    obj.Frot = obj.Frot + Frot_int;
                    obj.Crot = obj.Crot + Crot_int;
                    obj.Krot = obj.Krot + Krot_int;
                    obj.Gr1  = obj.Gr1  + Gr1_int;
                    obj.Gr2  = obj.Gr2  + Gr2_int;
                    obj.dGr1 = obj.dGr1 + dGr1_int;
                end
            end

            if obj.isRotationFrameEnabled()
                obj.Ct = obj.Ct + obj.Crot;
                obj.Kt = obj.Kt + obj.Krot;
            end
        end

        function updateSectionProperties(obj, xi)
            [S, ~, ~, ~] = Base.genShapeF(xi);

            obj.RM(1) = S(1,:) * obj.RM1_preIntp;
            obj.RM(2) = S(1,:) * obj.RM2_preIntp;
            obj.RM(3) = S(1,:) * obj.RM3_preIntp;
            obj.KI = S(1,:) * obj.KI_preIntp;

            if strcmpi(obj.sectionInputMode, 'direct')
                L = [S(1,1), S(1,obj.numDOFperN + 1), ...
                     S(1,2 * obj.numDOFperN + 1)];
                obj.rhoA = L * obj.rhoA_node;
                obj.Jrho = zeros(3,3);
                obj.CA_mB = zeros(3,3);
                obj.CI_mB = zeros(3,3);
                for n = 1:obj.numNodes
                    obj.Jrho = obj.Jrho + L(n) * obj.Jrho_node(:,:,n);
                    obj.CA_mB = obj.CA_mB + L(n) * obj.CA_mB_node(:,:,n);
                    obj.CI_mB = obj.CI_mB + L(n) * obj.CI_mB_node(:,:,n);
                end
                return
            end

            obj.E = S(1,:) * obj.E_preIntp;
            obj.Gs = S(1,:) * obj.G_preIntp;
            obj.density = S(1,:) * obj.density_preIntp;
            obj.Area = S(1,:) * obj.A_preIntp;
            obj.ks(1) = S(1,:) * obj.ks1_preIntp;
            obj.ks(2) = S(1,:) * obj.ks2_preIntp;
            obj.I1 = S(1,:) * obj.I1_preIntp;
            obj.I2 = S(1,:) * obj.I2_preIntp;
            obj.I3 = S(1,:) * obj.I3_preIntp;
            obj.rhoA = obj.Area * obj.density;
            obj.CA_mB = diag([obj.E * obj.Area; ...
                              obj.Gs * obj.Area * obj.ks(1); ...
                              obj.Gs * obj.Area * obj.ks(2)]);
            obj.CI_mB = diag([obj.Gs * obj.I1; ...
                              obj.E * obj.I2; ...
                              obj.E * obj.I3]);
            obj.Jrho = obj.density * diag([obj.I1; obj.I2; obj.I3]) * obj.KI;
        end

        function setRotationFrame(obj, rotationFrame)
            obj.rotationFrame = rotationFrame;
        end

        function Fg = genBodyForce(obj, S)
            g_gp = zeros(3,1);
            for n = 1:obj.numNodes
                g_node = obj.nodes{n}.g;
                if iscell(g_node)
                    g_node = g_node{1};
                end
                if obj.isRotationFrameEnabled()
                    g_node = obj.rotationFrame.Lambda.' * g_node(:);
                end
                g_gp = g_gp + S(:, obj.numDOFperN * (n - 1) + (1:3)) * g_node(:);
            end
            Fg = obj.rhoA * S.' * g_gp;
        end

        function tf = isRotationFrameEnabled(obj)
            tf = isstruct(obj.rotationFrame) && isfield(obj.rotationFrame, 'enabled') ...
                && obj.rotationFrame.enabled;
        end

        function [Frot, Crot, Krot, Gr1, Gr2, dGr1] = genRotFrameTerms(obj, S, N, epsilon, weightJac, consistentTangent)
            if nargin < 6
                consistentTangent = false;
            end
            nDim = obj.numDOF;
            % Frot = zeros(nDim,1);
            % Crot = zeros(nDim,nDim);
            % Krot = zeros(nDim,nDim);
            % Gr1  = zeros(nDim,nDim);
            % Gr2  = zeros(nDim,nDim);
            % dGr1 = zeros(nDim,nDim);

            Omega    = obj.rotationFrame.Omega(:);
            dotOmega = obj.rotationFrame.dotOmega(:);

            OmegaHat    = Base.genSkew(Omega);
            dotOmegaHat = Base.genSkew(dotOmega);

            dotEpsilon = N * obj.dq_elm;
            G          = Base.genG(epsilon);
            dotG       = Base.genG(dotEpsilon);
            Lambda     = Base.genA(epsilon);
            dA_de = Base.genElmdA_de(epsilon);
            dotA = sum(dA_de.*reshape(dotEpsilon,1,1,4),3);

            Gr1  = 2 * obj.rhoA * S.' * OmegaHat * S * weightJac;
            Gr2  = obj.rhoA * S.' * (OmegaHat.' * OmegaHat) * S * weightJac;
            dGr1 = 2 * obj.rhoA * S.' * dotOmegaHat * S * weightJac;

            a = dotG.' * obj.Jrho * Lambda.' * Omega ...
              + G.' * obj.Jrho * dotA.' * Omega ...
              + G.' * obj.Jrho * Lambda.' * dotOmega;
            Fi1 = 2 * N.' * a * weightJac;

            Fi2 = 0.5 * Gr1.' * obj.dq_elm + Gr2 * obj.q_elm;
            dG_de = Base.genElmdG_de();
            for j = 1:nDim
                gj = N(:,j);
                G_j = zeros(3,4);
                A_j = zeros(3,3);
                for alpha = 1:4
                    G_j = G_j + dG_de(:,:,alpha) * gj(alpha);
                    A_j = A_j + dA_de(:,:,alpha) * gj(alpha);
                end

                Fi2(j) = Fi2(j) + 0.5 * ( ...
                    4 * dotEpsilon.' * G_j.' * obj.Jrho * Lambda.' * Omega + ...
                    4 * dotEpsilon.' * G.' * obj.Jrho * A_j.' * Omega + ...
                    Omega.' * A_j * obj.Jrho * Lambda.' * Omega + ...
                    Omega.' * Lambda * obj.Jrho * A_j.' * Omega) * weightJac;
            end

            Frot = 0.5 * (dGr1 * obj.q_elm + Gr1 * obj.dq_elm) + Fi1 - Fi2;

            if consistentTangent
                [dFi1_dq, dFi1_ddq, dFi2_dq, dFi2_ddq] = ...
                    obj.genRotFrameForceJacobians(S, N, epsilon, Gr1, Gr2, weightJac);
                Crot = 0.5 * Gr1 + dFi1_ddq - dFi2_ddq;
                Krot = 0.5 * dGr1 + dFi1_dq - dFi2_dq;
            else
                % Approximate effective tangent used by the time integrator.
                Crot = 0.5 * Gr1;
                Krot = 0.5 * dGr1 - Gr2;
            end
        end

        function [dFi1_dq, dFi1_ddq, dFi2_dq, dFi2_ddq] = genRotFrameForceJacobians(obj, ~, N, epsilon, Gr1, Gr2, weightJac)
            nDim = obj.numDOF;
            dFi1_dq  = zeros(nDim,nDim);
            dFi1_ddq = zeros(nDim,nDim);
            dFi2_dq  = Gr2;
            dFi2_ddq = 0.5 * Gr1.';

            Omega    = obj.rotationFrame.Omega(:);
            dotOmega = obj.rotationFrame.dotOmega(:);

            dotEpsilon = N * obj.dq_elm;
            G          = Base.genG(epsilon);
            dotG       = Base.genG(dotEpsilon);
            Lambda     = Base.genA(epsilon);

            dG_de   = Base.genElmdG_de();
            dA_de   = Base.genElmdA_de(epsilon);
            d2A_de2 = Base.genElmd2A_de2();

            % N has nonzero columns only at the quaternion DOFs. Work on
            % that 12-DOF block instead of evaluating 21-by-21 pairs whose
            % rotational corrections are identically zero.
            activeDofs = find(any(N~=0,1));
            Nq = N(:,activeDofs);
            nActive = numel(activeDofs);

            % Batch contractions G_j = dG/de:g_j and A_j = dA/de:g_j.
            % The third dimension enumerates the active element DOFs.
            Gcols = reshape(reshape(dG_de,12,4)*Nq,3,4,nActive);
            Acols = reshape(reshape(dA_de,9,4)*Nq,3,3,nActive);

            % Frequently used page contractions. U(:,j)=A_j'*Omega,
            % Ud(:,j)=A_j'*dotOmega and X(:,j)=G_j*dotEpsilon.
            U = reshape(sum(Acols.*reshape(Omega,3,1,1),1),3,nActive);
            Ud = reshape(sum(Acols.*reshape(dotOmega,3,1,1),1),3,nActive);
            X = reshape(sum(Gcols.*reshape(dotEpsilon,1,4,1),2),3,nActive);
            Gt = @(z) reshape(sum(Gcols.*reshape(z,3,1,1),1),4,nActive);

            J = obj.Jrho;
            lambdaOmega = Lambda.'*Omega;
            lambdaDotOmega = Lambda.'*dotOmega;
            dotA = sum(dA_de.*reshape(dotEpsilon,1,1,4),3);
            dotU = zeros(3,nActive);
            for alpha = 1:4
                for beta = 1:4
                    dotU = dotU + (d2A_de2(:,:,alpha,beta).'*Omega) ...
                        * (dotEpsilon(alpha)*Nq(beta,:));
                end
            end
            da_dq = dotG.'*J*U + Gt(J*dotA.'*Omega) + G.'*J*dotU ...
                +Gt(J*lambdaDotOmega) ...
                +G.'*J*Ud;
            da_ddq = Gt(J*lambdaOmega) + G.'*J*U;

            dFi1_dq(:,activeDofs) = 2*N.'*da_dq*weightJac;
            dFi1_ddq(:,activeDofs) = 2*N.'*da_ddq*weightJac;

            % Velocity derivative of Fi2. For row i and column j the
            % correction is 2*b_i'*g_j, hence one matrix product produces
            % the complete active block.
            Y = J*U;
            B = Gt(J*lambdaOmega)+G.'*Y;
            dFi2_ddq(activeDofs,activeDofs) = ...
                dFi2_ddq(activeDofs,activeDofs) ...
                +2*(B.'*Nq)*weightJac;

            % All A_ij pages are a single contraction of the constant
            % d2A/de2 tensor with every g_i kron g_j pair. kron(Nq,Nq)
            % orders the columns as i+(j-1)*nActive, matching MATLAB's
            % column-major reshape into an nActive-by-nActive block.
            pairWeights = kron(Nq,Nq);
            AijFlat = reshape(d2A_de2,9,16)*pairWeights;
            aijRight = 4*J.'*G*dotEpsilon ...
                +J*lambdaOmega+J.'*lambdaOmega;
            aijCoefficient = Omega*aijRight.';
            aijTerms = reshape(aijCoefficient(:).'*AijFlat, ...
                nActive,nActive);

            fi2DisplacementBlock = 4*(X.'*Y+Y.'*X) ...
                +aijTerms+U.'*(J+J.')*U;
            dFi2_dq(activeDofs,activeDofs) = ...
                dFi2_dq(activeDofs,activeDofs) ...
                +0.5*fi2DisplacementBlock*weightJac;
        end

        function get_consistentRotFrameTangents(obj)
            nDim   = obj.numDOF;
            obj.Frot = zeros(nDim,1);
            obj.Crot = zeros(nDim,nDim);
            obj.Krot = zeros(nDim,nDim);
            obj.Gr1  = zeros(nDim,nDim);
            obj.Gr2  = zeros(nDim,nDim);
            obj.dGr1 = zeros(nDim,nDim);

            obj.q_elm  = [];
            obj.dq_elm = [];

            for n = 1:obj.numNodes
                obj.q_elm  = [obj.q_elm;  obj.nodes{n}.qti ];
                obj.dq_elm = [obj.dq_elm; obj.nodes{n}.dqti];
            end

            [weights, points] = Base.gaussQuadrature(obj.numGauss);

            for g = 1 : size(points,1)
                xi = points(g);
                obj.updateSectionProperties(xi);
                [S, N, dS_dxi, ~] = Base.genShapeF(xi);
                Jac = dS_dxi * obj.q0_elm;
                Jac_det = norm(Jac);
                epsilon = N * obj.q_elm;

                [Frot_int, Crot_int, Krot_int, Gr1_int, Gr2_int, dGr1_int] = ...
                    obj.genRotFrameTerms(S, N, epsilon, weights(g) * Jac_det, true);

                obj.Frot = obj.Frot + Frot_int;
                obj.Crot = obj.Crot + Crot_int;
                obj.Krot = obj.Krot + Krot_int;
                obj.Gr1  = obj.Gr1  + Gr1_int;
                obj.Gr2  = obj.Gr2  + Gr2_int;
                obj.dGr1 = obj.dGr1 + dGr1_int;
            end
        end

        function KtInertia = genClassicInertiaKt(obj)
            nDim = obj.numDOF;
            KtInertia = zeros(nDim,nDim);

            obj.q_elm  = [];
            obj.dq_elm = [];

            for n = 1:obj.numNodes
                obj.q_elm  = [obj.q_elm;  obj.nodes{n}.qti ];
                obj.dq_elm = [obj.dq_elm; obj.nodes{n}.dqti];
            end

            [weights, points] = Base.gaussQuadrature(obj.numGauss);
            dG_de = Base.genElmdG_de();

            for g = 1 : size(points,1)
                xi = points(g);
                obj.updateSectionProperties(xi);
                [~, N, dS_dxi, ~] = Base.genShapeF(xi);
                Jac = dS_dxi * obj.q0_elm;
                Jac_det = norm(Jac);

                G_col = cell(nDim,1);
                for j = 1:nDim
                    gj = N(:,j);
                    G_j = zeros(3,4);
                    for alpha = 1:4
                        G_j = G_j + dG_de(:,:,alpha) * gj(alpha);
                    end
                    G_col{j} = G_j;
                end

                for i = 1:nDim
                    G_i = G_col{i};
                    for j = 1:nDim
                        G_j = G_col{j};
                        M_ij = 4 * N.' * (G_i.' * obj.Jrho * G_j + ...
                                          G_j.' * obj.Jrho * G_i) * N;
                        KtInertia(i,j) = KtInertia(i,j) ...
                            - 0.5 * obj.dq_elm.' * M_ij * obj.dq_elm ...
                            * weights(g) * Jac_det;
                    end
                end
            end
        end

        function get_separate_nonlin(obj)
            nDim   = obj.numDOF;
            obj.M0 = zeros(nDim,nDim);
            obj.Mnl= zeros(nDim,nDim);
            obj.Ct = zeros(nDim,nDim);
            obj.Kt = zeros(nDim,nDim);
            obj.Fi = zeros(nDim,1);
            obj.Fe = zeros(nDim,1);
            obj.Fe1= zeros(nDim,1);
            obj.Fe2= zeros(nDim,1);
            obj.Fe3= zeros(nDim,1);
            obj.Fe5= zeros(nDim,1);

            obj.q_elm  = [];
            obj.dq_elm = [];

            for n = 1:obj.numNodes
                obj.q_elm  = [obj.q_elm;  obj.nodes{n}.qti ];
                obj.dq_elm = [obj.dq_elm; obj.nodes{n}.dqti];     
            end

            [weights, points] = Base.gaussQuadrature(obj.numGauss);

            for g = 1 : size(points,1)
                xi = points(g);
                obj.updateSectionProperties(xi);
                [S, N, dS_dxi, dN_dxi] = Base.genShapeF(xi);
                Jac = dS_dxi * obj.q0_elm;
                Jac_det = norm(Jac);

                S_prim = dS_dxi / Jac_det;
                N_prim = dN_dxi / Jac_det;
                e  = N * obj.q_elm;
                e0 = N * obj.q0_elm;
                r_prim  = S_prim * obj.q_elm;
                r0_prim = S_prim * obj.q0_elm;
                e_prim  = N_prim * obj.q_elm;
                e0_prim = N_prim * obj.q0_elm;
                G      = Base.genG(e);    
                lambd  = Base.genA(e);
                G0     = Base.genG(e0);
                lambd0 = Base.genA(e0);

                c = Base.genc(lambd0, r0_prim);

                gam  = lambd.'  * r_prim;
                gam0 = lambd0.' * r0_prim;
                kap  = 2 * G  * e_prim;
                kap0 = 2 * G0 * e0_prim; 
                M0_int = obj.rhoA * (S.' * S);
                Mnl_int= 4 * N.' * G.' * obj.Jrho * G * N;
                Ct_int = Base.genElmCt(obj.dq_elm, N, obj.Jrho, G, obj.numDOF);
                Kt_int = Base.genElmKt(c, e, G, r_prim, N, e_prim, ...
                            lambd, S_prim, N_prim, obj.CA_mB, obj.CI_mB, gam, gam0,...
                            kap, kap0, obj.numDOF);
                Fi_int = Base.genElmFi(obj.dq_elm, N, obj.Jrho, G, obj.numDOF);
                Fe_int = Base.genElmFe(c, e, lambd, G, S_prim, r_prim, N, e_prim, ...
                            N_prim, obj.CA_mB, gam, gam0, obj.CI_mB, kap, kap0);
                Fe1_int= Base.genElmFe1(c, e, lambd, G, S_prim, r_prim, N, e_prim, ...
                            N_prim, obj.CI_mB, kap0);
                Fe2_int= Base.genElmFe2(c, e, lambd, G, S_prim, r_prim, N, e_prim, ...
                            N_prim, obj.CA_mB, gam0);
                Fe3_int= Base.genElmFe3(c, e, lambd, G, S_prim, r_prim, N, e_prim, ...
                            N_prim, obj.CI_mB, kap);
                Fe5_int= Base.genElmFe5(c, e, lambd, G, S_prim, r_prim, N, e_prim, ...
                            N_prim, obj.CA_mB, gam);

                obj.M0 = obj.M0 + M0_int * weights(g) * Jac_det;
                obj.Mnl= obj.Mnl+ Mnl_int* weights(g) * Jac_det;
                obj.Ct = obj.Ct + Ct_int * weights(g) * Jac_det;
                obj.Kt = obj.Kt + Kt_int * weights(g) * Jac_det;
                obj.Fi = obj.Fi + Fi_int * weights(g) * Jac_det;
                obj.Fe = obj.Fe + Fe_int * weights(g) * Jac_det;
                obj.Fe1= obj.Fe1+ Fe1_int* weights(g) * Jac_det;
                obj.Fe2= obj.Fe2+ Fe2_int* weights(g) * Jac_det;
                obj.Fe3= obj.Fe3+ Fe3_int* weights(g) * Jac_det;
                obj.Fe5= obj.Fe5+ Fe5_int* weights(g) * Jac_det;
            end
        end

        function get_M0(obj)
            nDim   = obj.numDOF;
            obj.M0 = zeros(nDim,nDim);

            obj.q_elm  = [];
            obj.dq_elm = [];

            for n = 1:obj.numNodes
                obj.q_elm  = [obj.q_elm;  obj.nodes{n}.qti ];
                obj.dq_elm = [obj.dq_elm; obj.nodes{n}.dqti];     
            end

            [weights, points] = Base.gaussQuadrature(obj.numGauss);

            for g = 1 : size(points,1)
                xi = points(g);
                obj.updateSectionProperties(xi);
                [S, ~, dS_dxi, ~] = Base.genShapeF(xi);
                Jac = dS_dxi * obj.q0_elm;
                Jac_det = norm(Jac);

                M0_int = obj.rhoA * (S.' * S);
                
                obj.M0 = obj.M0 + M0_int * weights(g) * Jac_det;
            end
        end

        function get_Mnl(obj)
            nDim   = obj.numDOF;
            obj.Mnl= zeros(nDim,nDim);

            obj.q_elm  = [];
            obj.dq_elm = [];

            for n = 1:obj.numNodes
                obj.q_elm  = [obj.q_elm;  obj.nodes{n}.qti ];
                obj.dq_elm = [obj.dq_elm; obj.nodes{n}.dqti];     
            end

            [weights, points] = Base.gaussQuadrature(obj.numGauss);

            for g = 1 : size(points,1)
                xi = points(g);
                obj.updateSectionProperties(xi);
                [~, N, dS_dxi, ~] = Base.genShapeF(xi);
                Jac = dS_dxi * obj.q0_elm;
                Jac_det = norm(Jac);

                e  = N * obj.q_elm;
                G      = Base.genG(e);    

                Mnl_int= 4 * N.' * G.' * obj.Jrho * G * N;

                obj.Mnl= obj.Mnl+ Mnl_int* weights(g) * Jac_det;
            end
        end

        function get_dMdq(obj)
            nDim      = obj.numDOF;
            obj.dM_dq = zeros(nDim,nDim,nDim);

            obj.q_elm  = [];
            obj.dq_elm = [];

            for n = 1:obj.numNodes
                obj.q_elm  = [obj.q_elm;  obj.nodes{n}.qti ];
                obj.dq_elm = [obj.dq_elm; obj.nodes{n}.dqti];     
            end

            [weights, points] = Base.gaussQuadrature(obj.numGauss);

            for g = 1 : size(points,1)
                xi = points(g);
                obj.updateSectionProperties(xi);
                [~, N, dS_dxi, ~] = Base.genShapeF(xi);
                Jac = dS_dxi * obj.q0_elm;
                Jac_det = norm(Jac);

                e  = N * obj.q_elm;
                G      = Base.genG(e);    

                dM_dq_int = Base.genElmdM_dq(N, obj.Jrho, G, obj.numDOF);

                obj.dM_dq = obj.dM_dq + dM_dq_int * weights(g) * Jac_det;
            end
        end

        function get_Ct(obj)
            nDim   = obj.numDOF;
            obj.Ct = zeros(nDim,nDim);

            obj.q_elm  = [];
            obj.dq_elm = [];

            for n = 1:obj.numNodes
                obj.q_elm  = [obj.q_elm;  obj.nodes{n}.qti ];
                obj.dq_elm = [obj.dq_elm; obj.nodes{n}.dqti];     
            end

            [weights, points] = Base.gaussQuadrature(obj.numGauss);

            for g = 1 : size(points,1)
                xi = points(g);
                obj.updateSectionProperties(xi);
                [~, N, dS_dxi, ~] = Base.genShapeF(xi);
                Jac = dS_dxi * obj.q0_elm;
                Jac_det = norm(Jac);

                e  = N * obj.q_elm;
                G      = Base.genG(e);    

                Ct_int = Base.genElmCt(obj.dq_elm, N, obj.Jrho, G, obj.numDOF);

                obj.Ct = obj.Ct + Ct_int * weights(g) * Jac_det;
            end

            if obj.isRotationFrameEnabled()
                obj.get_consistentRotFrameTangents();
                obj.Ct = obj.Ct + obj.Crot;
            end
        end

        function get_Kt(obj)
            nDim   = obj.numDOF;
            obj.Kt = zeros(nDim,nDim);

            obj.q_elm  = [];
            obj.dq_elm = [];

            for n = 1:obj.numNodes
                obj.q_elm  = [obj.q_elm;  obj.nodes{n}.qti ];
                obj.dq_elm = [obj.dq_elm; obj.nodes{n}.dqti];     
            end

            [weights, points] = Base.gaussQuadrature(obj.numGauss);

            for g = 1 : size(points,1)
                xi = points(g);
                obj.updateSectionProperties(xi);
                [~, N, dS_dxi, dN_dxi] = Base.genShapeF(xi);
                Jac = dS_dxi * obj.q0_elm;
                Jac_det = norm(Jac);

                S_prim = dS_dxi / Jac_det;
                N_prim = dN_dxi / Jac_det;
                e  = N * obj.q_elm;
                e0 = N * obj.q0_elm;
                r_prim  = S_prim * obj.q_elm;
                r0_prim = S_prim * obj.q0_elm;
                e_prim  = N_prim * obj.q_elm;
                e0_prim = N_prim * obj.q0_elm;
                G      = Base.genG(e);    
                lambd  = Base.genA(e);
                G0     = Base.genG(e0);
                lambd0 = Base.genA(e0);

                c = Base.genc(lambd0, r0_prim);

                gam  = lambd.'  * r_prim;
                gam0 = lambd0.' * r0_prim;
                kap  = 2 * G  * e_prim;
                kap0 = 2 * G0 * e0_prim; 

                Kt_int = Base.genElmExactKt(c, e, G, r_prim, N, e_prim, ...
                            lambd, S_prim, N_prim, obj.CA_mB, obj.CI_mB, gam, gam0,...
                            kap, kap0, obj.numDOF);

                obj.Kt = obj.Kt + Kt_int * weights(g) * Jac_det;
            end

            if obj.isRotationFrameEnabled()
                obj.get_consistentRotFrameTangents();
                obj.Kt = obj.Kt + obj.Krot;
            end
        end

        function get_Fi(obj)
            nDim   = obj.numDOF;
            obj.Fi = zeros(nDim,1);

            obj.q_elm  = [];
            obj.dq_elm = [];

            for n = 1:obj.numNodes
                obj.q_elm  = [obj.q_elm;  obj.nodes{n}.qti ];
                obj.dq_elm = [obj.dq_elm; obj.nodes{n}.dqti];     
            end

            [weights, points] = Base.gaussQuadrature(obj.numGauss);

            for g = 1 : size(points,1)
                xi = points(g);
                obj.updateSectionProperties(xi);
                [~, N, dS_dxi, ~] = Base.genShapeF(xi);
                Jac = dS_dxi * obj.q0_elm;
                Jac_det = norm(Jac);

                e  = N * obj.q_elm;
                G      = Base.genG(e);    

                Fi_int = Base.genElmFi(obj.dq_elm, N, obj.Jrho, G, obj.numDOF);

                obj.Fi = obj.Fi + Fi_int * weights(g) * Jac_det;
            end
        end

        function FrotPdimPart = get_FrotPdimPart(obj, axisR, mu0, ...
                stateOrder, parameterMultiIndex)
            %GET_FROTPDIMPART Rotating-frame force coefficient for PDIM.
            %
            % This is a dedicated semi-intrusive path. The complete Frot
            % used by time integration remains in genRotFrameTerms.
            nDim = obj.numDOF;
            FrotPdimPart = zeros(nDim,1);

            obj.q_elm  = [];
            obj.dq_elm = [];
            for n = 1:obj.numNodes
                obj.q_elm  = [obj.q_elm;  obj.nodes{n}.qti ];
                obj.dq_elm = [obj.dq_elm; obj.nodes{n}.dqti];
            end

            [weights, points] = Base.gaussQuadrature(obj.numGauss);
            for g = 1:size(points,1)
                xi = points(g);
                obj.updateSectionProperties(xi);
                [S, N, dS_dxi, ~] = Base.genShapeF(xi);
                Jac = dS_dxi * obj.q0_elm;
                Jac_det = norm(Jac);

                FrotPdimPart_int = Base.genElmFrotPdimPart( ...
                    obj.q_elm, obj.dq_elm, S, N, obj.Jrho, obj.rhoA, ...
                    axisR, mu0, stateOrder, parameterMultiIndex, nDim);
                FrotPdimPart = FrotPdimPart ...
                    + FrotPdimPart_int * weights(g) * Jac_det;
            end
        end

        function get_Fe(obj)
            nDim   = obj.numDOF;
            obj.Fe = zeros(nDim,1);

            obj.q_elm  = [];
            obj.dq_elm = [];

            for n = 1:obj.numNodes
                obj.q_elm  = [obj.q_elm;  obj.nodes{n}.qti ];
                obj.dq_elm = [obj.dq_elm; obj.nodes{n}.dqti];     
            end

            [weights, points] = Base.gaussQuadrature(obj.numGauss);

            for g = 1 : size(points,1)
                xi = points(g);
                obj.updateSectionProperties(xi);
                [~, N, dS_dxi, dN_dxi] = Base.genShapeF(xi);
                Jac = dS_dxi * obj.q0_elm;
                Jac_det = norm(Jac);

                S_prim = dS_dxi / Jac_det;
                N_prim = dN_dxi / Jac_det;
                e  = N * obj.q_elm;
                e0 = N * obj.q0_elm;
                r_prim  = S_prim * obj.q_elm;
                r0_prim = S_prim * obj.q0_elm;
                e_prim  = N_prim * obj.q_elm;
                e0_prim = N_prim * obj.q0_elm;
                G      = Base.genG(e);    
                lambd  = Base.genA(e);
                G0     = Base.genG(e0);
                lambd0 = Base.genA(e0);

                c = Base.genc(lambd0, r0_prim);

                gam  = lambd.'  * r_prim;
                gam0 = lambd0.' * r0_prim;
                kap  = 2 * G  * e_prim;
                kap0 = 2 * G0 * e0_prim; 

                Fe_int = Base.genElmFe(c, e, lambd, G, S_prim, r_prim, N, e_prim, ...
                            N_prim, obj.CA_mB, gam, gam0, obj.CI_mB, kap, kap0);

                obj.Fe = obj.Fe + Fe_int * weights(g) * Jac_det;
            end
        end

        function get_Fe1(obj)
            nDim   = obj.numDOF;
            obj.Fe1= zeros(nDim,1);

            obj.q_elm  = [];
            obj.dq_elm = [];

            for n = 1:obj.numNodes
                obj.q_elm  = [obj.q_elm;  obj.nodes{n}.qti ];
                obj.dq_elm = [obj.dq_elm; obj.nodes{n}.dqti];     
            end

            [weights, points] = Base.gaussQuadrature(obj.numGauss);

            for g = 1 : size(points,1)
                xi = points(g);
                obj.updateSectionProperties(xi);
                [~, N, dS_dxi, dN_dxi] = Base.genShapeF(xi);
                Jac = dS_dxi * obj.q0_elm;
                Jac_det = norm(Jac);

                S_prim = dS_dxi / Jac_det;
                N_prim = dN_dxi / Jac_det;
                e  = N * obj.q_elm;
                e0 = N * obj.q0_elm;
                r_prim  = S_prim * obj.q_elm;
                r0_prim = S_prim * obj.q0_elm;
                e_prim  = N_prim * obj.q_elm;
                G      = Base.genG(e);    
                lambd  = Base.genA(e);
                lambd0 = Base.genA(e0);

                c = Base.genc(lambd0, r0_prim);
                e0_prim = N_prim * obj.q0_elm;
                G0 = Base.genG(e0);
                kap0 = 2 * G0 * e0_prim;

                Fe1_int= Base.genElmFe1(c, e, lambd, G, S_prim, r_prim, N, e_prim, ...
                            N_prim, obj.CI_mB, kap0);

                obj.Fe1= obj.Fe1+ Fe1_int* weights(g) * Jac_det;
            end
        end

        function get_Fe2(obj)
            nDim   = obj.numDOF;
            obj.Fe2= zeros(nDim,1);

            obj.q_elm  = [];
            obj.dq_elm = [];

            for n = 1:obj.numNodes
                obj.q_elm  = [obj.q_elm;  obj.nodes{n}.qti ];
                obj.dq_elm = [obj.dq_elm; obj.nodes{n}.dqti];     
            end

            [weights, points] = Base.gaussQuadrature(obj.numGauss);

            for g = 1 : size(points,1)
                xi = points(g);
                obj.updateSectionProperties(xi);
                [~, N, dS_dxi, dN_dxi] = Base.genShapeF(xi);
                Jac = dS_dxi * obj.q0_elm;
                Jac_det = norm(Jac);

                S_prim = dS_dxi / Jac_det;
                N_prim = dN_dxi / Jac_det;
                e  = N * obj.q_elm;
                e0 = N * obj.q0_elm;
                r_prim  = S_prim * obj.q_elm;
                r0_prim = S_prim * obj.q0_elm;
                e_prim  = N_prim * obj.q_elm;
                G      = Base.genG(e);    
                lambd  = Base.genA(e);
                lambd0 = Base.genA(e0);

                c = Base.genc(lambd0, r0_prim);
                gam0 = lambd0.' * r0_prim;

                Fe2_int= Base.genElmFe2(c, e, lambd, G, S_prim, r_prim, N, e_prim, ...
                            N_prim, obj.CA_mB, gam0);

                obj.Fe2= obj.Fe2+ Fe2_int* weights(g) * Jac_det;
            end
        end

        function get_Fe3(obj)
            nDim   = obj.numDOF;
            obj.M0 = zeros(nDim,nDim);
            obj.Mnl= zeros(nDim,nDim);
            obj.Ct = zeros(nDim,nDim);
            obj.Kt = zeros(nDim,nDim);
            obj.Fi = zeros(nDim,1);
            obj.Fe = zeros(nDim,1);
            obj.Fe1= zeros(nDim,1);
            obj.Fe2= zeros(nDim,1);
            obj.Fe3= zeros(nDim,1);
            obj.Fe5= zeros(nDim,1);

            obj.q_elm  = [];
            obj.dq_elm = [];

            for n = 1:obj.numNodes
                obj.q_elm  = [obj.q_elm;  obj.nodes{n}.qti ];
                obj.dq_elm = [obj.dq_elm; obj.nodes{n}.dqti];     
            end

            [weights, points] = Base.gaussQuadrature(obj.numGauss);

            for g = 1 : size(points,1)
                xi = points(g);
                obj.updateSectionProperties(xi);
                [~, N, dS_dxi, dN_dxi] = Base.genShapeF(xi);
                Jac = dS_dxi * obj.q0_elm;
                Jac_det = norm(Jac);

                S_prim = dS_dxi / Jac_det;
                N_prim = dN_dxi / Jac_det;
                e  = N * obj.q_elm;
                e0 = N * obj.q0_elm;
                r_prim  = S_prim * obj.q_elm;
                r0_prim = S_prim * obj.q0_elm;
                e_prim  = N_prim * obj.q_elm;
                e0_prim = N_prim * obj.q0_elm;
                G      = Base.genG(e);    
                lambd  = Base.genA(e);
                G0     = Base.genG(e0);
                lambd0 = Base.genA(e0);

                c = Base.genc(lambd0, r0_prim);

                kap = 2 * G * e_prim;
                Fe3_int= Base.genElmFe3(c, e, lambd, G, S_prim, r_prim, N, e_prim, ...
                            N_prim, obj.CI_mB, kap);

                obj.Fe3= obj.Fe3+ Fe3_int* weights(g) * Jac_det;
            end
        end

        function get_Fe5(obj)
            nDim   = obj.numDOF;
            obj.Fe5= zeros(nDim,1);

            obj.q_elm  = [];
            obj.dq_elm = [];

            for n = 1:obj.numNodes
                obj.q_elm  = [obj.q_elm;  obj.nodes{n}.qti ];
                obj.dq_elm = [obj.dq_elm; obj.nodes{n}.dqti];     
            end

            [weights, points] = Base.gaussQuadrature(obj.numGauss);

            for g = 1 : size(points,1)
                xi = points(g);
                obj.updateSectionProperties(xi);
                [~, N, dS_dxi, dN_dxi] = Base.genShapeF(xi);
                Jac = dS_dxi * obj.q0_elm;
                Jac_det = norm(Jac);

                S_prim = dS_dxi / Jac_det;
                N_prim = dN_dxi / Jac_det;
                e  = N * obj.q_elm;
                e0 = N * obj.q0_elm;
                r_prim  = S_prim * obj.q_elm;
                r0_prim = S_prim * obj.q0_elm;
                e_prim  = N_prim * obj.q_elm;
                G      = Base.genG(e);    
                lambd  = Base.genA(e);
                lambd0 = Base.genA(e0);
                c = Base.genc(lambd0, r0_prim);

                gam  = lambd.'  * r_prim;

                Fe5_int= Base.genElmFe5(c, e, lambd, G, S_prim, r_prim, N, e_prim, ...
                            N_prim, obj.CA_mB, gam);

                obj.Fe5= obj.Fe5+ Fe5_int* weights(g) * Jac_det;
            end
        end
    end
end
