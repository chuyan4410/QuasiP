classdef HM3DCurvRodNode_quat < matlab.mixin.SetGetExactNames
    %HM3DCurvRodNode_quat : define the node class
    
    properties
        %------Constants---------------------------------------------------
        KI  = 1         % correction factor for mass and inertia moment
        numDOF = 7      % number of DOFs per node
        dofID           % DOF ID list
        Material        % Material object 
        crossSection    % CrossSection object 
        RM              % Residual magnetic flux density vector 
        Rot0            % Initial rotation matrix [3 x 3] 
        eInit           % Initial quaternion [4 x 1]  
        
        %------State quantities--------------------------------------------
        q               % generalized coordinates of nodes [7 x Steps]
        dq              % Linear velocity history [7 x Steps]
        ddq             % Linear acceleration history [7 x Steps]
        qti             % current generalized coordinates of nodes [7 x 1]
        dqti            % current linear velocity [7 x 1]
        ddqti           % current linear acceleration [7 x 1]
        
        %------Derived History---------------------------------------------
        axiVector       % Axial vector history (3 x Steps)
        lambd           % Effective rotation matrix list (Cell array of 3 x 3)
        H               % Quaternion to H matrix (Cell array)
        G               % Quaternion to G matrix (Cell array)
        W               % Material angular velocity (3 x Steps)
        dW              % Material angular acceleration (3 x Steps)
        w               % Spatial angular velocity (3 x Steps)
        dw              % Spatial angular acceleration (3 x Steps) 
        
        %------Loads list--------------------------------------------------
        AMList          % Applied magnetic flux density list (cell array 3 x 4)
        gList           % Body force list (cell array of 3 x 1)
        fList           % Distributed force list (cell array of 6 x 1)

        %------Current loads-----------------------------------------------
        AM              % Current applied magnetic flux density
        g               % Current body force
        f               % Current distributed force
    end
    
    methods
        function obj = HM3DCurvRodNode_quat(dofID, material, crossSection, ...
                RM, q0, dq0, ddq0, timenum, AMList, gList, fList, KI)
            % constant: initialization
            obj.KI           = KI;
            obj.dofID        = dofID;
            obj.Material     = material;
            obj.crossSection = crossSection;
            obj.RM           = RM;

            % rotation: cell: {e0, de0, dde0}
            obj.eInit  = q0(4:7);
            obj.Rot0   = Base.genA(q0(4:7));

            obj.AMList = AMList;
            obj.gList  = gList;
            obj.fList  = fList;  

            % state quantities: initialization
            obj.q        = zeros(obj.numDOF, timenum);
            obj.dq       = zeros(obj.numDOF, timenum);
            obj.ddq      = zeros(obj.numDOF, timenum);
            obj.q(:,1)   = q0;
            obj.dq(:,1)  = dq0;
            obj.ddq(:,1) = ddq0;
            
            % derived quantities: initialization
            obj.axiVector = cell(timenum,1);
            obj.lambd     = cell(timenum,1);
            obj.H         = cell(timenum,1);
            obj.G         = cell(timenum,1);
            obj.W         = cell(timenum,1);
            obj.dW        = cell(timenum,1);
            obj.w         = cell(timenum,1);
            obj.dw        = cell(timenum,1);

            obj.axiVector{1} = zeros(3,1);
            obj.lambd{1}     = eye(3);
            H0               = Base.genH(q0(4:7));
            G0               = Base.genG(q0(4:7));
            obj.H{1}         = H0;
            obj.G{1}         = G0;
            obj.W{1}         = 2 * G0 * dq0(4:7);
            obj.dW{1}        = 2 * G0 * ddq0(4:7);
            obj.w{1}         = 2 * H0 * dq0(4:7);
            obj.dw{1}        = 2 * H0 * ddq0(4:7);
        end
        
        function updateNode(obj, i, qti_node, dqti_node, ddqti_node)
            obj.AM = obj.AMList{i};
            obj.g  = obj.gList{i};
            obj.f  = obj.fList{i};

            % state quantities
            obj.q(:,i)   = qti_node;
            obj.dq(:,i)  = dqti_node;
            obj.ddq(:,i) = ddqti_node;
            obj.qti   = qti_node;
            obj.dqti  = dqti_node;
            obj.ddqti = ddqti_node;

            % derived quantities
            obj.axiVector{i} = Base.quat2theta(qti_node(4:7));
            obj.lambd{i}     = Base.genA(qti_node(4:7));
            obj.H{i}         = Base.genH(qti_node(4:7));
            obj.G{i}         = Base.genG(qti_node(4:7));
            obj.W{i}         = 2 * obj.G{i} * dqti_node(4:7);
            obj.dW{i}        = 2 * obj.G{i} * ddqti_node(4:7);
            obj.w{i}         = 2 * obj.H{i} * dqti_node(4:7);
            obj.dw{i}        = 2 * obj.H{i} * ddqti_node(4:7);
        end
    end
end

