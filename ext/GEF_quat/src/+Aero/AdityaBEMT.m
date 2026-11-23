classdef AdityaBEMT < matlab.mixin.SetGetExactNames
    %QUASISTEADYBEMT Quasi-steady BEMT loads for rotating-frame blade models.
    % The returned load is a generalized force vector in the same rotating
    % frame used by the quaternion beam state.

    properties
        stationRadius
        chord
        twist
        airfoil
        airfoilData
        aeroCenterOffset

        rho = 1.225
        bladeCount = 3
        rotorRadius = 63
        hubWind = 11.4
        omega = 12.1 * 2 * pi / 60
        pitch = 0
        numGauss = 3
        inductionTolerance = 1e-3
        maxInductionIterations = 50
        clampInduction = false
        axialInductionBounds = [-0.2, 0.95]
        angularInductionBounds = [-0.95, 0.95]
        includeTipLoss = true
        includeAeroCenterMoment = true
        loadScale = 1
        time = []
        primingTime = 0

        spanAxis = 1
        tangentialAxis = 2
        axialAxis = 3
    end

    methods
        function obj = AdityaBEMT(data, options)
            if nargin < 2
                options = struct();
            end
            obj.stationRadius = data.radius(:);
            obj.chord = data.chord(:);
            obj.twist = data.twist(:);
            obj.airfoil = data.airfoil(:);

            if isfield(data, 'airfoilData')
                obj.airfoilData = data.airfoilData;
            else
                obj.airfoilData = {};
            end

            if isfield(data, 'aeroCenterOffset')
                obj.aeroCenterOffset = data.aeroCenterOffset;
            else
                obj.aeroCenterOffset = zeros(numel(obj.stationRadius), 3);
            end

            optionNames = fieldnames(options);
            for k = 1:numel(optionNames)
                obj.(optionNames{k}) = options.(optionNames{k});
            end

            obj.validateInput();
        end

        function [forceSys, out] = evaluate(obj, solver, stepID)
            nDof = solver.dof_sys;
            forceSys = zeros(nDof, 1);

            gp = 0;
            out.gauss = obj.initPointOutput();
            out.station = obj.initPointOutput();

            [weights, points] = Base.gaussQuadrature(obj.numGauss);
            for e = 1:numel(solver.elements)
                elm = solver.elements{e};
                qElm = zeros(elm.numDOF, 1);
                dqElm = zeros(elm.numDOF, 1);
                for n = 1:elm.numNodes
                    idLocal = (n - 1) * elm.numDOFperN + (1:elm.numDOFperN);
                    qElm(idLocal) = elm.nodes{n}.qti;
                    dqElm(idLocal) = elm.nodes{n}.dqti;
                end

                forceElm = zeros(elm.numDOF, 1);
                for g = 1:numel(points)
                    xi = points(g);
                    [S, N, dS_dxi, ~] = Base.genShapeF(xi);
                    jac = norm(dS_dxi * qElm);
                    if jac <= eps
                        continue
                    end

                    section = obj.sectionState(elm, qElm, dqElm, S, N, obj.isPrimingStep(stepID));
                    aero = obj.sectionAero(section);

                    % shAlphaFrame = section.shadowFrame * obj.rotX((aero.alpha-section.thetaElastic));
                    shAlphaFrame = Base.genA(N * qElm) * obj.rotX((aero.alpha));
                    liftDensity = aero.lift * (-shAlphaFrame(:, obj.axialAxis));
                    dragDensity = aero.drag * (-shAlphaFrame(:, obj.tangentialAxis));
                    forceDensity = liftDensity + dragDensity;

                    momentDensity = -aero.pitchMoment * section.shadowFrame(:, obj.spanAxis);
                    if obj.includeAeroCenterMoment
                        rAC = section.shadowFrame * section.aeroCenterOffset;
                        momentDensity = momentDensity + cross(rAC, liftDensity);
                    end

                    H = Base.genH(section.qCurrent);
                    forceElm = forceElm + (S.' * forceDensity + 2 * N.' * H.' * momentDensity) ...
                        * weights(g) * jac;

                    gp = gp + 1;
                    out.gauss = obj.storePointOutput(out.gauss, gp, section, aero);
                    out.gauss.elementID(gp, 1) = e;
                    out.gauss.gaussID(gp, 1) = g;
                    out.gauss.xi(gp, 1) = xi;
                end

                forceSys(elm.dofID) = forceSys(elm.dofID) + forceElm;

                [Sstation, Nstation, ~, ~] = Base.genShapeF(0);
                stationSection = obj.sectionState(elm, qElm, dqElm, Sstation, Nstation, obj.isPrimingStep(stepID));
                stationAero = obj.sectionAero(stationSection);
                out.station = obj.storePointOutput(out.station, e, stationSection, stationAero);
                out.station.elementID(e, 1) = e;
                out.station.stationID(e, 1) = e;
                out.station.xi(e, 1) = 0;
            end

            scale = obj.currentLoadScale(stepID);
            forceSys = scale * forceSys;
            out.loadScale = scale;
            out.priming = obj.isPrimingStep(stepID);
            out.stepID = stepID;
        end
    end

    methods (Access = private)
        function validateInput(obj)
            n = numel(obj.stationRadius);
            if any([numel(obj.chord), numel(obj.twist), numel(obj.airfoil)] ~= n)
                error('Aero:InvalidInput', ...
                    'radius, chord, twist, and airfoil must have the same length.');
            end
            if size(obj.aeroCenterOffset, 1) ~= n || size(obj.aeroCenterOffset, 2) ~= 3
                error('Aero:InvalidInput', ...
                    'aeroCenterOffset must be an n-by-3 array.');
            end
        end

        function section = sectionState(obj, elm, qElm, dqElm, S, N, freezeAeroelasticFeedback)
            if nargin < 7
                freezeAeroelasticFeedback = false;
            end
            r0 = S * elm.q0_elm;
            r = S * qElm;
            v = S * dqElm;

            qRef = obj.normalizeQuat(N * elm.q0_elm);
            qCurrent = obj.normalizeQuat(N * qElm);
            qRel = Base.quatMultiply(Base.quatConjugate(qRef), qCurrent);
            qRel = obj.normalizeQuat(qRel);
            if qRel(1) < 0
                qRel = -qRel;
            end
            thetaRel = Base.quat2theta(qRel);

            section.radius = max(abs(r0(obj.spanAxis)), eps);
            section.position = r;
            if freezeAeroelasticFeedback
                section.velocity = zeros(3, 1);
            else
                section.velocity = v;
            end
            section.qCurrent = qCurrent;
            section.shadowFrame = Base.genA(qRef);
            if freezeAeroelasticFeedback
                section.thetaElastic = 0;
            else
                section.thetaElastic = thetaRel(obj.spanAxis);
            end
            section.chord = interp1(obj.stationRadius, obj.chord, section.radius, 'linear', 'extrap');
            section.twist = interp1(obj.stationRadius, obj.twist, section.radius, 'linear', 'extrap');
            section.airfoil = obj.nearestAirfoil(section.radius);
            section.aeroCenterOffset = obj.interpAeroCenter(section.radius);
        end

        function aero = sectionAero(obj, section)
            vAxial = section.velocity(obj.axialAxis);
            vTan = section.velocity(obj.tangentialAxis);

            a = 0;
            aPrime = 0;
            for iter = 1:obj.maxInductionIterations
                % Positive axial coordinate points toward the incoming wind.
                vRelAxial = obj.hubWind * (1 - a) + vAxial;
                vRelTan = obj.omega * section.radius * (1 + aPrime) + vTan;
                phi = atan2(vRelAxial, max(abs(vRelTan), eps));
                phiAbs = max(abs(phi), 1e-6);
                alpha = phi - deg2rad(section.twist + obj.pitch) - section.thetaElastic;

                [cl, cd, ~] = obj.airfoilCoefficients(section.airfoil, rad2deg(alpha));
                sigma = section.chord * obj.bladeCount / (2 * pi * section.radius);
                fTip = obj.tipLoss(section.radius, phiAbs);

                cn = cl * cos(phiAbs) + cd * sin(phiAbs);
                ct = cl * sin(phiAbs) - cd * cos(phiAbs);
                cn = obj.protectDenominator(cn);
                ct = obj.protectDenominator(ct);

                if a <= 0.3539
                    Hcorr = 1;
                else
                    Hcorr = 4 * a * (1 - a) / (0.6 + 0.61 * a + 0.79 * a^2);
                end

                aNew = 1 / (1 + 4 * fTip * sin(phiAbs)^2 / sigma / cn / Hcorr);
                aPrimeNew = 1 / (4 * fTip * sin(phiAbs) * cos(phiAbs) / sigma / ct - 1);
                aNew = real(aNew);
                aPrimeNew = real(aPrimeNew);
                if obj.clampInduction
                    aNew = min(max(aNew, obj.axialInductionBounds(1)), obj.axialInductionBounds(2));
                    aPrimeNew = min(max(aPrimeNew, obj.angularInductionBounds(1)), obj.angularInductionBounds(2));
                end

                if max(abs([aNew - a, aPrimeNew - aPrime])) < obj.inductionTolerance
                    a = aNew;
                    aPrime = aPrimeNew;
                    break
                end
                a = aNew;
                aPrime = aPrimeNew;
            end

            vRelAxial = obj.hubWind * (1 - a) + vAxial;
            vRelTan = obj.omega * section.radius * (1 + aPrime) + vTan;
            phi = atan2(vRelAxial, max(abs(vRelTan), eps));
            phiAbs = max(abs(phi), 1e-6);
            alpha = phi - deg2rad(section.twist + obj.pitch) - section.thetaElastic;
            [cl, cd, cm] = obj.airfoilCoefficients(section.airfoil, rad2deg(alpha));
            % fTip = obj.tipLoss(section.radius, phiAbs);

            vRel = hypot(vRelAxial, vRelTan);
            % lift = 0.5 * fTip * obj.rho * vRel^2 * section.chord * cl;
            % drag = 0.5 * fTip * obj.rho * vRel^2 * section.chord * cd;
            lift = 0.5 * obj.rho * vRel^2 * section.chord * cl;
            drag = 0.5 * obj.rho * vRel^2 * section.chord * cd;

            aero.lift = lift;
            aero.drag = drag;
            aero.normalForce = lift * cos(phiAbs) + drag * sin(phiAbs);
            aero.tangentialForce = lift * sin(phiAbs) - drag * cos(phiAbs);
            aero.pitchMoment = 0.5 * obj.rho * vRel^2 * section.chord^2 * cm;
            aero.alpha = alpha;
            aero.phi = phi;
            aero.a = a;
            aero.aPrime = aPrime;
            aero.cl = cl;
            aero.cd = cd;
            aero.cm = cm;
        end

        function fTip = tipLoss(obj, radius, phiAbs)
            if ~obj.includeTipLoss
                fTip = 1;
                return
            end
            exponent = -obj.bladeCount / 2 * (obj.rotorRadius - radius) / radius / sin(phiAbs);
            fTip = 2 / pi * acos(min(max(exp(exponent), 0), 1));
            fTip = max(fTip, 1e-3);
        end

        function id = nearestAirfoil(obj, radius)
            [~, k] = min(abs(obj.stationRadius - radius));
            id = obj.airfoil(k);
        end

        function offset = interpAeroCenter(obj, radius)
            offset = zeros(3, 1);
            for k = 1:3
                offset(k) = interp1(obj.stationRadius, obj.aeroCenterOffset(:,k), ...
                    radius, 'linear', 'extrap');
            end
        end

        function R = rotX(~, angle)
            c = cos(angle);
            s = sin(angle);
            R = [1, 0, 0;
                 0, c,-s;
                 0, s, c];
        end

        function [cl, cd, cm] = airfoilCoefficients(obj, airfoilID, alphaDeg)
            airfoilID = round(airfoilID);
            if airfoilID <= numel(obj.airfoilData) && ~isempty(obj.airfoilData{airfoilID})
                table = obj.airfoilData{airfoilID};
                cl = interp1(table(:,1), table(:,2), alphaDeg, 'linear', 'extrap');
                cd = interp1(table(:,1), table(:,3), alphaDeg, 'linear', 'extrap');
                cm = interp1(table(:,1), table(:,4), alphaDeg, 'linear', 'extrap');
                return
            end

            if airfoilID == 1
                cl = 0;
                cd = 0.5;
                cm = 0;
            elseif airfoilID == 2
                cl = 0;
                cd = 0.35;
                cm = 0;
            else
                error('Aero:MissingAirfoilData', ...
                    'Missing coefficient table for airfoil ID %d.', airfoilID);
            end
        end

        function out = initPointOutput(~)
            out.elementID = [];
            out.gaussID = [];
            out.stationID = [];
            out.xi = [];
            out.radius = [];
            out.alpha = [];
            out.phi = [];
            out.thetaElastic = [];
            out.normalForce = [];
            out.tangentialForce = [];
            out.pitchMoment = [];
            out.lift = [];
            out.drag = [];
            out.axialInduction = [];
            out.angularInduction = [];
            out.cl = [];
            out.cd = [];
            out.cm = [];
        end

        function out = storePointOutput(~, out, idx, section, aero)
            out.radius(idx, 1) = section.radius;
            out.alpha(idx, 1) = aero.alpha;
            out.phi(idx, 1) = aero.phi;
            out.thetaElastic(idx, 1) = section.thetaElastic;
            out.normalForce(idx, 1) = aero.normalForce;
            out.tangentialForce(idx, 1) = aero.tangentialForce;
            out.pitchMoment(idx, 1) = aero.pitchMoment;
            out.lift(idx, 1) = aero.lift;
            out.drag(idx, 1) = aero.drag;
            out.axialInduction(idx, 1) = aero.a;
            out.angularInduction(idx, 1) = aero.aPrime;
            out.cl(idx, 1) = aero.cl;
            out.cd(idx, 1) = aero.cd;
            out.cm(idx, 1) = aero.cm;
        end

        function scale = currentLoadScale(obj, stepID)
            if isa(obj.loadScale, 'function_handle')
                scale = obj.loadScale(stepID);
            elseif isscalar(obj.loadScale)
                scale = obj.loadScale;
            else
                k = min(stepID, numel(obj.loadScale));
                scale = obj.loadScale(k);
            end
        end

        function tf = isPrimingStep(obj, stepID)
            tf = false;
            if obj.primingTime <= 0
                return
            end
            if isempty(obj.time)
                scale = obj.currentLoadScale(stepID);
                tf = scale < 1;
                return
            end
            k = min(stepID, numel(obj.time));
            tf = obj.time(k) < obj.primingTime;
        end
    end

    methods (Access = private)
        function q = normalizeQuat(~, q)
            q = q(:);
            q = q / norm(q);
        end

        function val = protectDenominator(~, val)
            if abs(val) < 1e-8
                val = sign(val + eps) * 1e-8;
            end
        end
    end

    methods (Static)
        function data = load_NREL5MW_airfoilTables(folder)
            data = cell(8, 1);
            files = { ...
                '', ...
                '', ...
                'AfDu40_A17.txt', ...
                'AfDu35_A17.txt', ...
                'AfDu30_A17.txt', ...
                'AfDu25_A17.txt', ...
                'AfDu21_A17.txt', ...
                'AfNACA64_A17.txt'};
            for k = 1:numel(files)
                if isempty(files{k})
                    continue
                end
                fileName = fullfile(folder, files{k});
                data{k} = readmatrix(fileName, 'FileType', 'text');
            end
        end
    end
end
