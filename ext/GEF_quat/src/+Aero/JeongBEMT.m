classdef JeongBEMT < matlab.mixin.SetGetExactNames
    %JEONGBEMT BEM loads following the Jeong/Li style aeroelastic model.
    % The model returns generalized loads in the rotating beam frame. The
    % BEM normal/tangential load signs are mapped through forceSign.

    properties
        stationRadius
        chord
        twist
        airfoil
        airfoilData

        rho = 1.225
        bladeCount = 3
        rotorRadius = 63
        hubRadius = 1.5
        hubWind = 11.4
        omega = 12.1 * 2 * pi / 60
        pitch = 0

        numGauss = 3
        inductionMode = 'quasiSteadyIterative'
        inductionTolerance = 1e-4
        maxInductionIterations = 80
        inductionRelaxation = 1.0
        clampInduction = false
        axialInductionBounds = [-0.2, 0.95]
        angularInductionBounds = [-0.95, 0.95]

        includeTipLoss = true
        includeHubLoss = true
        includeGlauert = true
        glauertCriticalInduction = 0.2
        minimumLossFactor = 1e-3
        minimumVelocity = 1e-8

        loadScale = 1
        time = []
        primingTime = 0

        spanAxis = 1
        tangentialAxis = 2
        axialAxis = 3
        airfoilInterpolationMode = 'nearest'
        forceSign = struct('axial', -1, 'tangential', 1)
        pitchMomentSign = 1
        thetaDeformationSign = 1
        tangentialVelocitySign = 1
        momentDirectionMode = 'referenceSpan'
    end

    methods
        function obj = JeongBEMT(data, options)
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

                    section = obj.sectionState(elm, qElm, dqElm, S, N, ...
                        obj.isPrimingStep(stepID));
                    aero = obj.sectionAero(section);
                    [forceDensity, momentDensity] = obj.loadDensity(section, aero);

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
                stationSection = obj.sectionState(elm, qElm, dqElm, ...
                    Sstation, Nstation, obj.isPrimingStep(stepID));
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
                    'radius, chord, and airfoil must have the same length.');
            end

            validModes = {'quasiSteadyIterative'};
            if ~any(strcmpi(obj.inductionMode, validModes))
                error('Aero:InvalidInductionMode', ...
                    'The first JeongBEMT version supports only quasiSteadyIterative.');
            end

            validAirfoilModes = {'nearest', 'blend'};
            if ~any(strcmpi(obj.airfoilInterpolationMode, validAirfoilModes))
                error('Aero:InvalidAirfoilInterpolationMode', ...
                    'airfoilInterpolationMode must be nearest or blend.');
            end

            if ~isfield(obj.forceSign, 'axial') || ~isfield(obj.forceSign, 'tangential')
                error('Aero:InvalidForceSign', ...
                    'forceSign must define axial and tangential entries.');
            end

            if ~isscalar(obj.pitchMomentSign) || ~isnumeric(obj.pitchMomentSign)
                error('Aero:InvalidPitchMomentSign', ...
                    'pitchMomentSign must be a numeric scalar.');
            end

            if ~isscalar(obj.thetaDeformationSign) || ~isnumeric(obj.thetaDeformationSign)
                error('Aero:InvalidThetaDeformationSign', ...
                    'thetaDeformationSign must be a numeric scalar.');
            end

            if ~isscalar(obj.tangentialVelocitySign) || ~isnumeric(obj.tangentialVelocitySign)
                error('Aero:InvalidTangentialVelocitySign', ...
                    'tangentialVelocitySign must be a numeric scalar.');
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
                section.thetaDeformation = 0;
            else
                section.velocity = v;
                section.thetaDeformation = thetaRel(obj.spanAxis);
            end
            section.qCurrent = qCurrent;
            section.referenceFrame = Base.genA(qRef);
            section.spanDirection = section.referenceFrame(:, obj.spanAxis);
            section.chord = interp1(obj.stationRadius, obj.chord, section.radius, 'linear', 'extrap');
            section.twist = interp1(obj.stationRadius, obj.twist, section.radius, 'linear', 'extrap');
            [section.airfoilLeft, section.airfoilRight, section.airfoilWeight] = ...
                obj.airfoilBlend(section.radius);
            section.airfoil = (1 - section.airfoilWeight) * section.airfoilLeft ...
                + section.airfoilWeight * section.airfoilRight;
        end

        function aero = sectionAero(obj, section)
            vAxial = section.velocity(obj.axialAxis);
            vTan = section.velocity(obj.tangentialAxis);

            Wtan = 0;
            Waxial = 0;
            iterationCount = 0;

            for iter = 1:obj.maxInductionIterations
                iterationCount = iter;
                aeroTrial = obj.aeroFromInducedVelocity(section, vAxial, vTan, Wtan, Waxial);

                aForGlauert = obj.axialInductionFromVelocity(Waxial);
                fGlauert = obj.glauertCorrection(aForGlauert);
                denom = obj.inductionDenominator(section.radius, Waxial, fGlauert, aeroTrial.lossFactor);

                WtanNew = obj.bladeCount * aeroTrial.lift * sin(aeroTrial.phiAbs) / denom;
                WaxialNew = -obj.bladeCount * aeroTrial.lift * cos(aeroTrial.phiAbs) / denom;

                if obj.clampInduction
                    aNew = min(max(obj.axialInductionFromVelocity(WaxialNew), ...
                        obj.axialInductionBounds(1)), obj.axialInductionBounds(2));
                    aPrimeNew = min(max(obj.angularInductionFromVelocity(section.radius, WtanNew), ...
                        obj.angularInductionBounds(1)), obj.angularInductionBounds(2));
                    WaxialNew = -aNew * max(abs(obj.hubWind), obj.minimumVelocity);
                    WtanNew = aPrimeNew * max(abs(obj.omega * section.radius), obj.minimumVelocity);
                end

                WtanNext = (1 - obj.inductionRelaxation) * Wtan + obj.inductionRelaxation * WtanNew;
                WaxialNext = (1 - obj.inductionRelaxation) * Waxial + obj.inductionRelaxation * WaxialNew;

                if max(abs([WtanNext - Wtan, WaxialNext - Waxial])) < obj.inductionTolerance
                    Wtan = WtanNext;
                    Waxial = WaxialNext;
                    break
                end

                Wtan = WtanNext;
                Waxial = WaxialNext;
            end

            aero = obj.aeroFromInducedVelocity(section, vAxial, vTan, Wtan, Waxial);
            aero.inducedTangential = Wtan;
            aero.inducedAxial = Waxial;
            aero.a = obj.axialInductionFromVelocity(Waxial);
            aero.aPrime = obj.angularInductionFromVelocity(section.radius, Wtan);
            aero.glauertCorrection = obj.glauertCorrection(aero.a);
            aero.iterationCount = iterationCount;
        end

        function aero = aeroFromInducedVelocity(obj, section, vAxial, vTan, Wtan, Waxial)
            vRelAxial = obj.hubWind + vAxial + Waxial;
            vRelTan = obj.omega * section.radius ...
                + obj.tangentialVelocitySign * vTan + Wtan;
            vRelTanSafe = obj.withMinimumMagnitude(vRelTan);

            phi = atan2(vRelAxial, vRelTanSafe);
            phiAbs = max(abs(phi), 1e-6);
            thetaTotal = deg2rad(section.twist + obj.pitch) ...
                + obj.thetaDeformationSign * section.thetaDeformation;
            alpha = phi - thetaTotal;

            [cl, cd, cm, airfoilSec] = obj.airfoilCoefficientsAtRadius( ...
                section.radius, rad2deg(alpha));

            vRel = hypot(vRelAxial, vRelTan);
            lift = 0.5 * obj.rho * vRel^2 * section.chord * cl;
            drag = 0.5 * obj.rho * vRel^2 * section.chord * cd;

            aero.lift = lift;
            aero.drag = drag;
            aero.normalForce = lift * cos(phiAbs) + drag * sin(phiAbs);
            aero.tangentialForce = lift * sin(phiAbs) - drag * cos(phiAbs);
            aero.pitchMoment = 0.5 * obj.rho * vRel^2 * section.chord^2 * cm;
            aero.alpha = alpha;
            aero.phi = phi;
            aero.phiAbs = phiAbs;
            aero.thetaTotal = thetaTotal;
            aero.vRel = vRel;
            aero.vRelAxial = vRelAxial;
            aero.vRelTangential = vRelTan;
            aero.tipLoss = obj.tipLoss(section.radius, phiAbs);
            aero.hubLoss = obj.hubLoss(section.radius, phiAbs);
            aero.lossFactor = max(aero.tipLoss * aero.hubLoss, obj.minimumLossFactor);
            aero.cl = cl;
            aero.cd = cd;
            aero.cm = cm;
            aero.airfoil = airfoilSec.id;
            aero.airfoilLeft = airfoilSec.left;
            aero.airfoilRight = airfoilSec.right;
            aero.airfoilWeight = airfoilSec.weight;
        end

        function [forceDensity, momentDensity] = loadDensity(obj, section, aero)
            forceDensity = zeros(3, 1);
            forceDensity(obj.axialAxis) = obj.forceSign.axial * aero.normalForce;
            forceDensity(obj.tangentialAxis) = obj.forceSign.tangential * aero.tangentialForce;

            switch lower(obj.momentDirectionMode)
                case 'referencespan'
                    momentDirection = section.spanDirection;
                case 'axiscomponent'
                    momentDirection = zeros(3, 1);
                    momentDirection(obj.spanAxis) = 1;
                otherwise
                    error('Aero:InvalidMomentDirectionMode', ...
                        'momentDirectionMode must be referenceSpan or axisComponent.');
            end
            momentDensity = obj.pitchMomentSign * aero.pitchMoment * momentDirection;
        end

        function denom = inductionDenominator(obj, radius, Waxial, fGlauert, lossFactor)
            effectiveAxial = abs(obj.hubWind + fGlauert * Waxial);
            effectiveAxial = max(effectiveAxial, obj.minimumVelocity);
            denom = 4 * pi * obj.rho * radius * lossFactor * effectiveAxial;
            denom = max(denom, obj.minimumVelocity);
        end

        function a = axialInductionFromVelocity(obj, Waxial)
            a = -Waxial / max(abs(obj.hubWind), obj.minimumVelocity);
            a = real(a);
        end

        function aPrime = angularInductionFromVelocity(obj, radius, Wtan)
            aPrime = Wtan / max(abs(obj.omega * radius), obj.minimumVelocity);
            aPrime = real(aPrime);
        end

        function fg = glauertCorrection(obj, a)
            if ~obj.includeGlauert
                fg = 1;
                return
            end

            ac = obj.glauertCriticalInduction;
            if a <= ac || abs(a) < obj.minimumVelocity
                fg = 1;
            else
                ratio = ac / a;
                fg = ratio * (2 - ratio);
            end
        end

        function fTip = tipLoss(obj, radius, phiAbs)
            if ~obj.includeTipLoss
                fTip = 1;
                return
            end
            exponent = -obj.bladeCount / 2 * (obj.rotorRadius - radius) ...
                / max(radius * sin(phiAbs), obj.minimumVelocity);
            fTip = 2 / pi * acos(min(max(exp(exponent), 0), 1));
            fTip = max(fTip, obj.minimumLossFactor);
        end

        function fHub = hubLoss(obj, radius, phiAbs)
            if ~obj.includeHubLoss
                fHub = 1;
                return
            end
            if radius <= obj.hubRadius
                fHub = obj.minimumLossFactor;
                return
            end
            exponent = -obj.bladeCount / 2 * (radius - obj.hubRadius) ...
                / max(obj.hubRadius * sin(phiAbs), obj.minimumVelocity);
            fHub = 2 / pi * acos(min(max(exp(exponent), 0), 1));
            fHub = max(fHub, obj.minimumLossFactor);
        end

        function id = nearestAirfoil(obj, radius)
            [~, k] = min(abs(obj.stationRadius - radius));
            id = obj.airfoil(k);
        end

        function [idLeft, idRight, weightRight] = airfoilBlend(obj, radius)
            if radius <= obj.stationRadius(1)
                idLeft = obj.airfoil(1);
                idRight = idLeft;
                weightRight = 0;
                return
            end

            if radius >= obj.stationRadius(end)
                idLeft = obj.airfoil(end);
                idRight = idLeft;
                weightRight = 0;
                return
            end

            kLeft = find(obj.stationRadius <= radius, 1, 'last');
            kRight = kLeft + 1;
            span = obj.stationRadius(kRight) - obj.stationRadius(kLeft);
            weightRight = (radius - obj.stationRadius(kLeft)) / span;
            idLeft = obj.airfoil(kLeft);
            idRight = obj.airfoil(kRight);
        end

        function [cl, cd, cm, airfoil] = airfoilCoefficientsAtRadius(obj, radius, alphaDeg)
            switch lower(obj.airfoilInterpolationMode)
                case 'nearest'
                    id = obj.nearestAirfoil(radius);
                    [cl, cd, cm] = obj.airfoilCoefficients(id, alphaDeg);
                    airfoil.id = id;
                    airfoil.left = id;
                    airfoil.right = id;
                    airfoil.weight = 0;
                case 'blend'
                    [idLeft, idRight, weightRight] = obj.airfoilBlend(radius);
                    [clLeft, cdLeft, cmLeft] = obj.airfoilCoefficients(idLeft, alphaDeg);
                    [clRight, cdRight, cmRight] = obj.airfoilCoefficients(idRight, alphaDeg);
                    cl = (1 - weightRight) * clLeft + weightRight * clRight;
                    cd = (1 - weightRight) * cdLeft + weightRight * cdRight;
                    cm = (1 - weightRight) * cmLeft + weightRight * cmRight;
                    airfoil.id = (1 - weightRight) * idLeft + weightRight * idRight;
                    airfoil.left = idLeft;
                    airfoil.right = idRight;
                    airfoil.weight = weightRight;
                otherwise
                    error('Aero:InvalidAirfoilInterpolationMode', ...
                        'airfoilInterpolationMode must be nearest or blend.');
            end
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
            out.thetaDeformation = [];
            out.thetaTotal = [];
            out.normalForce = [];
            out.tangentialForce = [];
            out.pitchMoment = [];
            out.lift = [];
            out.drag = [];
            out.axialInduction = [];
            out.angularInduction = [];
            out.inducedAxial = [];
            out.inducedTangential = [];
            out.vRel = [];
            out.vRelAxial = [];
            out.vRelTangential = [];
            out.tipLoss = [];
            out.hubLoss = [];
            out.lossFactor = [];
            out.glauertCorrection = [];
            out.iterationCount = [];
            out.cl = [];
            out.cd = [];
            out.cm = [];
            out.airfoil = [];
            out.airfoilLeft = [];
            out.airfoilRight = [];
            out.airfoilWeight = [];
        end

        function out = storePointOutput(~, out, idx, section, aero)
            out.radius(idx, 1) = section.radius;
            out.alpha(idx, 1) = aero.alpha;
            out.phi(idx, 1) = aero.phi;
            out.thetaElastic(idx, 1) = section.thetaDeformation;
            out.thetaDeformation(idx, 1) = section.thetaDeformation;
            out.thetaTotal(idx, 1) = aero.thetaTotal;
            out.normalForce(idx, 1) = aero.normalForce;
            out.tangentialForce(idx, 1) = aero.tangentialForce;
            out.pitchMoment(idx, 1) = aero.pitchMoment;
            out.lift(idx, 1) = aero.lift;
            out.drag(idx, 1) = aero.drag;
            out.axialInduction(idx, 1) = aero.a;
            out.angularInduction(idx, 1) = aero.aPrime;
            out.inducedAxial(idx, 1) = aero.inducedAxial;
            out.inducedTangential(idx, 1) = aero.inducedTangential;
            out.vRel(idx, 1) = aero.vRel;
            out.vRelAxial(idx, 1) = aero.vRelAxial;
            out.vRelTangential(idx, 1) = aero.vRelTangential;
            out.tipLoss(idx, 1) = aero.tipLoss;
            out.hubLoss(idx, 1) = aero.hubLoss;
            out.lossFactor(idx, 1) = aero.lossFactor;
            out.glauertCorrection(idx, 1) = aero.glauertCorrection;
            out.iterationCount(idx, 1) = aero.iterationCount;
            out.cl(idx, 1) = aero.cl;
            out.cd(idx, 1) = aero.cd;
            out.cm(idx, 1) = aero.cm;
            out.airfoil(idx, 1) = aero.airfoil;
            out.airfoilLeft(idx, 1) = aero.airfoilLeft;
            out.airfoilRight(idx, 1) = aero.airfoilRight;
            out.airfoilWeight(idx, 1) = aero.airfoilWeight;
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

        function q = normalizeQuat(~, q)
            q = q(:);
            q = q / norm(q);
        end

        function value = withMinimumMagnitude(obj, value)
            if abs(value) < obj.minimumVelocity
                value = sign(value + eps) * obj.minimumVelocity;
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
