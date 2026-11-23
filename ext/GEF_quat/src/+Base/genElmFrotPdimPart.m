function Fpart = genElmFrotPdimPart(q_elm, dq_elm, S, N, Jrho, ...
        rhoA, axisR, mu0, stateOrder, parameterMultiIndex, numDOF)
%GENELMFROTPDIMPART Homogeneous rotating-frame force coefficient for PDIM.
%
% The fixed-axis rotating-frame parameters are
%
%   mu = [omega; dotOmega] = mu0 + delta,
%
% and parameterMultiIndex = [ellOmega, ellDotOmega] selects the
% coefficient of deltaOmega^ellOmega * deltaDotOmega^ellDotOmega.
% Fpart is homogeneous of degree stateOrder in [q_elm; dq_elm].
%
% This routine is intentionally independent of the time-integration path.
% It does not assemble the complete Frot and does not use the prescribed
% rotating-frame quaternion.

axisR = axisR(:);
mu0 = mu0(:);
parameterMultiIndex = parameterMultiIndex(:).';

epsilon    = N * q_elm;
dotEpsilon = N * dq_elm;
G          = Base.genG(epsilon);
Lambda     = Base.genA(epsilon);

switch stateOrder
    case 1
        OmegaHatUnit = Base.genSkew(axisR);
        Gr1Unit = 2 * rhoA * S.' * OmegaHatUnit * S;
        Gr2Unit = rhoA * S.' * (OmegaHatUnit.' * OmegaHatUnit) * S;

        % Keep the two gyroscopic contributions in their original form
        % instead of relying on Gr1Unit being exactly skew-symmetric.
        Fomega = 0.5 * Gr1Unit * dq_elm ...
               - 0.5 * Gr1Unit.' * dq_elm;
        FdotOmega = 0.5 * Gr1Unit * q_elm;
        Fomega2 = -Gr2Unit * q_elm;

    case 3
        dotG = Base.genG(dotEpsilon);
        Fomega = 2 * N.' * (dotG.' * Jrho * Lambda.' * axisR);
        FdotOmega = 2 * N.' * (G.' * Jrho * Lambda.' * axisR);
        Fomega2 = zeros(numDOF,1);

        dG_de = Base.genElmdG_de();
        dA_de = Base.genElmdA_de(epsilon);
        dotA = sum(dA_de.*reshape(dotEpsilon,1,1,4),3);
        Fomega = Fomega + 2*N.'*G.'*Jrho*dotA.'*axisR;
        for j = 1:numDOF
            shapeColumn = N(:,j);
            G_j = zeros(3,4);
            A_j = zeros(3,3);
            for component = 1:4
                G_j = G_j + dG_de(:,:,component) * shapeColumn(component);
                A_j = A_j + dA_de(:,:,component) * shapeColumn(component);
            end

            Fomega(j) = Fomega(j) - 0.5 * ( ...
                4 * dotEpsilon.' * G_j.' * Jrho * Lambda.' * axisR + ...
                4 * dotEpsilon.' * G.'   * Jrho * A_j.'     * axisR);

            Fomega2(j) = -0.5 * ( ...
                axisR.' * A_j   * Jrho * Lambda.' * axisR + ...
                axisR.' * Lambda * Jrho * A_j.'   * axisR);
        end

    otherwise
        error('RotFramePDIM:UnsupportedStateOrder', ...
            'The rotating-frame PDIM force supports state orders 1 and 3 only.');
end

omega0 = mu0(1);
dotOmega0 = mu0(2);

if isequal(parameterMultiIndex, [0,0])
    Fpart = omega0 * Fomega + dotOmega0 * FdotOmega ...
        + omega0^2 * Fomega2;
elseif isequal(parameterMultiIndex, [1,0])
    Fpart = Fomega + 2 * omega0 * Fomega2;
elseif isequal(parameterMultiIndex, [0,1])
    Fpart = FdotOmega;
elseif isequal(parameterMultiIndex, [2,0])
    Fpart = Fomega2;
elseif sum(parameterMultiIndex) <= 2
    % The current structural model has no omega*dotOmega or dotOmega^2
    % rotating-frame force terms.
    Fpart = zeros(numDOF,1);
else
    error('RotFramePDIM:UnsupportedParameterOrder', ...
        'The rotating-frame force is at most quadratic in the parameters.');
end
end
