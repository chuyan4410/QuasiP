function theta = quat2theta(q, eps)
%QUAT2THETA Quaternion to axis vector
%
%   theta = quat2theta(q)
%   theta = quat2theta(q, eps)
%
%   q   : quaternion [4x1], q = [q0; q1; q2; q3]
%   eps : tolerance (optional)
%
%   theta : rotation vector [3x1]

    if nargin < 2
        eps = 1e-12;
    end

    % Normalize quaternion
    q = q / norm(q);

    % Rotation angle
    a = 2 * acos(q(1));

    if a <= eps
        theta = [0; 0; 0];
    else
        theta = q(2:4) * a / sin(a/2);
    end
end
