function [S, N, dS_dxi, dN_dxi] = genShapeF(xi)
s1  = xi * (xi - 1) / 2;
s2  = 1 - xi^2;
s3  = xi * (xi + 1) / 2;
ds1 = xi - 0.5;
ds2 = -2 * xi;
ds3 = xi + 0.5;
S = [s1*eye(3) zeros(3,4) s2*eye(3) zeros(3,4) s3*eye(3) zeros(3,4)];
N = [zeros(4,3) s1*eye(4) zeros(4,3) s2*eye(4) zeros(4,3) s3*eye(4)];
dS_dxi = [ds1*eye(3) zeros(3,4) ds2*eye(3) zeros(3,4) ds3*eye(3) zeros(3,4)];
dN_dxi = [zeros(4,3) ds1*eye(4) zeros(4,3) ds2*eye(4) zeros(4,3) ds3*eye(4)];
end