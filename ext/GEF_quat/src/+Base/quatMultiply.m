function e_result = quatMultiply(e1, e2)
% quaternion multiplication
% inut  ：q = [w, x, y, z] or  [w; x; y; z]
% output：q_result = q1 * q2 = [w, x, y, z]
% q1 * q2 = [w1*w2 - dot(v1,v2), w1*v2 + w2*v1 + cross(v1,v2)]

e1 = e1(:);
e2 = e2(:);
    
if length(e1) ~= 4 || length(e2) ~= 4
    error('input must be 4-D vector [w, x, y, z]');
end
    
w1 = e1(1); v1 = e1(2:4);
w2 = e2(1); v2 = e2(2:4);
    
w = w1 * w2 - dot(v1, v2);
v = w1 * v2 + w2 * v1 + cross(v1, v2);

e_result = [w; v];
end