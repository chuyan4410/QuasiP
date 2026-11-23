function c = genc(A0, r0_prim)
e1 = A0.' * r0_prim;%eq.(21)
e  = [1;0;0];
c  = 1/(e1'*e);
end