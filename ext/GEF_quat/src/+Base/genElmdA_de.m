function dA_de = genElmdA_de(e)
e0 = e(1);
e1 = e(2);
e2 = e(3);
e3 = e(4);
dA_de = zeros(3,3,4);
dA_de(:,:,1) = 2 * [ e0 -e3  e2
                     e3  e0 -e1
                    -e2  e1  e0];
dA_de(:,:,2) = 2 * [ e1  e2  e3
                     e2 -e1 -e0
                     e3  e0 -e1];
dA_de(:,:,3) = 2 * [-e2  e1  e0
                     e1  e2  e3
                    -e0  e3 -e2];
dA_de(:,:,4) = 2 * [-e3 -e0  e1
                     e0 -e3  e2
                     e1  e2  e3];
end