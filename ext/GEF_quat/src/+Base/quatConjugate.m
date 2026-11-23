function eStar = quatConjugate(e)
eStar = zeros(4,1);
eStar(1)   =  e(1);
eStar(2:4) = -e(2:4);
end