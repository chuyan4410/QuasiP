function [weights,points] = gaussQuadrature(num_points)
switch num_points
    case 1
        points  = 0;
        weights = 2;
    case 2
        points  = [-sqrt(1/3); 
                   sqrt(1/3)];
        weights = [1; 1];
    case 3
        points  = [-sqrt(3/5); 
                    0; 
                    sqrt(3/5)];
        weights = [5/9; 
                   8/9; 
                   5/9];
    case 4
        points  = [-sqrt(3/7 + sqrt(120)/35); 
                   -sqrt(3/7 - sqrt(120)/35);
                    sqrt(3/7 - sqrt(120)/35); 
                    sqrt(3/7 + sqrt(120)/35)];
        weights = [1/2 - 5/(3*sqrt(120)); 
                   1/2 + 5/(3*sqrt(120));
                   1/2 + 5/(3*sqrt(120)); 
                   1/2 - 5/(3*sqrt(120))];
    case 5
        points  = [-(sqrt(5 + 2*sqrt(10/7)))/3; 
                   -(sqrt(5 - 2*sqrt(10/7)))/3; 
                    0;
                    (sqrt(5 - 2*sqrt(10/7)))/3; 
                    (sqrt(5+2*sqrt(10/7)))/3];
        weights = [ (322 - 13*sqrt(70))/900; 
                    (322 + 13*sqrt(70))/900; 
                    128/225;
                    (322 + 13*sqrt(70))/900; 
                    (322 - 13*sqrt(70))/900];
end