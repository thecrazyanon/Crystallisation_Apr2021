function RHOmix = calcDensity(D1, D2, D3)
persistent D1_BETA
persistent D2_BETA
persistent constantDensity

switch nargin
    case 3
        D1_BETA = D1;
        D2_BETA = D2;
        constantDensity = D3;
    case 1
        if (D1 == 1) % Density of Pure solid phase
            RHOmix = D1_BETA + D2_BETA; 
        else % Density of liquid phase
            if (constantDensity == 1)
                RHOmix = D1_BETA;
            else
                RHOmix = D1_BETA + D2_BETA*D1;
            end
        end            
end
end