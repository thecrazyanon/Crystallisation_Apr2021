function G = rateGrowth(G1, G2, G3, G4)
persistent G1_BETA
persistent G2_BETA
persistent G3_BETA
persistent G4_BETA
switch nargin
    case 4
        G1_BETA = G1;
        G2_BETA = G2;
        G3_BETA = G3;
        G4_BETA = G4;
    case 2
        % Calculate growth rate give mass fraction (G1) and Temperature (G2)
        delConc = calcConc(G1) - calcSol(G2);
        if (delConc > 0)&&(G4_BETA > G2)
            G = G1_BETA*exp(-(G2_BETA/G2))*((delConc)^G3_BETA);
            % Growth on both sides
            G = 2*G;
        else
            G = 0;
        end
end
end