function CONC_SOL = calcSol(S1, S2)
persistent S1_BETA
persistent S2_BETA

switch nargin
    case 2
        S1_BETA = S1;
        S2_BETA = S2;   
    case 1
        xD_sol = S1_BETA*exp(S2_BETA*S1);
        yD_sol = xD_sol/(1 + xD_sol);
        CONC_SOL = calcConc(yD_sol);
end

end