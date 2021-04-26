function CONC = calcConc(C, dummy)
% Dummy variable added for ease of initialization
persistent CONC_BETA

switch nargin
    case 2
        CONC_BETA = C;
    case 1
        CONC = CONC_BETA*C*calcDensity(C);
end

end