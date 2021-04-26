function CPmix = calcSpcHeat(CP1, CP2)
persistent CP1_BETA
persistent CP2_BETA

switch nargin
    case 2
        CP1_BETA = CP1;
        CP2_BETA = CP2;
    case 1
        CPmix = CP1_BETA + CP2_BETA*CP1;
end

end