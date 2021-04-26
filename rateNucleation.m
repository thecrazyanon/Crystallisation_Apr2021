function B0 = rateNucleation(N1, N2, N3, N4, N5)
persistent N1_BETA
persistent N2_BETA
persistent N3_BETA
persistent N4_BETA
persistent N5_BETA 

switch nargin
    case 5
        N1_BETA = N1;
        N2_BETA = N2;
        N3_BETA = N3;
        N4_BETA = N4;
        N5_BETA = N5;
    case 3
        % Calculate primary and secondary nucleation rates given mass
        % fraction (N1), Temperature (N2) and second moment (N3)
        delConc = calcConc(N1) - calcSol(N2);
        if ((delConc > 0)&&(N5_BETA > N2))
            B0pri = N1_BETA*((delConc)^N2_BETA);
            B0sec = N3_BETA*N3*((delConc)^N4_BETA);
            B0 = B0pri + B0sec;
        else
            B0 = 0;
        end
end

end