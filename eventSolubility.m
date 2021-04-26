function [position,isterminal,direction] = eventSolubility(t, y, inletState, modeContinuous, nTanks, ND_BETA)
under trial function not ready
delConc = calcConc(G1) - calcSol(G2);
crit_num_dens = 100;
position = y(1) - crit_num_dens; % The value that we want to be zero
isterminal = 1;  % Halt integration 
direction = 0;   % The zero can be approached from either direction
end

