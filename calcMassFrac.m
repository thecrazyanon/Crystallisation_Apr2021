function massFrac = calcMassFrac(y, nVars, numTank)
massFrac = y(nVars*(numTank-1)+8)/y(nVars*(numTank-1)+6);
end