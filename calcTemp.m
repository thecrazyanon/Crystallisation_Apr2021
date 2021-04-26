function unitTemp = calcTemp(y,nVars,numTank)
massFrac = calcMassFrac(y, nVars, numTank);
unitTemp = (y(nVars*(numTank-1)+9)/(y(nVars*(numTank-1)+6)*calcSpcHeat(massFrac) + y(nVars*(numTank-1)+7)*calcSpcHeat(1)));
end