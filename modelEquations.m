function dY = modelEquations(t, Y, inletState, modeContinuous, nTanks, ND_BETA)
[nVar, ~] = size(Y);
nVar = nVar/nTanks;
dY = zeros(nTanks*nVar,1);

BETA_1 = ND_BETA(1);
BETA_2 = ND_BETA(2);
BETA_3 = ND_BETA(3);
BETA_4 = ND_BETA(4);

massFracUnit = zeros(nTanks, 1);
tempUnit = zeros(nTanks, 1);
RHOmixUnit = zeros(nTanks,1);
CPmixUnit = zeros(nTanks,1);
outFlowUnit = zeros(nTanks,1);
holdUpUnit = zeros(nTanks,1);
rateMCUnit = zeros(nTanks,1);

inletM0 = inletState(1);
inletM1 = inletState(2);
inletM2 = inletState(3);
inletM3 = inletState(4);
inletM4 = inletState(5);
inletMassFrac = inletState(6);
inletTemp = inletState(7);
inletFlow = inletState(8);
tempEnv = inletState(9);

tempJacket = funcTempJacket(t);

inletRHOmix = calcDensity(inletMassFrac);
inletCPmix = calcSpcHeat(inletMassFrac);

% Solid phase properties
RHOsol = calcDensity(1);
CPsol = calcSpcHeat(1);

inletHoldUp = BETA_1*inletM3;

for i = 1:nTanks
    massLiqUnit = Y(6+(i-1)*nVar);
    dissMassLiqUnit = Y(8+(i-1)*nVar);
    massSolUnit = Y(7+(i-1)*nVar);
    heatUnit = Y(9+(i-1)*nVar);
    M2Unit = Y(3+(i-1)*nVar);
    M3Unit = Y(4+(i-1)*nVar);
    volUnit = Y(10+(i-1)*nVar);
    
    massFracUnit(i) = dissMassLiqUnit/massLiqUnit;
    RHOmixUnit(i) = calcDensity(massFracUnit(i));
    CPmixUnit(i) = calcSpcHeat(massFracUnit(i));
    holdUpUnit(i) = BETA_1*M3Unit;
    
    tempUnit(i) = heatUnit/(massLiqUnit*CPmixUnit(i) + massSolUnit*CPsol);
    rateMCUnit(i) = 3*BETA_1*rateGrowth(massFracUnit(i), tempUnit(i))*M2Unit*calcDensity(1)*volUnit;
end

dmassFracUnitdt = zeros(nTanks,1);
BETA_D2 = calcDensity(1) - calcDensity(0);
if (modeContinuous == 1)
    % Calculating outlet velocity for each tank
    for i = 1:nTanks
        massLiqUnit = Y(6+(i-1)*nVar);
        if (i == 1)
            dmassFracUnitdt(i) = (1/massLiqUnit)*((inletFlow*RHOmixUnit(i)*(1-inletHoldUp)*(inletMassFrac-massFracUnit(i))) ...
                                   -rateMCUnit(i)*(1-massFracUnit(i)));
            outFlowUnit(i) = (inletFlow*(inletRHOmix/RHOmixUnit(i))*(1-inletHoldUp)+inletFlow*inletHoldUp)...
                                   -(BETA_D2*massLiqUnit/(RHOmixUnit(i)^2))*dmassFracUnitdt(i) ...
                                   -(rateMCUnit(i)*((1/RHOmixUnit(i))-(1/RHOsol)));
        else
            dmassFracUnitdt(i) = (1/massLiqUnit)*((outFlowUnit(i-1)*RHOmixUnit(i)*(1-holdUpUnit(i-1))*(massFracUnit(i-1)-massFracUnit(i))) ...
                                   -rateMCUnit(i)*(1-massFracUnit(i)));
            outFlowUnit(i) = (outFlowUnit(i-1)*(RHOmixUnit(i-1)/RHOmixUnit(i))*(1-holdUpUnit(i-1))+outFlowUnit(i-1)*holdUpUnit(i-1))...
                                   -(BETA_D2*massLiqUnit/(RHOmixUnit(i)^2))*dmassFracUnitdt(i) ...
                                   -(rateMCUnit(i)*((1/RHOmixUnit(i))-(1/RHOsol)));                               
        end
    end
end

for i = 1:nTanks
    M0Unit = Y(1 + (i-1)*nVar)/Y(10 + (i-1)*nVar);
    M1Unit = Y(2 + (i-1)*nVar)/Y(10 + (i-1)*nVar);
    M2Unit = Y(3 + (i-1)*nVar)/Y(10 + (i-1)*nVar);
    M3Unit = Y(4 + (i-1)*nVar)/Y(10 + (i-1)*nVar);    
    M4Unit = Y(5 + (i-1)*nVar)/Y(10 + (i-1)*nVar); 
    MLiqUnit = Y(6 + (i-1)*nVar)/Y(10 + (i-1)*nVar); 
    rateNuclUnit = rateNucleation(massFracUnit(i), tempUnit(i), M2Unit);
    rateGrowthUnit = rateGrowth(massFracUnit(i), tempUnit(i));
    volUnit = Y(10 + (i-1)*nVar);
    if (i == 1)
        % M0
        dY(1 + (i-1)*nVar) = rateNuclUnit*volUnit + modeContinuous*(inletFlow*inletM0-outFlowUnit(i)*M0Unit);
        % M1
        dY(2 + (i-1)*nVar) = rateGrowthUnit*M0Unit*volUnit + modeContinuous*(inletFlow*inletM1-outFlowUnit(i)*M1Unit);
        % M2
        dY(3 + (i-1)*nVar) = 2*rateGrowthUnit*M1Unit*volUnit + modeContinuous*(inletFlow*inletM2-outFlowUnit(i)*M2Unit);
        % M3
        dY(4 + (i-1)*nVar) = 3*rateGrowthUnit*M2Unit*volUnit + modeContinuous*(inletFlow*inletM3-outFlowUnit(i)*M3Unit);
        % M4
        dY(5 + (i-1)*nVar) = 4*rateGrowthUnit*M3Unit*volUnit + modeContinuous*(inletFlow*inletM4-outFlowUnit(i)*M4Unit);
        % Mliq
        dY(6 + (i-1)*nVar) = modeContinuous*(inletFlow*inletRHOmix*(1-inletHoldUp) - outFlowUnit(i)*RHOmixUnit(i)*(1-holdUpUnit(i))) - rateMCUnit(i);
        % MSol
        dY(7 + (i-1)*nVar) = modeContinuous*(inletFlow*RHOsol*inletHoldUp - outFlowUnit(i)*RHOsol*holdUpUnit(i)) + rateMCUnit(i);
        % MLiqDiss
        dY(8 + (i-1)*nVar) = modeContinuous*(inletFlow*inletRHOmix*(1-inletHoldUp)*inletMassFrac - outFlowUnit(i)*RHOmixUnit(i)*(1-holdUpUnit(i))*massFracUnit(i)) ...
                         - rateMCUnit(i);
        % Heat Balance
        dY(9 + (i-1)*nVar) = modeContinuous*(inletFlow*inletRHOmix*(1-inletHoldUp)*inletCPmix*inletTemp - outFlowUnit(i)*RHOmixUnit(i)*(1-holdUpUnit(i))*CPmixUnit(i)*tempUnit(i)) ...
                        + modeContinuous*(inletFlow*RHOsol*inletHoldUp*CPsol*inletTemp - outFlowUnit(i)*RHOsol*holdUpUnit(i)*CPsol*tempUnit(i)) ...
                        + BETA_2*rateMCUnit(i) + BETA_3*(tempJacket-tempUnit(i)) - BETA_4*(tempUnit(i)-tempEnv);
        % Volume Balance
        if (modeContinuous == 1)
            dY(10 + (i-1)*nVar) = 0;
        else
            % Batch mode
            dmassFracdt = -(rateMCUnit(i)*(1-massFracUnit(i)))/MLiqUnit;
            dY(10 + (i-1)*nVar) = -(rateMCUnit(i)/RHOmixUnit(i))-((BETA_D2*MLiqUnit)/(RHOmixUnit(i)^2))*dmassFracdt + rateMCUnit(i)/RHOsol;
        end
    else
        M0PrevUnit = Y(1 + (i-2)*nVar)/Y(10 + (i-2)*nVar);
        M1PrevUnit = Y(2 + (i-2)*nVar)/Y(10 + (i-2)*nVar);
        M2PrevUnit = Y(3 + (i-2)*nVar)/Y(10 + (i-2)*nVar);
        M3PrevUnit = Y(4 + (i-2)*nVar)/Y(10 + (i-2)*nVar);    
        M4PrevUnit = Y(5 + (i-2)*nVar)/Y(10 + (i-2)*nVar); 
        % M0
        dY(1 + (i-1)*nVar) = rateNuclUnit*volUnit + modeContinuous*(outFlowUnit(i-1)*M0PrevUnit-outFlowUnit(i)*M0Unit);
        % M1
        dY(2 + (i-1)*nVar) = rateGrowthUnit*M0Unit*volUnit + modeContinuous*(outFlowUnit(i-1)*M1PrevUnit-outFlowUnit(i)*M1Unit);
        % M2
        dY(3 + (i-1)*nVar) = 2*rateGrowthUnit*M1Unit*volUnit + modeContinuous*(outFlowUnit(i-1)*M2PrevUnit-outFlowUnit(i)*M2Unit);
        % M3
        dY(4 + (i-1)*nVar) = 3*rateGrowthUnit*M2Unit*volUnit + modeContinuous*(outFlowUnit(i-1)*M3PrevUnit-outFlowUnit(i)*M3Unit);
        % M4
        dY(5 + (i-1)*nVar) = 4*rateGrowthUnit*M3Unit*volUnit + modeContinuous*(outFlowUnit(i-1)*M4PrevUnit-outFlowUnit(i)*M4Unit);
        % Mliq
        dY(6 + (i-1)*nVar) = modeContinuous*(outFlowUnit(i-1)*RHOmixUnit(i-1)*(1-holdUpUnit(i-1)) - outFlowUnit(i)*RHOmixUnit(i)*(1-holdUpUnit(i))) - rateMCUnit(i);
        % MSol
        dY(7 + (i-1)*nVar) = modeContinuous*(outFlowUnit(i-1)*RHOsol*holdUpUnit(i-1) - outFlowUnit(i)*RHOsol*holdUpUnit(i)) + rateMCUnit(i);
        % MLiqDiss
        dY(8 + (i-1)*nVar) = modeContinuous*(outFlowUnit(i-1)*RHOmixUnit(i-1)*(1-holdUpUnit(i-1))*massFracUnit(i-1) ...
                        - outFlowUnit(i)*RHOmixUnit(i)*(1-holdUpUnit(i))*massFracUnit(i)) ...
                        - rateMCUnit(i);
        % Heat Balance
        dY(9 + (i-1)*nVar) = modeContinuous*(outFlowUnit(i-1)*RHOmixUnit(i-1)*(1-holdUpUnit(i-1))*CPmixUnit(i-1)*tempUnit(i-1) ...
                        - outFlowUnit(i)*RHOmixUnit(i)*(1-holdUpUnit(i))*CPmixUnit(i)*tempUnit(i) ...
                        + outFlowUnit(i-1)*RHOsol*holdUpUnit(i-1)*CPsol*tempUnit(i-1) ...
                        - outFlowUnit(i)*RHOsol*holdUpUnit(i)*CPsol*tempUnit(i)) ...
                        + BETA_2*rateMCUnit(i) + BETA_3*(tempJacket-tempUnit(i)) - BETA_4*(tempUnit(i)-tempEnv);
        % Volume Balance
        if (modeContinuous == 1)
            dY(10 + (i-1)*nVar) = 0;
        else
            % Batch mode
            dmassFracdt = -(rateMCUnit(i)*(1-massFracUnit(i)))/MLiqUnit;
            dY(10 + (i-1)*nVar) = -(rateMCUnit(i)/RHOmixUnit(i))-((BETA_D2*MLiqUnit)/(RHOmixUnit(i)^2))*dmassFracdt + rateMCUnit(i)/RHOsol;
        end
    end
end


end