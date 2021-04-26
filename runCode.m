clear; 
clc;
close all;

% Input for non-dimensional groups
ND_BETA = xlsread('Input.xlsx','Non-Dimensional Groups','C2:C5');

% Input for non-dimensional groups for constitutive laws
ND_LAW_BETA = xlsread('Input.xlsx','Non-Dimensional Groups','C8:C15');

% Input for non-dimensional groups for auxilliary equations
ND_AUX_BETA = xlsread('Input.xlsx','Non-Dimensional Groups','C19:C25');

% Input for solver settings
solverSettings = xlsread('Input.xlsx','Solver Settings','C2:C10');

% Initialize reference variables
varRef = xlsread('Input.xlsx','Ref','C2:C9'); 
LRef = varRef(1);
VRef = varRef(2);
M0Ref = varRef(3);
tRef = varRef(4);
RhoRef = varRef(5);
CPRef = varRef(6);
TRef = varRef(7);
ConcRef = varRef(8); 

% Initialize Solver Setting Variables
nTanks = solverSettings(1);
tStart = solverSettings(2);
tEnd = solverSettings(3);
tReport = solverSettings(4);
errTolRel = solverSettings(5);
modeSwitchTol = solverSettings(6);
isDimensional = solverSettings(7);
constantDensity = solverSettings(8);
tankDisplay = solverSettings(9);

% Initialize Auxiliary Law functions
calcDensity(ND_AUX_BETA(2), ND_AUX_BETA(3), constantDensity);
calcSpcHeat(ND_AUX_BETA(4), ND_AUX_BETA(5));
calcConc(ND_AUX_BETA(1), 'init');
calcSol(ND_AUX_BETA(6), ND_AUX_BETA(7));
funcTempJacket(tRef, TRef);

trial = ConcRef*calcSol(1);

% Initialize Constitutive Law functions
rateGrowth(ND_LAW_BETA(1), ND_LAW_BETA(2), ND_LAW_BETA(3), ND_LAW_BETA(8));
rateNucleation(ND_LAW_BETA(4), ND_LAW_BETA(5), ND_LAW_BETA(6), ND_LAW_BETA(7), ND_LAW_BETA(8));

icRangeEnd = char(66 + nTanks);
icRange = strcat('C2:',icRangeEnd,'11');
initY = xlsread('Input.xlsx','Initial Conditions',icRange);

inletState = xlsread('Input.xlsx','System Properties','C10:C18'); 

% Minimum flow rate of 1ml/min for continuous mode of operation 
if(inletState(8) < modeSwitchTol)
    modeContinuous = 0;
else
    modeContinuous = 1;
end

tspan = tStart:tReport:tEnd;
% options = odeset('RelTol', errTolRel, 'Events', @eventSolubility);
options = odeset('RelTol', errTolRel);

[T, Y] = ode23s(@modelEquations, tspan, initY, options, inletState, modeContinuous, nTanks, ND_BETA);

[nReadings, nVars] = size(Y);
nVars = nVars/nTanks;

% Module to supply experimental data for comparison with model predictions
experimentalData;

% Module to process data for analysis
postProcessing;
