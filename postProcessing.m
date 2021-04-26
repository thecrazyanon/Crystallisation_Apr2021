% Result variables
outConc = zeros(nReadings,nTanks);
outTemp = zeros(nReadings,nTanks);
outM0 = zeros(nReadings,nTanks);
outM1 = zeros(nReadings,nTanks);
outM2 = zeros(nReadings,nTanks);
outM3 = zeros(nReadings,nTanks);
outM4 = zeros(nReadings,nTanks);
var = zeros(nReadings,nTanks);
volTank = zeros(nReadings,nTanks);
holdUp = zeros(nReadings,nTanks);
d32 = zeros(nReadings,nTanks);
d10 = zeros(nReadings,nTanks);
SS = zeros(nReadings,nTanks);
SSR = zeros(nReadings,nTanks);
yield = zeros(nReadings,nTanks);
solConc = zeros(nReadings,nTanks);

jacketTemp = zeros(nReadings,1);

inletMassFrac = inletState(6);
inletConc = calcConc(inletMassFrac);

ssEPS = 1e-9;
tSS = 0;
oldConc = 0;
oldTime = 0;
T = (tRef^isDimensional)*T;
for i = 1:nReadings
    for j = 1:nTanks
        numTank = j;
        massFrac = calcMassFrac(Y(i,:),nVars,numTank);
        volTank(i,j) = Y(i,nVars*(j-1)+10);
        outConc(i,j) = (ConcRef^isDimensional)*calcConc(massFrac); 
        outTemp(i,j) = (TRef^isDimensional)*calcTemp(Y(i,:),nVars,numTank);
        SS(i,j) = (ConcRef^isDimensional)*(calcConc(massFrac) - calcSol(outTemp(i,j)/TRef));
        SSR(i,j) = (calcConc(massFrac) - calcSol(outTemp(i,j)/TRef))/calcSol(outTemp(i,j)/TRef);
        yield(i,j) = 100*((inletConc - calcConc(massFrac))/(inletConc - calcSol(outTemp(i,j)/TRef)));
        solConc(i,j) = (ConcRef^isDimensional)*calcSol(outTemp(i,j)/TRef);
        outM0(i,j) = ((M0Ref/VRef)^isDimensional)*Y(i,nVars*(j-1)+1)/volTank(i,j);
        outM1(i,j) = (((M0Ref*LRef)/VRef)^isDimensional)*Y(i,nVars*(j-1)+2)/volTank(i,j);        
        outM2(i,j) = (((M0Ref*(LRef^2))/VRef)^isDimensional)*Y(i,nVars*(j-1)+3)/volTank(i,j);
        outM3(i,j) = (((M0Ref*(LRef^3))/VRef)^isDimensional)*Y(i,nVars*(j-1)+4)/volTank(i,j);        
        outM4(i,j) = (((M0Ref*(LRef^4))/VRef)^isDimensional)*Y(i,nVars*(j-1)+5)/volTank(i,j);           
        holdUp(i,j) = (Y(i,nVars*(j-1)+7)/calcDensity(1))/((Y(i,nVars*(j-1)+7)/calcDensity(1))+(Y(i,nVars*(j-1)+6)/calcDensity(massFrac)));
        d32(i,j) = 1e6*outM3(i,j)/outM2(i,j);
        d10(i,j) = outM1(i,j)/outM0(i,j);
        var(i,j) = 1e6*sqrt((outM2(i,j)/outM0(i,j)) - (d10(i,j)^2));
        d10(i,j) = 1e6*d10(i,j);
    end
        if (oldConc ~= 0)
            dConcdt = abs((outConc(i,numTank)-oldConc)/(T(i)-oldTime));
            if (dConcdt < ssEPS)&&(tSS == 0)
                tSS = T(i);
            end
        end
        oldTime = T(i);
        oldConc = outConc(i,numTank);
    jacketTemp(i) = (TRef^isDimensional)*funcTempJacket(T(i)/(tRef^isDimensional));
end

% Normalize M0
dispM0 = outM0(:,tankDisplay);
% disp1M0 = outM0(:,tankDisplay-1);
for i = 1:nReadings
dispM0(i) = (dispM0(i)-dispM0(1))/(dispM0(nReadings)-dispM0(1));
% disp1M0(i) = (disp1M0(i)-disp1M0(1))/(disp1M0(nReadings)-disp1M0(1));
end

% yyaxis 'left';
% figure(1);
figure('units','normalized','outerposition',[0 0 1 1])
subplot(2,2,1);
% plot(T, outTemp(:,tankDisplay), 'r',T, outTemp(:,tankDisplay-1), 'g', T, jacketTemp, 'b');
% plot(T, var(:,tankDisplay), 'r',T, var(:,tankDisplay-1), 'g');
% plot(T, var(:,tankDisplay), 'r');
plot(t_exp, Tr_exp, 'x', T, outTemp(:,tankDisplay), 'r',T, outTemp(:,tankDisplay-1), 'g', T, jacketTemp, 'b');
% ylim([270 350]);
xlim([0 T(nReadings)]);
% ylim([0 1.5]) ;
xlabel('Time/s'); 
ylabel('Temperature (K)');
% 
% yyaxis 'right';
% plot(T, ss);

% xlim([0 T(nReadings)]);

hold off;
% figure(2);
subplot(2,2,2);
% plot(T, dispM0, 'r', T, disp1M0, 'g');
plot(T, dispM0, 'r');
xlim([0 T(nReadings)]);
xlabel('Time/s'); 
ylabel('Normalized 0th Moment');
% [AXA H1A H2A] = plotyy(t_exp,Tr_exp, t_exp,counts_exp);
% xlim(AXA(1),[0 T(nReadings)]);
% xlim(AXA(2),[0 T(nReadings)]);
% ylim(AXA(2),[0 1.5]);
% hold(AXA(1));
% hold(AXA(2));
% plot(AXA(1), T, outTemp(:,tankDisplay), 'r');
% plot(AXA(2), T, dispM0, 'r');

% hold off;
% figure(3);
% [AXB H1B H2B] = plotyy(T,d32, t_exp, counts_exp);
% xlim(AXB(1),[0 T(nReadings)]);
% xlim(AXB(2),[0 T(nReadings)]);
% hold(AXB(1));
% plot(AXB(1), t_exp, d_exp, 'r');

% figure(3);
% subplot(2,2,3);
% plot(T,d32,t_exp_d32, d32_PSD, 'r');
% xlim([0 T(nReadings)]);

subplot(2,2,3);
% plot(T,d10(:, tankDisplay), 'r', T,d10(:, tankDisplay-1), 'g');
plot(T,d10(:, tankDisplay), 'r');
xlim([0 T(nReadings)]);
xlabel('Time/s'); 
ylabel('Average diameter (micron)');

subplot(2,2,4);
% plot(T,outConc(:,tankDisplay-1),'g', T,outConc(:,tankDisplay),'r', t_conc_exp, conc_exp, '-o', T, solConc(:,tankDisplay), 'b');
% plot(T,outConc(:,tankDisplay-1),'g', T,outConc(:,tankDisplay),'r', T, solConc(:,tankDisplay), 'b');
plot(T,outConc(:,tankDisplay),'r', T, solConc(:,tankDisplay), 'b');
% plot(T,SSR(:,tankDisplay-1),'g', T,SSR(:,tankDisplay),'r');
plot(T,SSR(:,tankDisplay),'r');
xlim([0 T(nReadings)]);
xlabel('Time/s'); 
ylabel('SSR');

results = [d10(nReadings,tankDisplay) var(nReadings, tankDisplay) yield(nReadings, tankDisplay) tSS];
if (modeContinuous ~= 1)
    % Output for Batch
    V = VRef*Y(:,10);
    resultSummary = [T outM0 outM1 outM2 outM3 outM4 outConc outTemp V holdUp]; 
else
    % Output for continuous
    resultSummary = [T outM0 outM1 outM2 outM3 outM4 outConc outTemp holdUp]; 
end
% output = [T,jacketTemp,outTemp(:,tankDisplay),Tr_exp]; 
% output = [T Y];
% output_exp = [t_exp(1:25:2301) counts_exp(1:25:2301)];
% output_sim = [T dispM0];
% output_d32_exp = [t_exp_d32 d32_PSD]; 
% output_d32_sim = [T d32];
hold off;