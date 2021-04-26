function dY = analyticalTanks(t, Y, nTanks)
dY = zeros(nTanks,1);
Cin = 0.07;
RT = 900;
RT_tanks = RT/nTanks; 

for i = 1:nTanks
    if (i == 1)
        dY(i) = (Cin - Y(1))/RT_Tanks;
    else
        dY(i) = (Y(i-1) - Y(i))/RT_Tanks;
    end
    
end

end