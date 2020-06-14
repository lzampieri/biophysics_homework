function depolarization()
    
    % System parameters; [ K+ Na+ ]
    P = [ 10e-3 0.2e-3 ];
    cap = 1e-2;
    R = 10.6e-6;
    A = 1400e-12;
    V = 5000e-18;
    Cin = [ 155 12]; % mol / m^3
    Cout = [ 4 145]; % mol / m^3
    DeltaV = -90e-3;
    z = [ 1 1 ];
    e = 1.6e-19;
    kBT = 4.27e-21; % @36°C
    Na = 6.022e23;
    
    % Print initial fluxes
    disp("K+ Na+")
    disp("Initial fluxes:")
    disp( P.*z*e*DeltaV/kBT .* ( Cin.*exp( z*e*DeltaV/kBT ) - Cout ) ./ ( 1 - exp( z*e*DeltaV/kBT ) ) )
    
    % Simulation parameters
    timestep = 1e-9;
    tottime = 0.02;
    numstep = floor( tottime / timestep );
    C = zeros(numstep,2);
    C(1,:) = Cin;
    I = zeros(numstep,2);
    pot = zeros(numstep,1);
    pot(1) = DeltaV;
    fluxes = [ 0 0 ];
    
    % Factors to reduce computation complexity
    fluxfact = P.*z*e/kBT;
    expfact = z*e/kBT;
    cfact = A * timestep / V;
    ifact =  z*e*Na*A;
    potfact = e*Na / cap * timestep;
    
    % Simulation
    for i = 2:numstep
        fluxes = fluxfact .* pot(i-1) .* ( C(i-1,:).*exp( expfact * pot(i-1) ) - Cout ) ./ ( 1 - exp( expfact * pot(i-1) ) );
        C(i,:) = C(i-1,:) + fluxes * cfact ;
        I(i,:) = ifact.* fluxes ;
        pot(i) = pot(i-1) + potfact * sum(  z.* fluxes ) ;
    end
    
    % plots
    timesteps = 0:timestep:(numstep-1)*timestep;
    opengl software;
    figure();
    plot(timesteps,pot*1e3);
    hold on;
    title("Potential and concentrations inside the cell");
    xlabel("Time [s]");
    ylabel("Potential [mV]");
    yyaxis right;
    ylabel("Concentration [mM]");
    plot(timesteps,C);
    legend('\Delta V','K+','Na+');
    hold off;
    
    figure();
    plot(timesteps,[I, sum(I,2)] * 1e6);
    title("Current inside the cell");
    xlabel("Time [s]");
    ylabel("Current [\mu A]");
    legend('K+','Na+','Total');
    
    disp("Equilibrium values:");
    disp("K+ Na+");
    disp("Concentrations [mM]:");
    disp(C(numstep,:));
    disp("Potential [V]:");
    disp(pot(numstep));
    disp("Final fluxes:")
    disp( P.*z*e*pot(numstep)/kBT .* ( C(numstep,:).*exp( z*e*pot(numstep)/kBT ) - Cout ) ./ ( 1 - exp( z*e*pot(numstep)/kBT ) ) )

end