...
    % System parameters;
    Qm = A*DeltaV*cap - sum( Cin ) * Na * e * V;
    disp("Initial charge [pC]: ");
    disp(A*DeltaV*cap*1e12);
...
    
...
    % Factors to reduce computation complexity
    potfact =  e * Na * V / ( A * cap );
    potfactadd = Qm  / ( A * cap );
    
    % Simulation
    for i = 2:numstep
        fluxes = fluxfact .* pot(i-1) .* ( C(i-1,:).*exp( expfact * pot(i-1) ) - Cout ) ./ ( 1 - exp( expfact * pot(i-1) ) );
        C(i,:) = C(i-1,:) + fluxes * cfact ;
        I(i,:) = ifact .* fluxes ;
        pot(i) = sum( C(i,:) ) * potfact + potfactadd ;
    end
    
...
    
...
    % plots
    figure();
    plot(timesteps, (sum( C,2 ) * e * Na * V + Qm ) * 1e12 );
    title("Charge inside the cell");
    xlabel("Time [s]");
    ylabel("Charge [pC]");
    
    disp("Equilibrium values:");
    disp("Charge [pC]:");
    disp((sum( C(numstep,:) ) * e * Na * V + Qm ) * 1e12 )
...