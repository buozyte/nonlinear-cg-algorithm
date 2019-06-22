function [fin, x_steps, counter] = modPRalg(maxsteps, x0, eps, beta, sigma, delta1, delta2, f, g)

    % Test, ob alle Parameter im passanden Bereich sind
    if (eps<0 || beta<=0 || beta>=1 || sigma<=0 || sigma>=1 || delta1<=0 || delta1>=1 || delta2 <=1)
        error("one of your chosen parameters does not fit the requierements")
    end
    
    % Solange Abbruchbedingung nicht erfuellt ist, 
    % soll x_(k+1) berechnet werden
    x = x0;
    d = -g(x);
    counter = 0;
    
    x_steps = x';
    
    while (norm(g(x)) > eps)
        counter = counter + 1;
        
        if (counter > maxsteps)
            error("A minimum couldn't be found in " + maxsteps + " steps.")
        end
        
        rho = norm(g(x)' * d) / (norm(d)^2);
        
        % Berechnung der x_k+1 und d_k+1
        [x, d] = stepsizePR(maxsteps, d, x, rho, beta, sigma, delta1, delta2, f, g);
        
        x_steps = [x_steps; x'];
    end
    fin = x;
end