function [fin, x_steps, counter] = algFRHS(maxsteps, x0, eps, eta, f, g, alg)

    % Test, ob alle Parameter im passanden Bereich sind
    if (eps<0 || eta<=1/2 || eta>=1)
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
        
        % Berechnung der Schrittweite
        t = strongerPW(d, f, g, maxsteps, x);
        
        oldx = x;
        x = x + t*d;
        
        switch alg
            case "FR"
                beta = norm(g(x))^2 / norm(g(oldx))^2;
            case "HS"
                beta = ((g(x)-g(oldx))'*g(x))/((g(x)-g(oldx))' *d);
            case "PR"
                beta = (g(x)'*(g(x) - g(oldx))) / (norm(g(oldx))^2);
            case "Grad"
                beta = 0;
        end
        
        d = -g(x) + beta * d;
        
        x_steps = [x_steps; x'];
    end
    fin = x;
end

