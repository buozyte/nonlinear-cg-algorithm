function [xfin, dfin] = stepsizeCG(maxsteps, d, x, rho, beta, sigma, delta1, delta2, f, g)
    l = 0;
    t = rho;    % beta^l = beta^0 = 1
    
    % Berechnung von x_(k+1) und d_(k+1) fuer l=0
    xnext = x + t * d;
    betaPR = (g(xnext)' * (g(xnext) - g(x))) / (norm(g(x))^2);
    dnext = -g(xnext) + betaPR * d;
    
    % Evaluierung der Bedingungen (a) und (b)
    a = (f(xnext) <= (f(x) - sigma * (t^2) * (norm(d)^2)));
    b1 = ((-delta2 * (norm(g(xnext))^2)) <= (g(xnext)' * dnext));
    b2 = ((g(xnext)' * dnext) <= (-delta1 * (norm(g(xnext))^2)));
    
    while (~a || ~b1 || ~b2)
        l = l + 1;
        
        if (l > maxsteps)
            error("A fitting stepsize couldn't be found in " + maxsteps + " steps.")
        end
        
        t = rho * (beta^l);
        
        % Berechnung von x_(k+1) und d_(k+1)
        xnext = x + t * d;
        betaPR = (g(xnext)' * (g(xnext) - g(x))) / (norm(g(x))^2);
        dnext = -g(xnext) + betaPR * d;
        
        
        % Falls man statt Polak-Ribiere Fletcher-Reeves oder
        % Hestenes-Stiefel verwenden moechte, waehlt man ein anderes beta
        
        %betaFR = (norm(g(xnext))^2) / (norm(g(x))^2);
        %betaHS = (g(xnext)' * (g(xnext)-g(x))) / ((g(xnext)-g(x))' * d);
        
        
        % Evaluierung der Bedingungen (a) und (b)
        a = (f(xnext) <= (f(x) - sigma * (t^2) * (norm(d)^2)));
        b1 = ((-delta2 * (norm(g(xnext))^2)) <= (g(xnext)' * dnext));
        b2 = ((g(xnext)' * dnext) <= (-delta1 * (norm(g(xnext))^2)));
    end
    
    xfin = xnext;
    dfin = dnext;
end