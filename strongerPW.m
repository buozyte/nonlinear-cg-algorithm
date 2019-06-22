function [t] = strongerPW(d, f, g, maxsteps, x)


    % phi(t) = f(x + td)
    % psi(t) = phi(t) - phi(0) - sigma*t*g(x)'d

    
    sigma = 10^(-4);            % 0 < sigma < 1/2
    rho = 0.01;                 % sigma < rho < 1
    
    
    % PHASE A
    t = 1;                      % t0 > 0
    gamma = 1.1;                % gamma > 1
    
    counter = 0;
    re1 = true;
    while re1
    	if (counter > maxsteps)
            error("A fitting stepsize couldn't be found in " + maxsteps + " steps in phase A.")
    	end
            
        phi = f(x + t*d);
        dphi = g(x + t*d)'*d;
        psi = phi - f(x) - sigma*t*(g(x)'*d);
        dpsi = dphi - sigma*(g(x)'*d);
            
        if psi >= 0
            a = 0;
            b = t;
            re1 = false;
        elseif psi < 0 && (norm(dpsi) <= (rho - sigma)*norm(g(x)'*d))
            return
        elseif psi < 0 && dpsi > 0
            a = t;
            b = 0;
            re1 = false;
        elseif psi < 0 && dpsi < 0
            t = gamma*t;
            counter = counter + 1;
        end
    end
    
    
    
    % PHASE B
    tau1 = 0.1; 
    tau2 = 0.1;        % 0 < tau1, tau2 <= 1/2
    
    counter = 0;
    while true
        if (counter > maxsteps)
        	error("A fitting stepsize couldn't be found in " + maxsteps + " steps in phase B.")
        end
        
        t = (a + tau1*(b - a) + b - tau2*(b - a))/2;
        
        phi = f(x + t*d);
        dphi = g(x + t*d)'*d;
        psi = phi - f(x) - sigma*t*(g(x)'*d);
        psia = f(x + a*d) - f(x) - sigma*a*(g(x)'*d);
        dpsi = dphi - sigma*(g(x)'*d);
        
        if psi >= psia
            b=t;
        elseif psi < psia && norm(dpsi) <= (rho - sigma)*norm(g(x)'*d)
            return
        elseif psi < psia && (t - a)*dpsi < 0
            a = t;
        elseif psi < psia && (t - a)*dpsi > 0
            b = a;
            a = t;
        end
        counter = counter + 1;
    end
end

