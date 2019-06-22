function [sigmam] = strongPW(gamma, eta, d, f, g, maxsteps, x)
    counter = 0;
    sigmam = 1;
    a = (f(x + sigmam*d) <= (sigmam*gamma*(g(x)')*d));
    if a
        p = ((g(x+sigmam*d)')*d >= eta*(g(x)')*d);
        if p
            return
        else
            sigmap = 2;
            a = (f(x + sigmap*d) <= (sigmap*gamma*(g(x)')*d));
            while a
                counter = counter + 1;
                if (counter > maxsteps)
                    error("A fitting stepsize couldn't be found in " + maxsteps + " steps.")
                end
                sigmap = 2*sigmap;
                a = (f(x + sigmap*d) <= (sigmap*gamma*(g(x)')*d));
            end
            sigmam = sigmap/2;
        end
    else 
        while ~a
            counter = counter + 1;
            if (counter > maxsteps)
                error("A fitting stepsize couldn't be found in " + maxsteps + " steps.")
            end
            sigmam = sigmam/2;
            %f(x + sigmam*d)
            %(sigmam*gamma*(g(x)')*d)
            a = (f(x + sigmam*d) <= (sigmam*gamma*(g(x)')*d));
        end
        sigmap = 2*sigmam;
    end

    counter = 0;
    p = ((g(x+sigmam*d)')*d >= eta*(g(x)')*d);
    while ~p
        if (counter > maxsteps)
            error("A fitting stepsize couldn't be found in " + maxsteps + " steps.")
        end
        sigma = (sigmam + sigmap)/2;
        a = (f(x + sigma*d) <= (sigma*gamma*(g(x)')*d));
        if a
            sigmam = sigma;
        else
            sigmap = sigma;
        end
        p = ((g(x+sigmam*d)')*d >= eta*(g(x)')*d);
    end
end

