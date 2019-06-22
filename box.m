function [y] = box(x)
    if (length(x) ~= 3)
        error("Wrong dimension")
    else
        y = 0;
        x1 = x(1);
        x2 = x(2);
        x3 = x(3);
        for i = 1:3
            t = 0.1 * i;
            y = y + (exp(-t * x1) - exp(-t * x2) - ...
                                        x3 * (exp(-t) - exp(-10 * t)))^2;
        end    
    end
end

