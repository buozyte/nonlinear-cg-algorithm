function [y] = gradbox(x)
    if (length(x) ~= 3)
        error("Wrong dimension")
    else
        y1 = 0;
        y2 = 0;
        y3 = 0;
        x1 = x(1);
        x2 = x(2);
        for i = 1:3
            t = 0.1 * i;
            y1 = y1 - t * exp(-t * x1);
            y2 = y2 + t * exp(-t * x2);
            y3 = y3 - (exp(-t) - exp(-10 *t));
        end
        y = [y1; y2; y3];
    end
end

