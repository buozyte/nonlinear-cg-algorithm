function [y] = gradros(x)
    if length(x) ~= 2
        error("Wrong dimension")
    else
        y1 = 2 * x(1) - 400 * x(1) * (x(2) - x(1)^2) - 2;
        y2 = 200 * (x(2) - x(1)^2);
        y = [y1; y2];
    end
end

