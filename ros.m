function [y] = ros(x)
    if length(x) ~= 2
        error("Wrong dimension")
    else
        y = 100 * (x(2) - x(1)^2)^2 + (1 - x(1))^2;
    end
end

