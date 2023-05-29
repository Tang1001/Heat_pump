function [T, Y] = ode4(fun, tspan, y0, h)
    T = tspan(1):h:tspan(2);
    n = length(T);
    Y = zeros(n, length(y0));
    Y(1,:) = y0';

    for i = 1:(n-1)
        k1 = fun(T(i), Y(i,:))';
        k2 = fun(T(i)+h/2, Y(i,:)+h*k1/2)';
        k3 = fun(T(i)+h/2, Y(i,:)+h*k2/2)';
        k4 = fun(T(i)+h, Y(i,:)+h*k3)';
        Y(i+1,:) = Y(i,:) + h/6 * (k1 + 2*k2 + 2*k3 + k4);
    end
end