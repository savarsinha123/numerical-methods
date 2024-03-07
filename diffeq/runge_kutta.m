% solving y' = y sin^2(x) using RK4
y1 = @(x, y0) sin(x)^2 * y0;
sol = rk4(y1, [0; 1], 0.01, 1000);

figure
plot(sol(1, :), sol(2, :))

% function declarations
function m = rk4(y, initial, h, n)
    values = zeros(2, n);
    values(:, 1) = initial;
    for i = 1:n
        k1 = y(values(1, i), values(2, i));
        k2 = y(values(1, i) + h/2, values(2, i) + h * k1/2);
        k3 = y(values(1, i) + h/2, values(2, i) + h * k2/2);
        k4 = y(values(1, i) + h, values(2, i) + h * k3);
        xnew = values(1, i) + h;
        ynew = values(2, i) + h/6 * (k1 + 2 * k2 + 2 * k3 + k4);
        values(:, i + 1) = [xnew; ynew];
    end
    m = values;
end