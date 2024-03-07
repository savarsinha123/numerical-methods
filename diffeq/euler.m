% example calculation with y' = y - x, y(0) = 2
y = @(x, y0) y0 - x;
sol1 = forwardEuler(y, [0; 2], 0.01, 1000);

% plotting against actual solution y = e^x + x + 1
figure
plot(sol1(1, :), sol1(2, :))
hold on
x1 = linspace(0, 10, 1000);
y1 = exp(x1) + x1 + 1;
plot(x1, y1)

hold off

% example calculation with z'' = -z, z(0) = 1, z'(0) = 0
z = @(x, y0, y1) -y0;
sol2 = forwardEuler(z, [0; 1; 0], 0.01, 5000);

% plotting against actual solution y = cos(x)
figure
plot(sol2(1, :), sol2(2, :))
hold on
x2 = linspace(0, 50, 1000);
y2 = cos(x2);
plot(x2, y2)

hold off

% backward euler for y' = y - x, y(10) = e^10 + 10 + 1
sol3 = backwardsEuler(y, [10; exp(10) + 10 + 1], 0.01, 1000);
figure
plot(sol3(1, :), sol3(2, :))
hold on
x3 = linspace(10, 20, 1000);
y3 = exp(x3) + x3 + 1;
plot(x3, y3)

% forward euler
function m = forwardEuler(y, initial, h, n)
    values = zeros(length(initial), n);
    values(:, 1) = initial;
    for i = 1:n
        args = num2cell(values(:, i));
        eval = y(args{:});
        slopes = [values(3:end, i); eval];
        next = [values(1, i) + h; values(2:end, i) + h * slopes];
        values(:, i + 1) = next;
    end
    m = values;
end

% backward euler
function m = backwardsEuler(y, initial, h, n)
    values = zeros(2, n);
    values(:, 1) = initial;
    for i = 1:n
        eval = y(values(1, i), values(2, i));
        xnew = values(1, i) + h;
        ystep = @(y1) y1 - h * y(xnew, y1) - eval;
        ynew = newtonRaphson(ystep, eval, 1e-05);
        values(:, i + 1) = [xnew; ynew];
    end
    m = values;
end

% same Newton-Raphson as specified in findingzeros.m
function r = newtonRaphson(y, x, threshold)
    % calculate y'(x) numerically
    deriv = @(x) (y(x + 1e-6) - y(x)) / 1e-6;
    while abs(y(x)) > threshold
        x = x - y(x)/deriv(x);
    end
    r = x;
end