% function to find the zeros of
y = @(x) x^3-3*x-10;

% examples with Newton-Raphson and Muller
disp("The root found by Newton-Raphson is x = " + ...
    newtonRaphson(y, 1.1, eps));
disp("The root found by Muller's method is x = " + ...
    muller(y, -1, -2, -3));
disp("The 3 roots of y are given by")
disp((find_zeros(y, 3)))

% example of finding some of the roots for a transcendental expression
z = @(x) x * exp(x) - 10 * x^2;
disp("5 roots of the function z are")
disp((find_zeros(z, 5)))

% function declarations
function r = newtonRaphson(y, x, threshold)
    % calculate y'(x) numerically
    deriv = @(x) (y(x + 1e-6) - y(x)) / 1e-6;
    while abs(y(x)) > threshold
        x = x - y(x)/deriv(x);
    end
    r = x;
end

function r = muller(y, x1, x2, x3)
    r = muller_rec(y, [x1; x2; x3]);
end

function r = muller_rec(y, x)
    % fit quadratic through three guesses
    xmatrix = [(x - x(3)).^2, (x - x(3)), ones(3, 1)];
    yvalues = arrayfun(y, x);
    coeff = xmatrix \ yvalues;
    discrim = sqrt(coeff(2)^2 - 4 * coeff(1) * coeff(3));

    % choose sign for denominator that maximizes its magnitude
    if abs(coeff(2) + discrim) > abs(coeff(2) - discrim)
       xnew = x(3) - 2 * coeff(3) / (coeff(2) + discrim);
    else
       xnew = x(3) - 2 * coeff(3) / (coeff(2) - discrim);
    end

    % if the accuracy threshold is met, return, else recurse with new guess
    if abs(y(xnew)) > 1e-5
        r = muller_rec(y, [x(2:3); xnew]); 
    else
        r = xnew;
    end
end

function r = find_zeros(y, n)
   r = find_zeros_rec(y, n, []);
end

function r = find_zeros_rec(y, n, l)
    % recursively find all roots
    if n == 0
        r = l;
    else
        root = muller(y, -1.2382394234, 1.23042234, 0.09234234);
        y1 = @(x) y(x)/(x - root);
        l = [l; root];
        r = find_zeros_rec(y1, n - 1, l);
    end
end