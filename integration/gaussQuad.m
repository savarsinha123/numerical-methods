function I = gaussQuad(y, a, b, n)
    ynew = @(x) 0.5*(b-a)*y(0.5*(b-a)*x + 0.5*(a + b));
    P = legendreP2(n);
    roots = find_zeros(P, n);
    Pderiv = @(x) (P(x + 1e-6) - P(x)) ./ 1e-6;
    weights = 2./((1 - roots.^2).*(Pderiv(roots)).^2);
    I = sum(weights .* arrayfun(ynew, roots));
end

function P = legendreP2(n)
    if n == 0
        P = @(x) 1;
    elseif n == 1
        P = @(x) x;
    else
        P1 = legendreP2(n - 1);
        P2 = legendreP2(n - 2);
        P = @(x) ((2.*n-1).*x.*P1(x) - (n-1).*P2(x))./n;
    end
end

function r = newtonRaphson(y, x, threshold)
    % calculate y'(x) numerically
    deriv = @(x) (y(x + 1e-6) - y(x)) / 1e-6;
    while abs(y(x)) > threshold
        x = x - y(x)/deriv(x);
    end
    r = x;
end

function r = find_zeros(y, n)
   r = find_zeros_rec(y, n, []);
end

function r = find_zeros_rec(y, n, l)
    % recursively find all roots
    if n == 0
        r = l;
    else
        root = newtonRaphson(y, 0.5, 1e-8);
        y1 = @(x) y(x)/(x - root);
        l = [l; root];
        r = find_zeros_rec(y1, n - 1, l);
    end
end