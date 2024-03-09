function I = gaussQuad(y, a, b, n)
    %{
        gaussQuad: This function evalues the integral of y over the bounds 
        defined by a and b using n total nodes
        
        Parameters:
            y: integrand
            a: lower bounds
            b: upper bounds
            n: number of nodes

        Returns:
            I: integral value
    %}

    % calculating weights
    P = legendreP2(n);
    roots = find_zeros(P, n);
    Pderiv = @(x) (P(x + 1e-6) - P(x)) ./ 1e-6;
    weights = 2./((1 - roots.^2).*(Pderiv(roots)).^2);

    % combine weights in multidimensional case
    N = nargin(y);
    all_weights = 1;
    for i = 1:N
        all_weights = kron(all_weights, weights);
    end

    % compute modified inputs for each function in order to change bound
    % to [-1, 1] from [a, b]
    x = cell(1, N);
    for i = 1:N
        transform = @(x) 0.5 * (b(i) - a(i)) * x + 0.5 * (b(i) + a(i));
        x{i} = transform(roots);
    end
    ntuples = mat2cell(cartesian(x{:})', ones(N, 1));

    % evaluate integral
    diffs = 0.5*(b - a);
    I = prod(diffs) .* sum(all_weights .* y(ntuples{:})');
end

function C = cartesian(varargin)
    % calculates Cartesian product
    args = varargin;
    n = nargin;

    [F{1:n}] = ndgrid(args{:});

    for i=n:-1:1
        G(:,i) = F{i}(:);
    end

    C = G;
end

function P = legendreP2(n)
    % recursively calcluate Legendre polynomials
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