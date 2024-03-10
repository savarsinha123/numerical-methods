function m = backwardEuler(f, initial, h, n)
    %{
        BACKWARDEULER: Implements backward euler method to compute solution
        to IVP
        
        Parameters:
            f: implicit differential expression
            initial: initial conditions provided as column vector
            h: step size
            n: number of steps

        Returns:
            matrix containing discrete data points that approximate
            solution
    %}
    values = zeros(2, n);
    values(:, 1) = initial;
    for i = 1:n
        eval = f(values(1, i), values(2, i));
        xnew = values(1, i) + h;
        ystep = @(y1) y1 - h * f(xnew, y1) - eval;
        ynew = newtonRaphson(ystep, eval, 1e-05);
        values(:, i + 1) = [xnew; ynew];
    end
    m = values;
end

% same Newton-Raphson as specified in finding_zeros
function r = newtonRaphson(y, x, threshold)
    % calculate y'(x) numerically
    deriv = @(x) (y(x + 1e-6) - y(x)) / 1e-6;
    while abs(y(x)) > threshold
        x = x - y(x)/deriv(x);
    end
    r = x;
end

