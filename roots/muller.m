function r = muller(y, x1, x2, x3, threshold)
    %{
        MULLER: Computes single solution to equation using Muller's method

        Parameters:
            y: Equation equal to zero
            x1: First guess
            x2: Second guess
            x3: Third guess
            threshold: Margin of error in function evaluation (how close to
            zero does the evaluation have to be to stop recursing)

        Returns:
            Root computed
    %}
    r = muller_rec(y, [x1; x2; x3], threshold);
end

function r = muller_rec(y, x, threshold)
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
    if abs(y(xnew)) > threshold
        r = muller_rec(y, [x(2:3); xnew], threshold); 
    else
        r = xnew;
    end
end