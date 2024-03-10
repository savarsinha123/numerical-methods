function r = newtonRaphson(y, x, threshold)
    %{
        MULLER: Computes single solution to equation using Newton-Raphson

        Parameters:
            y: Equation equal to zero
            x: Initial guess
            threshold: Margin of error in function evaluation (how close to
            zero does the evaluation have to be to stop iterating)

        Returns:
            Root computed
    %}

    % calculate y'(x) numerically
    deriv = @(x) (y(x + 1e-6) - y(x)) / 1e-6;
    while abs(y(x)) > threshold
        x = x - y(x)/deriv(x);
    end
    r = x;
end
