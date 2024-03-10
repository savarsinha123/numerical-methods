function m = rk4(f, initial, h, n)
    %{
        RK4: Implements Runge-Kutta method to compute solution
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

    % calculating slopes
    for i = 1:n
        k1 = f(values(1, i), values(2, i));
        k2 = f(values(1, i) + h/2, values(2, i) + h * k1/2);
        k3 = f(values(1, i) + h/2, values(2, i) + h * k2/2);
        k4 = f(values(1, i) + h, values(2, i) + h * k3);
        xnew = values(1, i) + h;
        ynew = values(2, i) + h/6 * (k1 + 2 * k2 + 2 * k3 + k4);
        values(:, i + 1) = [xnew; ynew];
    end
    m = values;
end