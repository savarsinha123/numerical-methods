function m = forwardEuler(f, initial, h, n)
    %{
        FORWARDEULER: Implements forward euler method to compute solution
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
    values = zeros(length(initial), n);
    values(:, 1) = initial;
    for i = 1:n
        args = num2cell(values(:, i));
        eval = f(args{:});
        slopes = [values(3:end, i); eval];
        next = [values(1, i) + h; values(2:end, i) + h * slopes];
        values(:, i + 1) = next;
    end
    m = values;
end