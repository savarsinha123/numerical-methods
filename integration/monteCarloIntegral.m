function I = monteCarloIntegral(f, a, b, subset, n)
    %{
        monteCarloIntegral: This function evaluates the integral of f over
        the subset using the assumption that the function is enclosed
        within a rectangular region whose bounds are specified by a and b.
        A total of n points are sampled to calculate this integral.

        Parameters:
            f: Integrand
            a: Lower bounds of enclosing region
            b: Upper bounds of enclosing region
            subset: Region of integration inclosed within [a, b]
            n: Number of points sampled

        Returns:
            Estimated value of integral
    %}

    % sampling uniformly across [a, b] rectangular region
    N = nargin(f);
    points = a + (b - a) .* rand(n, N);
    point_args = mat2cell(points', ones(N, 1));

    % check how many points in domain and calculate volume
    mask = subset(point_args{:});
    new_points = points(mask, :);
    V = sum(mask)/n * prod(b - a);

    % evaluate average value and integral
    args = mat2cell(new_points', ones(N, 1));
    I = V*mean(f(args{:}));
end
