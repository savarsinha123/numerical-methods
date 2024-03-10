function x = idft(X)
    %{
        IDFT: Computes inverse discrete Fourier transform of dataset

        Parameters:
            x: Input frequency spectrum

        Returns:
            Time series corresponding to Fourier transform, where the
            "time" of each point is determined based on the resolution of
            the frequency spectrum (step size between frequencies).
    %}

    % define length of vector
    N = length(X);
    
    % define wiggle factor
    W = exp(2j * pi / N);
    
    % construct inverse DFT matrix
    [R, C] = meshgrid(1:N, 1:N);
    A = (1 / sqrt(N)) * W .^ ((R - 1) .* (C - 1));
    
    % perform matrix multiplication
    x = A * X;
end

