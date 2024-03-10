function X = dft(x)
    %{
        DFT: Computes discrete Fourier transform of dataset

        Parameters:
            x: Input time series/dataset

        Returns:
            Fourier amplitudes corresponding to each frequency, over a
            range of values corresponding to Nth fractions of the sampling
            rate.
    %}

    % define length of vector
    N = length(x);
    
    % define wiggle factor
    W = exp(-2j * pi / N);
    
    % construct DFT matrix
    [R, C] = meshgrid(1:N, 1:N);
    A = (1 / sqrt(N)) * W .^ ((R - 1) .* (C - 1));
    
    % perform matrix multiplication
    X = A * x;
end

