function X = ctfft(x)
    %{
        DFT: Computes fast Fourier transform of dataset via the
        Cooley-Tukey algorithm.

        Parameters:
            x: Input time series/dataset

        Returns:
            Fourier amplitudes corresponding to each frequency, over a
            range of values corresponding to Nth fractions of the sampling
            rate.
    %}

    N = length(x);
    new_size = 2^nextpow2(N);
    xnew = [zeros(new_size - N, 1); x];
    X = ctfftRec(xnew, new_size, 2) / sqrt(new_size);
end

function X = ctfftRec(x, N, s)
    if N == 1
        X = x;
    else
        % split the input sequence into even and odd indices
        x_even = x(1:s:end);
        x_odd = x((1 + s/2):s:end);
    
        % recursive DFT computation for both halves
        X_even = ctfftRec(x_even, N/2, 2);
        X_odd = ctfftRec(x_odd, N/2, 2);
    
        % combine DFTs of two halves into full DFT
        for k = 0:N/2-1
            twiddle = exp(-2j * pi / N * k);
            X_k = X_even(k + 1);
            X_k_N2 = twiddle * X_odd(k + 1);
            X(k + 1) = X_k + X_k_N2;
            X(k + 1 + N/2) = X_k - X_k_N2;
        end
    end
end