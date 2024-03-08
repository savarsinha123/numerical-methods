% sine wave
samples = 1000;
t = linspace(0, 2*pi, samples)';
x1 = sin(10*t);
X1 = dft(x1);
w = fourierFreqs(t);
figure
plot(w(1:samples/2), (abs(X1(1:samples/2))))

% square wave
x2 = square(10 * t);
X2 = dft(x2);
figure
plot(w(1:samples/2), (abs(X2(1:samples/2))))
hold on
plot(w(samples/100:samples/2), 20.1331./(0.1 * w(samples/100:samples/2)))
hold off

% recovering square wave from FT using inverse DFT
x2i = idft(X2);
figure
plot(t(1:samples), real(x2i(1:samples)))

% comparing accuracy of Cooley-Tukey FFT with DFT when size is not
% power of 2
x3 = sawtooth(10 * t);
X31 = ctfft(x3);
X32 = dft(x3);
figure
hold on
plot(w(1:samples/2), (abs(X31(1:samples/2))))
plot(w(1:samples/2), (abs(X32(1:samples/2))))
hold off


function X = dft(x)
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

function x = idft(X)
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

function w = fourierFreqs(t)
    sampling_rate = 1/(t(2) - t(1));
    N = length(t);
    w = 2*pi*sampling_rate*(0:N-1)/N;
end

function X = ctfft(x)
    N = length(x);
    new_size = 2^nextpow2(N);
    xnew = [zeros(new_size - N, 1); x];
    X = ctfft_rec(xnew, new_size, 2) / sqrt(new_size);
end

function X = ctfft_rec(x, N, s)
    if N == 1
        X = x;
    else
        % split the input sequence into even and odd indices
        x_even = x(1:s:end);
        x_odd = x((1 + s/2):s:end);

        % recursive DFT computation for both halves
        X_even = ctfft_rec(x_even, N/2, 2);
        X_odd = ctfft_rec(x_odd, N/2, 2);

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