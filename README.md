# numerical-methods

This repository contains a collection of implementations of various
numerical methods such as Newton-Raphson, backward Euler, and discrete
Fourier transform.

## Roots

###  Newton-Raphson

This root-finding method iteratively moves towards the root of a function 
starting from an initial guess by using the derivative of the function. If the
initial guess is real and the function itself maps reals to reals, then the 
roots found will only be restricted to the real domain. See here for more info:
https://en.wikipedia.org/wiki/Newton%27s_method

### Muller's Method

This root-finding algorithm locates the roots of a function by interpolating
a quadratic through three points of the given function and recursively
approximating the root of the function using the roots of this quadratic. Since
a quadratic function is used, this method can also find complex roots. See here 
for more info:
https://en.wikipedia.org/wiki/Muller%27s_method

## Differential Equations

### Forward Euler

Forward Euler uses linear approximations of the solution to the differential
equation at each point to predict the next point at the next step. This is a
method of order 1. See here for more info:
https://en.wikipedia.org/wiki/Euler_method

### Backward Euler

Backward Euler uses the evaluation of y' at the next step to compute the next
value of y. This is a method of order 1. See here for more info:
https://en.wikipedia.org/wiki/Backward_Euler_method 

### Runge-Kutta (RK4)

Runge-Kutta, in partiular the RK4 method, utilizes four different slopes to
predict the next value of y at the next step. This is a fourth order method. See
here for more info:
https://en.wikipedia.org/wiki/Runge%E2%80%93Kutta_methods

## Fourier Transform

### Discrete Fourier Transform (DFT)

DFT represents a discretized version of the Fourier integral and is often
encoded in a unitary matrix. See here for more info:
https://en.wikipedia.org/wiki/Discrete_Fourier_transform
https://en.wikipedia.org/wiki/DFT_matrix

### Fast Fourier Transform (FFT)

FFT refers to any methods of calculating DFT that have complexity $O(n log(n))$
rather than $O(n^2)$, as is typically standard for simpler DFT methods. The 
discrete_fourier_transform.m file contains an implementation of the Cooley-Tukey
Algorithm. See here for more info:
https://en.wikipedia.org/wiki/Fast_Fourier_transform
https://en.wikipedia.org/wiki/Cooley%E2%80%93Tukey_FFT_algorithm

## Integration

### Gauss-Legendre Quadrature

Gauss-Legendre Quadrature approximates an integral over the interval $[-1, 1]$ 
by approximating the integral as a weighted sum of function evaluations of the
nth Legendre polynomial roots, where larger $n$ results in more accurate
approximations. See here for more info:
https://en.wikipedia.org/wiki/Gauss%E2%80%93Legendre_quadrature

### Monte-Carlo Integration

Monte-Carlo integration refers to techniques which involve randomly sampling
from the integration domain, computing the MLE for the expected value of the
function, and then multiplying this mean by the volume of the subset. See here 
for more info:
https://en.wikipedia.org/wiki/Monte_Carlo_integration