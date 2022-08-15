% Finding the wavelength, given the value of d, T

function [L] = get_wavelength(d,T)
g = 9.81;
L = 1.56 * T^2;
e = 10;
while abs(e) > 1e-05
    e = L - (g*T^2/2/pi)*tanh(2*pi*d/L);
    L = L - e;
end
end