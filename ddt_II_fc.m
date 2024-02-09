function  [theta_dot, psi_dot] = ddt_II_fc(theta, psi,sigma,w,what,z)

N = length(theta);

alpha = (theta + z*psi)/2;
beta = z*(N/(N-1))*(theta - sum(theta)/N)/2 - psi;

xalpha = sum(exp(1i*alpha))/N;
xbeta = sum(exp(1i*beta))/N;

theta_dot = w + (N/(N-1))*sigma*imag(xalpha.*exp(-1i*alpha));
psi_dot = what - (N/(N-1))*sigma*imag(xbeta)/2 + (N/(N-1))*sigma*imag(exp(1i*beta))/2;

end
