function  [theta_dot, psi_dot] = ddt_II(theta, psi,sigma, a,k,w,what,z)

alpha = (theta + z*psi)/2;
beta = z*(theta - a*theta./k)/2 - psi;

xalpha = a*exp(1i*alpha);
xbeta = a*exp(1i*beta);
theta_dot = w + sigma*imag(xalpha.*exp(-1i*alpha))./k;

psi_dot = what - sigma*imag(xbeta)./(2*k) + sigma*imag(exp(1i*beta))/2;

end
