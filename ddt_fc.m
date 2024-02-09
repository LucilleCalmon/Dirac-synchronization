function  [theta_dot, psi_dot] = ddt_fc(theta,psi,sigma,w,what)
N = length(w);

alpha = theta + psi/2;
beta = (N/(N-1))*(theta-sum(theta)/N) - psi;  

xalpha = sum(exp(1i*alpha))/N;
xbeta = sum(exp(1i*beta))/N;

theta_dot = w + sigma*imag(xalpha*exp(-1i*alpha));
psi_dot = what - sigma*imag(xbeta) - sigma*imag(exp(-1i*beta));

end
