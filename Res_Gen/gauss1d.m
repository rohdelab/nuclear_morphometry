function [y] = gauss1d(x,mu,sig)

b=(x-mu); b=b.^2;
c=2*sig*sig;

r=b/c;

y=exp(-r);
end