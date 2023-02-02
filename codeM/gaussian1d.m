function res = gaussian1d(x,sigma,mu)

res = 1/(sigma*sqrt(2*pi))*exp( -(x-mu).^2/(2*sigma^2));