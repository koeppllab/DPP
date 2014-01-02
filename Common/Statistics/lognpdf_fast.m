function [p, erf] = lognpdf_fast(x, mu, sigma)


    erf = power(log(x) - mu, 2); 
    p = 1/(x*sigma*sqrt(2*pi))*exp(-erf / (2*sigma^2));

end