function [x] = lognrnd_fast(mu, sigma)
    z = randn;
    x = exp(mu + sigma*z);
end