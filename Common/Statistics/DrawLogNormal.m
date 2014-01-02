%% Sampling from a bivariate Log-normal distribution

function [X] = DrawLogNormal(mu, sigma, N)

    % For convenience in the PDF evaluation, we use a reparameterized
    % version of the bivariate log-normal distribution!

    %sigma remains the same
    %mu is transformed
    mu = mu + diag(sigma) + [sigma(1, 2); sigma(2, 1)];
    
    K = length(mu);

    Y = lognrnd(0, 1, K, N);

    logY = log(Y);
    
    logX = sigma^(0.5) * logY + repmat(mu, 1, N);

    X = exp(logX);
end