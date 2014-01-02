function [X] = DrawNormal(mu, sigma, N)

    K = length(mu);

    Y = normrnd(0, 1, K, N);
    
    X = sigma^(0.5) * Y + repmat(mu, 1, N);
   
end