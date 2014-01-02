function propDens = NormalPDF(mu, sigma, x)

    K = length(mu);

    if (det(sigma) <= 0)
        propDens = 0;
    else
        propDens = 1 /( (2*pi)^(K/2) * sqrt(det(sigma)) ) * ...
            exp(-0.5 * (x - mu)' * inv(sigma) * (x - mu));
    end

end