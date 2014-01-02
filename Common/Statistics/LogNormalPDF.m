function propDens = LogNormalPDF(mu, sigma, x)

    K = length(mu);

    if (det(sigma) <= 0)
        propDens = 0;
    elseif (sum(x<=0) > 0)
        propDens = 0;
    else
        logx = log(x);
        fastDet = sigma(1,1)*sigma(2,2) - sigma(1,2)*sigma(2,1);
        propDens = 1 /( (2*pi)^(K/2) * sqrt(fastDet) ) * ...
            exp(-0.5 * (logx - mu)' * (sigma \ (logx - mu)));
    end

end