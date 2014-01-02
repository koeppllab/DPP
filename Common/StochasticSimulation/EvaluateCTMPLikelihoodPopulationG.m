function [L] = EvaluateCTMPLikelihoodPopulationG(R, GInt, populationRateIdx,...
                                                    a, b, u, v, k)

    r = R(populationRateIdx);
    gInt = GInt(populationRateIdx);

    if (nargin == 5)
       u = 0;
       v = 0;
       k = 0;
       
       morphFeatures = 0;
    else
       morphFeatures = 1; 
    end
    
    if (morphFeatures == 1)
    	L = sum(a.*log(b) + gammaln(a+r+u) - gammaln(a) - (a+r+u).*log(b+gInt+k*v) ...
                  + (u-1).*log(k) - gammaln(u) + u.*log(v));
    else
        L = sum(a.*log(b) + gammaln(a+r) - gammaln(a) - (a+r).*log(b + gInt));
    end

end