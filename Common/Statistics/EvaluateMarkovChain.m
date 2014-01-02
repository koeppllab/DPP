function [h, bins, mmseEstimate, mapEstimate] = EvaluateMarkovChain(chain, startIndex, numBins)


    values = chain(startIndex:end);
    n = length(values);
  
    % compute mmse
    mmseEstimate = mean(values);

    % compute sample std. dev.
    sigma = std(values);
    
    if (nargin < 3)
        
        % optimal bin width (Scott's choice, see Wikipedia: Histogram)
        dOpt = 3.5*sigma / (n)^(1/3);
        
        numBins = (max(values) - min(values)) / dOpt;
    end
    
    
    % compute histogram;
    [h, bins] = hist(values, numBins);
    deltaT = (max(bins) - min(bins))/numBins;
    h = h/sum(h*deltaT);
    
    % compute MAP estimate
    [~, maxIndex] = max(h);
    mapEstimate = bins(maxIndex);
end