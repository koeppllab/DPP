function pPosteriorFinalized = FinalizeParticleDistribution(pPosterior, currentTime, burnIn)

    pPosteriorFinalized = pPosterior;
    pPosteriorFinalized.CurrentTime = currentTime;
    
    K = size(pPosterior.Weights, 2);
    
    startIndex = burnIn + 1;
    
    pPosteriorFinalized.IndexChain = pPosterior.IndexChain(startIndex:end);
    pPosteriorFinalized.PathLikelihoodChain = pPosterior.PathLikelihoodChain(:, startIndex:end);
    pPosteriorFinalized.MeasurementLikelihoodChain = pPosterior.MeasurementLikelihoodChain(:, startIndex:end);
    pPosteriorFinalized.MeasurementLikelihoodChainPrior = pPosterior.MeasurementLikelihoodChainPrior(:, startIndex:end);
    pPosteriorFinalized.Weights = ones(1, K - burnIn) / (K - burnIn);
    pPosteriorFinalized.RandomEffectHyperChain = pPosterior.RandomEffectHyperChain(:, startIndex:end);
    pPosteriorFinalized.MorphHyperChain = pPosterior.MorphHyperChain(:, startIndex:end);
    pPosteriorFinalized.CurrentStates = pPosterior.CurrentStates(:, startIndex:end, :);
    
    pPosteriorFinalized.NumParticles = size(pPosteriorFinalized.Weights, 2);

    pPosteriorFinalized.Paths = pPosterior.Paths(:, startIndex:end);
    
    pPosteriorFinalized.MeasurementVarianceChain = ...
        pPosterior.MeasurementVarianceChain(:, startIndex:end);
   
    pPosteriorFinalized.ErrorFunctionChain = ...
        pPosterior.MeasurementVarianceChain(:, startIndex:end);
end