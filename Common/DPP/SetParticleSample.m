function [pDistNew] = SetParticleSample(pDist, particle, idx)

    pDist.IndexChain(idx) = particle.Index;
    pDist.RateChain(:, idx) = particle.Rates;
    pDist.PathLikelihoodChain(:, idx) = particle.PathLikelihood;
    pDist.MeasurementLikelihoodChain(:, idx) = particle.MeasurementLikelihood;
    pDist.MeasurementLikelihoodChainPrior(:, idx) = particle.MeasurementLikelihoodPrior;
    pDist.CurrentStates(:, idx, :) = particle.CurrentState;
    pDist.Weights(idx) = particle.Weight;
    pDist.Paths(:, idx) = particle.Paths;
    
    pDist.RandomEffectChain(:, idx, :) = particle.RatesRandomEffect;
    pDist.RandomEffectAChain(:, idx, :) = particle.RandomEffectA;
    pDist.RandomEffectBChain(:, idx, :) = particle.RandomEffectB;
    pDist.RandomEffectHyperChain(:, idx) = particle.RandomEffectHyper;
    pDist.MorphHyperChain(:, idx) = particle.MorphHyper;
    
    pDist.MeasurementVarianceChain(:, idx) = particle.MeasurementVariance;
    pDist.ErrorFunctionChain(:, idx) = particle.ErrorFunction;
    
    
    pDistNew = pDist;
    
end