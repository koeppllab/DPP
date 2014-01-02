function [particle] = GetParticleSample(pDist, idx)

    particle.Index = pDist.IndexChain(idx);

    particle.PathLikelihood = pDist.PathLikelihoodChain(:, idx);
    particle.MeasurementLikelihood = pDist.MeasurementLikelihoodChain(:, idx);
    particle.MeasurementLikelihoodPrior = pDist.MeasurementLikelihoodChainPrior(:, idx);
    particle.Weight = pDist.Weights(idx);

    particle.CurrentState = pDist.CurrentStates(:, idx, :);
    particle.Paths = pDist.Paths(:, idx);

    particle.RandomEffectHyper = pDist.RandomEffectHyperChain(:, idx, :);
    particle.MorphHyper = pDist.MorphHyperChain(:, idx, :);
    
    particle.MeasurementVariance = pDist.MeasurementVarianceChain(idx);
    
    particle.ErrorFunction = pDist.ErrorFunctionChain(idx);
end