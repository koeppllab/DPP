function [particle, particleIdx] = DrawFromParticleDistribution(pDist, K)

    numParticles = length(pDist.NumParticles);
    
    weights = pDist.Weights / sum(pDist.Weights);
    particleIdx = DrawSample(K, weights');
    
    particle = cell(K, 1);
    
    for l=1:K
        idx = particleIdx(l);
        particle{l} = GetParticleSample(pDist, idx);
    end

end