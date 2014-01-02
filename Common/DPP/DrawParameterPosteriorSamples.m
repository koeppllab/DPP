function [C, HyperParameters, MorphParameters] = DrawParameterPosteriorSamples(pDist, ...
    aPrior, bPrior, idx, N, cellIdx)


    if (nargin < 6)
        cellIdx = 1:size(pDist.Paths, 1);
    end
    
    if (N<1)
       delta = floor(1 / N); 
       M = pDist.NumParticles;
    else
       delta = 1; 
       M = floor(N*pDist.NumParticles);
    end
    
    K = floor(N*pDist.NumParticles);
          
    C = zeros(length(idx), K);
    HyperParameters = zeros(2, K);
    MorphParameters = zeros(2, K);
     
    
    for i=1:delta:M
       
        particle = GetParticleSample(pDist, mod(i-1, pDist.NumParticles)+1);%particles{i};
        

        index = floor((i-1)/delta) + 1;
        if (~isempty(idx))
            aData = zeros(length(idx), 1);
            bData = zeros(length(idx), 1);
            for k=1:length(particle.Paths)
                if (sum(cellIdx == k)>0)
                    aData = aData + particle.Paths{k}.a(idx);
                    bData = bData + particle.Paths{k}.b(idx);
                end
            end
            
            C(:, index) = SampleParameterConditional(aData + aPrior,...
                bData + bPrior);
        else
            C = [];
        end
        
        HyperParameters(:, index) = particle.RandomEffectHyper;
        MorphParameters(:, index) = particle.MorphHyper;
        
    end

end