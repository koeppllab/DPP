function [rndIdx, p] = SampleParticleGivenParameterME(pDist, N, options)

    K = pDist.NumParticles;

    weights = zeros(K, 1);

    
    if (strcmp(options.Type, 'Rate'))
        
        APrior = options.ParameterPrior(:, 1);
        BPrior = options.ParameterPrior(:, 2);
        
        A = pDist.AChain+APrior(:, ones(1, K));
        B = pDist.BChain+BPrior(:, ones(1, K));
        

        
        APrior = option.ParameterPrior(1);
        BPrior = option.ParameterPrior(2);
        
        LogWeightsFromParameters = GammaLogPDF(options.Parameter, A+APrior, B+BPrior);
        LogWeightsFromParameters(and(A == 0, B == 0)) = 0;
%         
        weights = log(pDist.Weights) + sum(LogWeightsFromParameters, 1);
%         
%         
%         % random effect params
%         
        APrior = options.HyperParameter(1); %pDist.RandomEffectHyperChain(1, :);
        BPrior = options.HyperParameter(2); %pDist.RandomEffectHyperChain(2, :);
        
        for l=1:pDist.NumCells
            A = pDist.RandomEffectAChain(:, :, l) + APrior;
            B = pDist.RandomEffectBChain(:, :, l) + BPrior;
        
            randomEffect = options.RandomEffect(:, :, l);
            
            %LogPrior = GammaLogPDF(randomEffect, APrior, BPrior);
            
            LogWeightsFromRandomEffects = GammaLogPDF(randomEffect, ...
                A, B);
            LogWeightsFromRandomEffects(and(A == 0, B == 0)) = 0;
            
            weights = weights + sum(LogWeightsFromRandomEffects, 1);% + sum(LogPrior, 1);
        end
        
%         for l=1:pDist.NumCells
%            a = pDist.Paths{l, k}.a;
%            b = 
%         end
        
%         for k=1:K
%             particle = GetParticleSample(pDist, k);
%             A = particle.A;
%             B = particle.B;
%             
%             if (sum(A==0)==0 && sum(B==0)==0)
%                 weights(k) = log(pDist.Weights(k)) + sum(GammaPDF(options.Parameter, A, B));
%             else
%                weights(k) = log(pDist.Weights(k)); 
%             end
%         end

        wSum = sum(exp(weights));

        weights = exp(weights - log(wSum));

        rndIdx = DrawSample(N, weights');

        zeroIdx = rndIdx == 0;
        rdVec = unidrnd(K, 1, length(zeroIdx));
        rndIdx(zeroIdx) = rdVec;

        p = weights(rndIdx);
    else
        rndIdx = 0;
        p = 0;
    end

end