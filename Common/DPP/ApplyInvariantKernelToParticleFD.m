function [newParticle] = ApplyInvariantKernelToParticleFD(particle, model, cells, N)

numCells = length(cells);

Paths = particle.Paths;

randomEffectIdx = model.RandomEffect == 1;
randomEffectKineticIdx = model.RandomEffectIndex(randomEffectIdx);

% for k=1:numCells
%     particle.Paths{k}.a = lognrnd(log(particle.Paths{k}.a), 0.02);
%     particle.Paths{k}.b = lognrnd(log(particle.Paths{k}.b), 0.02);
% end


if (~isempty(randomEffectKineticIdx))
    
    PathsTmp = Paths;
    
    for j=1:N
        sigmaProp = 0.3; %hardcoded. make variable in future version
        
        hyperParams(1) = particle.RandomEffectHyper(1);
        hyperParams(2) = particle.RandomEffectHyper(2);
        
        if (model.MorphologicalFeatures == 1)
            hyperParams(3) = particle.MorphHyper(1);
            hyperParams(4) = particle.MorphHyper(2);
        end
        
        proposal = zeros(size(hyperParams));
        
        forwardProb = 1;
        backwardProb = 1;
        
        for m=1:length(hyperParams)
            
            logMu = log(hyperParams(m));
            proposal(m) = lognrnd_fast(logMu, sigmaProp);
            
            forwardProb = forwardProb * lognpdf_fast(proposal(m), log(hyperParams(m)), sigmaProp);
            backwardProb = backwardProb * lognpdf_fast(hyperParams(m), log(proposal(m)), sigmaProp);
            
        end
        
        
        
        LNewV = zeros(1, numCells);
        LOldV = zeros(1, numCells);
        
        for k=1:numCells
            
            if (sum(PathsTmp{k}.a == 0) == model.NumReactions)
                continue;
            end
            
            if (model.MorphologicalFeatures == 1)
                LNewV(k) = ...
                    EvaluateHeterogeneousPathLikelihood(PathsTmp{k}.a, ...
                    PathsTmp{k}.b, ...
                    randomEffectKineticIdx, ...
                    proposal(1),...
                    proposal(2),...
                    proposal(3),...
                    proposal(4),...
                    1 / cells{k}.MorphologicalFeatureValue);
                LOldV(k) = ...
                    EvaluateHeterogeneousPathLikelihood(PathsTmp{k}.a, ...
                    PathsTmp{k}.b, ...
                    randomEffectKineticIdx, ...
                    hyperParams(1), ...
                    hyperParams(2),...
                    hyperParams(3),...
                    hyperParams(4),...
                    1 / cells{k}.MorphologicalFeatureValue);
            else
                LNewV(k) = ...
                    EvaluateHeterogeneousPathLikelihood(PathsTmp{k}.a, ...
                    PathsTmp{k}.b, ...
                    randomEffectKineticIdx, ...
                    proposal(1),...
                    proposal(2));
                LOldV(k) = ...
                    EvaluateHeterogeneousPathLikelihood(PathsTmp{k}.a, ...
                    PathsTmp{k}.b, ...
                    randomEffectKineticIdx, ...
                    hyperParams(1), ...
                    hyperParams(2));
            end
        end
        
        LNewHyper = sum(LNewV);
        LOldHyper = sum(LOldV);
        
        if (sum(sum(model.HyperPriorSigma))>0)
            
            Mu = model.HyperPriorMu;
            Sigma = model.HyperPriorSigma;
            
            priorNew = LogNormalPDF(Mu, Sigma, ...
                [proposal(1); proposal(2)]);
            priorOld = LogNormalPDF(Mu, Sigma, ...
                [hyperParams(1); hyperParams(2)]);
            
        else
            priorNew = 1;
            priorOld = 1;
        end
        
        if (model.MorphologicalFeatures == 1 ...
                && sum(sum(model.MorphPriorSigma)) > 0)
            Mu = model.MorphPriorMu;
            Sigma = model.MorphPriorSigma;
            
            priorNew = priorNew * ...
                LogNormalPDF(Mu, Sigma, [proposal(3); proposal(4)]);
            priorOld = priorOld * ...
                LogNormalPDF(Mu, Sigma, [hyperParams(3); hyperParams(4)]);
            
        end
        
        
        accept = min(1, exp(LNewHyper - LOldHyper) * backwardProb / forwardProb ...
            * priorNew / priorOld);
        
        
        
        if (rand <= accept)
            %fprintf('A');
            
            hyperParams = proposal;
            
        else
            %fprintf('R');
            
        end
        
        
        particle.RandomEffectHyper(1) = hyperParams(1);
        particle.RandomEffectHyper(2) = hyperParams(2);
        
        if (model.MorphologicalFeatures == 1)
            particle.MorphHyper(1) = hyperParams(3);
            particle.MorphHyper(2) = hyperParams(4);
        end
        
    end
    
end

%% Estimate measurement noise parameters

if (model.EstimateMeasurementNoise == 1)
    
    %ErrorFunction = zeros(1, numCells);
    
    ErrorFunction = [];
    for k=1:numCells
        ErrorFunction = [ErrorFunction; Paths{k}.Erf(:)];
    end
    
    %ErrorFunction = Paths{1}.Erf;
    
    n = length(ErrorFunction);
    NoiseA = model.NoiseA + n / 2;
    NoiseB = model.NoiseB + sum(ErrorFunction) / 2;
    
    particle.MeasurementVariance = sqrt(1 / gamrnd(NoiseA, 1 / NoiseB));
    
end



newParticle = particle;


end