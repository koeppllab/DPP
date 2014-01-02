function [pPosterior] = UpdatePosteriorSingleCell(pDist, model, targetCells, targetCellIdx, options)

callbackProvided = ~isempty(options.CallBackIteration);



if (length(options.M)==length(targetCells))
    M = options.M(targetCellIdx);
else
    M = options.M(1);
end

if (options.TimeIndex > 1)
   M = ceil(M * options.ParticleReduce); 
end

targetCell = targetCells{targetCellIdx};

maxTime = max(targetCell.MeasurementTime(max(options.TimeIndex)));

numCells = length(targetCells);

pPosterior = InitializeParticleDistributionME(M, model, 1);
currentParticle = GetParticleSample(pDist, 1);
acceptedCount = 0;

Paths = cell(1);
c = model.c;


allSimulatedCells = cell(1, M);

aPrior = model.RandomEffectA;
bPrior = model.RandomEffectB;

measurementVariance = model.MeasurementVariance;

randomEffectIdx = model.RandomEffect == 1;
randomEffectKineticIdx = model.RandomEffectIndex(randomEffectIdx);

sampleGrid = linspace(pDist.CurrentTime, ...
    model.MeasurementTime(max(options.TimeIndex)), options.StorePathNumPoints);

blockSize = options.BlockSize;
numBlocks = ceil(M / blockSize);
rest = mod(M, blockSize);
currentBlockSize = blockSize;




oldPSize = pDist.NumParticles;

if (M < oldPSize)
   delta = floor(oldPSize / M); 
   K = oldPSize;
else
   delta = 1; 
   K = M;
end

particleIdx = mod(0:delta:K, oldPSize) + 1;




%if (oldPSize < M)
%   particleIdx = [particleIdx, randi(oldPSize, 1, M-oldPSize)];
%end

timeOffset = pDist.CurrentTime;



for u=1:numBlocks
    
    if (u==numBlocks && rest > 0)
        currentBlockSize = rest;
    end
    
    if (callbackProvided == 1)
        options.CallBackIteration();
    end
    
    [particles, particleIdx] = DrawFromParticleDistribution(pDist, currentBlockSize);
    
    %particles = CopyParticles(pDist, particleIdx, currentBlockSize);
    

    %% Pre-simulate particles
    parfor l=1:currentBlockSize

        currParticleIdx = (u - 1) * blockSize + l;
        particleIdxTmp = particleIdx;
        
        optionsTmp = options;
        modelTmp = model;
        
        particleTmp = particles{l};
        prevIdx = particleIdx(l); %particleIdx(l);
        
        %% Resample static parameters.
        particleTmp = ApplyInvariantKernelToParticle(particleTmp, ...
            modelTmp, targetCells, optionsTmp.SubMCMCRuns, targetCellIdx);
        

        
        %% Resample Paths
        
        aPriorTmp = aPrior;
        bPriorTmp = bPrior;
        
        
        aPriorTmp(randomEffectIdx) = particleTmp.RandomEffectHyper(1);
        bPriorTmp(randomEffectIdx) = particleTmp.RandomEffectHyper(2);
        
        simulatedCell = ForwardSimulateSingleCell(particleTmp, model, 1, c, ...
            aPriorTmp, bPriorTmp, timeOffset, maxTime, options, targetCellIdx, targetCell);
        
        allSimulatedCells{l}.cell = simulatedCell;
        allSimulatedCells{l}.particle = particleTmp;
        allSimulatedCells{l}.Index = prevIdx;


        
    end
    
    for l=1:currentBlockSize
        currParticleIdx = (u - 1) * blockSize + l;
        
        MeasurementLikelihood = currentParticle.MeasurementLikelihood;
        simulatedCell = allSimulatedCells{l}.cell;
        previousPathIdx = allSimulatedCells{l}.Index;
        particle = allSimulatedCells{l}.particle;
        
        
        if (model.EstimateMeasurementNoise == 1)
            measurementVariance = particle.MeasurementVariance;
        end
        
        
        xPath = simulatedCell.xPath;
        
        %xSampled = xPath(:, end);
        xSampled = simulatedCell.xSampled;
        
        prediction = model.OutputWeightMatrix' * xSampled;
        measurement = targetCell.Measurement(:, options.TimeIndex);
        
        [L, erf] = EvaluateMeasurementLikelihood(model, ...
            measurement, prediction, measurementVariance);
        
        LNew = sum(sum(log(L)));
        
        
        
        MeasurementLikelihoodPrior = LNew;
        
        if (currParticleIdx==1)
            alpha = 1;
        elseif (simulatedCell.Reject == 1)
            alpha = 0;
        else
            LOld = MeasurementLikelihood;
            alpha = min(1, exp(LNew - LOld));
        end
        
        if (rand<=alpha)
            acceptedCount = acceptedCount + 1;
            MeasurementLikelihood = LNew;
            
            CurrentState = particle.CurrentState;
            CurrentState(:, 1, 1) = xSampled(:, end);
            
            Paths = particle.Paths;
            
            Paths{1}.deltamRNA = [particle.Paths{1}.deltamRNA simulatedCell.deltamRNA];
            Paths{1}.a = simulatedCell.a;
            Paths{1}.b = simulatedCell.b;
            Paths{1}.LastA = simulatedCell.LastA;
            Paths{1}.LastB = simulatedCell.LastB;
            Paths{1}.PriorCellsA = particle.Paths{1}.PriorCellsA;
            Paths{1}.PriorCellsB = particle.Paths{1}.PriorCellsB;
            Paths{1}.xPath = simulatedCell.xPath;
            Paths{1}.tPath = simulatedCell.tPath;
            
            Paths{1}.Erf = [particle.Paths{1}.Erf, erf];
            
            %ErrorFunction(k) = sum(sum(Paths{k}.Erf));
            
            particleIndex = particle.Index;
            Paths{1}.PreviousPathIndex = previousPathIdx;
            
            RandomEffectHyper = particle.RandomEffectHyper;
            MorphHyper = particle.MorphHyper;
            MeasurementVariance = particle.MeasurementVariance;
            %fprintf('+');
        else
            
            %fprintf('-');
        end
        
        %fprintf('(particle %d) \n', currParticleIdx);
        
        
        
        %% Set newly sampled states and parameters
        currentParticle.MeasurementLikelihood = MeasurementLikelihood;
        currentParticle.MeasurementLikelihoodPrior = MeasurementLikelihoodPrior;
        currentParticle.CurrentState = CurrentState;
        currentParticle.Paths = Paths;
        currentParticle.Index = particleIndex;
        currentParticle.RandomEffectHyper = RandomEffectHyper;
        currentParticle.MorphHyper = MorphHyper;
        currentParticle.MeasurementVariance = MeasurementVariance;
        %currentParticle = ApplyInvariantKernelToParticle(currentParticle, model, targetCells, options.SubMCMCRuns);
        
        
        %pPosterior = SetParticleSample(pPosterior, currentParticle, currParticleIdx);
        
        % HACK: Much faster than using the function!
        idx = currParticleIdx;
        pPosterior.IndexChain(idx) = currentParticle.Index;

        pPosterior.MeasurementLikelihoodChain(:, idx) = currentParticle.MeasurementLikelihood;
        pPosterior.MeasurementLikelihoodChainPrior(:, idx) = currentParticle.MeasurementLikelihoodPrior;
        pPosterior.CurrentStates(:, idx, 1) = currentParticle.CurrentState;
        pPosterior.Weights(idx) = currentParticle.Weight;
        pPosterior.Paths(:, idx) = currentParticle.Paths;

        pPosterior.RandomEffectHyperChain(:, idx) = currentParticle.RandomEffectHyper;
        pPosterior.MorphHyperChain(:, idx) = currentParticle.MorphHyper;

        pPosterior.MeasurementVarianceChain(:, idx) = currentParticle.MeasurementVariance;
        pPosterior.ErrorFunctionChain(:, idx) = currentParticle.ErrorFunction;
        
        

        %% Data Plotting
        
        if (mod(currParticleIdx, options.UpdateLength) == 0 && options.Plot == 2)
            axes(options.FigureHandle(1));
            p = plot(1:currParticleIdx, model.OutputWeightMatrix' * ...
                pPosterior.CurrentStates(:, 1:currParticleIdx, targetCellIdx));
            
            axes(options.FigureHandle(2));
            p = plot(1:currParticleIdx, pPosterior.MeasurementVarianceChain(1:currParticleIdx));
            
            axes(options.FigureHandle(3));
            p = plot(1:currParticleIdx, pPosterior.RandomEffectHyperChain(1, 1:currParticleIdx));
            
            axes(options.FigureHandle(4));
            p = plot(1:currParticleIdx, pPosterior.RandomEffectHyperChain(2, 1:currParticleIdx));
            
            axes(options.FigureHandle(5));
            p = plot(1:currParticleIdx, pPosterior.MorphHyperChain(1, 1:currParticleIdx));
            
            axes(options.FigureHandle(6));
            p = plot(1:currParticleIdx, pPosterior.MorphHyperChain(2, 1:currParticleIdx));
            
            
            drawnow;
            %fprintf('\n');
        end
    end
    
    fprintf('Finished block %d (%d/%d particles) out of %d.\n', u, ...
        currentBlockSize, acceptedCount, numBlocks);
end



if (options.BurnIn < 1)
   burnIn = ceil(M * options.BurnIn); 
else
   burnIn = options.BurnIn; 
end

pPosterior = FinalizeParticleDistribution(pPosterior, maxTime, burnIn);

pPosterior.ModelEvidence = [pDist.ModelEvidence, mean(exp(sum(pPosterior.MeasurementLikelihoodChainPrior, 1)))];

fprintf('\n');


end
