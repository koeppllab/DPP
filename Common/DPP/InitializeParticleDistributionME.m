function pDist = InitializeParticleDistributionME(numParticles, model, numCells)

    randomEffectCount = length(model.RandomEffectIndex);
    
    pDist.NumParticles = numParticles;
    
    pDist.RandomEffectIndex = model.RandomEffectIndex;
    pDist.NumCells = numCells;
    
    pDist.IndexChain = (1:numParticles)';
    
    
    pDist.ModelEvidence = [];

    if (sum(sum(model.HyperPriorSigma))>0)
       pDist.RandomEffectHyperChain = DrawLogNormal(model.HyperPriorMu, model.HyperPriorSigma, numParticles);
       %pDist.RandomEffectHyperChain(1, :) = lognrnd(model.HyperPrior(1,1), model.HyperPrior(1, 2), 1, numParticles);
       %pDist.RandomEffectHyperChain(2, :) = lognrnd(model.HyperPrior(2,1), model.HyperPrior(2, 2), 1, numParticles);
    else
        pDist.RandomEffectHyperChain = repmat([1; 100], 1, numParticles); %0.2*rand(2, numParticles);
    end
    
    if (sum(sum(model.MorphPriorSigma))>0)
       pDist.MorphHyperChain = DrawLogNormal(model.MorphPriorMu, model.MorphPriorSigma, numParticles);
       %pDist.RandomEffectHyperChain(1, :) = lognrnd(model.HyperPrior(1,1), model.HyperPrior(1, 2), 1, numParticles);
       %pDist.RandomEffectHyperChain(2, :) = lognrnd(model.HyperPrior(2,1), model.HyperPrior(2, 2), 1, numParticles);
    else
        pDist.MorphHyperChain = repmat([1; 100], 1, numParticles); %0.2*rand(2, numParticles);
    end
    
    
    pDist.PathLikelihoodChain = zeros(numCells, numParticles);
    pDist.MeasurementLikelihoodChain = zeros(numCells, numParticles);
    pDist.MeasurementLikelihoodChainPrior = zeros(numCells, numParticles);
    pDist.Weights = ones(1, numParticles) / numParticles;
    
    pDist.CurrentStates = model.X0(:, ones(1, numParticles), ones(1, numCells));
    
    pDist.CurrentTime = 0;
    
    pathStruct.xPath = zeros(size(model.X0, 1), 0);
    pathStruct.tPath = [];
    pathStruct.rPath = [];
    pathStruct.a = zeros(size(model.c));
    pathStruct.b = zeros(size(model.c));
    pathStruct.PreviousPathIndex = -1;
    pathStruct.Erf = zeros(model.OutputDim, 0);
    pathStruct.PriorCellsA = [];
    pathStruct.PriorCellsB = [];
    pathStruct.deltamRNA = [];
    pathStruct.StateMeasurements = pathStruct.xPath;
    pathStruct.SumStatA = pathStruct.a;
    pathStruct.SumStatB = pathStruct.b;

    pDist.Paths = repmat({pathStruct}, pDist.NumCells, pDist.NumParticles);

    

    if (isfield(model, 'NoiseA') && isfield(model, 'NoiseB'))
        pDist.MeasurementVarianceChain = sqrt(1 ./ gamrnd(model.NoiseA, 1/model.NoiseB, 1, numParticles));
    else
        pDist.MeasurementVarianceChain = model.MeasurementVariance(ones(1, numParticles));
    end
    
    pDist.ErrorFunctionChain = zeros(1, numParticles);
end