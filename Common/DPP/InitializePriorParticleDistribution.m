function pDist = InitializePriorParticleDistribution(pDistOld, model, cellIdx)


    pDist = pDistOld;
    
    numSpecies = model.NumSpecies;
    numParticles = pDist.NumParticles;
    
    if (cellIdx == 1)
       return; 
    end
    
    
    for i=1:size(pDist.Paths, 2)
       
        pDist.Paths{1, i}.xPath = zeros(numSpecies, 0);
        pDist.Paths{1, i}.tPath = [];
        pDist.Paths{1, i}.rPath = [];
        pDist.Paths{1, i}.PriorCellsA = [pDist.Paths{1, i}.PriorCellsA, pDist.Paths{1, i}.a];
        pDist.Paths{1, i}.PriorCellsB = [pDist.Paths{1, i}.PriorCellsB, pDist.Paths{1, i}.b];
        pDist.Paths{1, i}.SumStatA = pDist.Paths{1, i}.SumStatA + pDist.Paths{1, i}.a;
        pDist.Paths{1, i}.SumStatB = pDist.Paths{1, i}.SumStatB + pDist.Paths{1, i}.b;
        pDist.Paths{1, i}.a = zeros(size(model.c));
        pDist.Paths{1, i}.b = zeros(size(model.c));
    end
    
    pDist.CurrentStates = model.X0(:, ones(1, numParticles), 1);
    pDist.CurrentTime = 0;
    
end