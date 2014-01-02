function simulatedCells = ForwardSimulateState(particle, model, numCells, c, aPrior,...
    bPrior, startTime, endTime, targetCells)

    N = 50000;


    if (nargin == 8)
       targetCells = {};
       morph = 0;
    elseif (nargin == 9)
        if (model.MorphologicalFeatures == 1)
            morph = 1; 
        else
            morph = 0;
        end
    end

    CurrentState = particle.CurrentState;

    Pre = model.Pre;
    Post = model.Post;


    simulatedCells = cell(numCells, 1);

    
    for k=1:numCells

        X0 = CurrentState(:, 1, k);

        currentA = particle.Paths{k}.a;
        currentB = particle.Paths{k}.b;
        
        aPriorTmp = aPrior;
        bPriorTmp = bPrior;
        
        if (morph==1)
           aPriorTmp = aPriorTmp + particle.MorphHyper(1);
           bPriorTmp = bPriorTmp + ...
               particle.MorphHyper(2) / targetCells{k}.MorphologicalFeatureValue;
        end

        [xPath, tPath, rPath, ~, G] = ...
            SimulatePopulationTimeVaryingSyncNR(X0, c, Pre, Post, N, ...
            endTime, ...
            model.InputParams.inputRateIndex,...
            model.InputFun, model.InputParams, ...
            model.InputParams.inputTimes, startTime, ...
            model.RandomEffectIndex, aPrior, bPrior, currentA, currentB);

        simulatedCells{k}.xPath = xPath;
        simulatedCells{k}.tPath = tPath;
        simulatedCells{k}.rPath = rPath;
        simulatedCells{k}.G = G;
        simulatedCells{k}.EndPoint = xPath(:, end);
        
        [a, b] = CalculatePathStatistics(tPath, rPath, G);
        simulatedCells{k}.LastA = a;
        simulatedCells{k}.LastB = b;
        simulatedCells{k}.a = a + currentA;
        simulatedCells{k}.b = b + currentB;
        
        %fprintf('*');
    end

end