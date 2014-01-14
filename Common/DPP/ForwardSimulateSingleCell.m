function simulatedCell = ForwardSimulateSingleCell(particle, model, numCells, c, aPrior,...
    bPrior, startTime, endTime, options, targetCellIdx, targetCell)

N = options.N;


if (nargin == 10)
    targetCell = {};
    morph = 0;
elseif (nargin == 11)
    if (model.MorphologicalFeatures == 1)
        morph = 1;
    else
        morph = 0;
    end
end

CurrentState = particle.CurrentState;

sampleGrid = linspace(startTime, endTime, options.StorePathNumPoints);

Pre = model.Pre;
Post = model.Post;


X0 = CurrentState(:, 1, 1);

currentA = particle.Paths{1}.a;
currentB = particle.Paths{1}.b;


aPriorTmp = aPrior;
bPriorTmp = bPrior;

kineticRateIdx = ~model.RandomEffect;
kineticRateTotIdx = model.RandomEffectIndex(kineticRateIdx);
randomEffectIdx = logical(model.RandomEffect);

aPriorTmp(kineticRateIdx) = aPriorTmp(kineticRateIdx) + particle.Paths{1}.SumStatA(kineticRateTotIdx);
bPriorTmp(kineticRateIdx) = bPriorTmp(kineticRateIdx) + particle.Paths{1}.SumStatB(kineticRateTotIdx);


if (morph==1)
    aPriorTmp(randomEffectIdx) = aPriorTmp(randomEffectIdx) + particle.MorphHyper(1);
    bPriorTmp(randomEffectIdx) = bPriorTmp(randomEffectIdx) + ...
        particle.MorphHyper(2) / targetCell.MorphologicalFeatureValue;
end


[xPath, tPath, rPath, G] = ...
    SimulateMarginalProcess(X0, c, Pre,...
    Post, N, endTime, model.InputParams.inputRateIndex, ...
    model.InputParams.inputTimes, model.InputParams.inputLevels, ...
    model.InputParams.inputTimes, ...
    startTime, model.RandomEffectIndex, ...
    aPriorTmp, bPriorTmp,...
    currentA(model.RandomEffectIndex), ...
    currentB(model.RandomEffectIndex));

if (tPath(end) ~= endTime)
    simulatedCell.Reject = 1;
else
    simulatedCell.Reject = 0;
end

if (options.StorePaths > 0)
    
    if (options.StorePathNumPoints > 0)
        xSampled = ...
            SampleCTMPPathGrid(xPath, tPath, sampleGrid);
        tSampled = sampleGrid;
    else
        xSampled = xPath;
        tSampled = tPath;
    end
    
    if (sum(options.StorePathCellIdx == targetCellIdx) == 1)
        if (options.StorePaths == 1)
            tSampled = ...
                [particle.Paths{1}.tPath(1:end-1), tSampled];
            xSampled = ...
                [particle.Paths{1}.xPath(:, 1:end-1), xSampled];
        end
    end
else
    xSampled = xPath(:, end);
    tSampled = tPath(:, end);
end


xSampledMeasurement = SampleCTMPPathGrid(xPath, tPath, model.MeasurementTime(options.TimeIndex));


simulatedCell.StateMeasurements = ...
    [particle.Paths{1}.StateMeasurements, xSampledMeasurement];

simulatedCell.xPath = xSampled;
simulatedCell.tPath = tSampled;
simulatedCell.EndPoint = xSampled;
simulatedCell.xSampled = xSampledMeasurement;


[a, b] = CalculatePathStatistics(tPath, rPath, G);
simulatedCell.LastA = a;
simulatedCell.LastB = b;
simulatedCell.a = a + currentA;
simulatedCell.b = b + currentB;

%fprintf('*');
end

