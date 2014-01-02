%% Create transcriptional input function (double pulse)


addpath('../../Common/');
addpath('../../Common/Statistics');
addpath('../../Common/StochChemKin');
addpath('../../Common/StochChemKin/SMC_MixedEffect');
addpath('../../Models');


modelData = load('SimpleModel.mat');
kineticModel = modelData.kineticModel;

inputData = load('InputSimpleModel.mat');
inputParams = inputData.inputParams;
inputParams.inputRateIndex = kineticModel.InputRateIdx;

T = 100*60;

measuredSpecies = 'Protein';
W = strcmp(measuredSpecies, kineticModel.SpeciesNames)';


MeasurementTime = linspace(0, T, 25);
MeasurementTime = MeasurementTime(2:end);

%% initialize alogorithm

targetReactionIdx = [2; 3; 4];

aPrior = [1; 2; 3];
bPrior = [50; 200; 100];

randomEffect = [0; 0; 0];

model = InitializeMEGSM(kineticModel, MeasurementTime, ...
    W, [0.15, 2, 0.03], 'logn', targetReactionIdx, randomEffect, aPrior, bPrior, ...
    [log(1);log(1000)], [2, 0; 0, 2], [], ...
    [], [], @GetPieceWiseConstantInput, inputParams);

numCells = 5;

cells = SimulateCells(model, numCells);

save MeasurementsSimpleModel.mat cells;
