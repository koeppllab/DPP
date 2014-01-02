function model = InitializeMEGSM(kineticModel, ...
    MeasurementTime, ...
    OutputWeightMatrix, ...
    MeasurementVariance, ...
    MeasurementDensity, ...
    randomEffectParams, ...
    randomEffect, ...
    randomEffectA,...
    randomEffectB,...
    HyperPriorMu, HyperPriorSigma, ...
    MorphHyper, MorphPriorMu, MorphPriorSigma, ...
    inputFun, inputParams)

    model.KineticModel = kineticModel;

    model.X0 = kineticModel.X0;
    model.Pre = kineticModel.Pre;
    model.Post = kineticModel.Post;
    model.S = kineticModel.S;
    model.NumReactions = kineticModel.NumReactions;
    model.NumSpecies = kineticModel.NumSpecies;
    model.c = kineticModel.c;
    model.OutputWeightMatrix = OutputWeightMatrix;
    model.OutputDim = size(OutputWeightMatrix, 2);
    
    if (length(MeasurementVariance) == 1)
        model.MeasurementVariance = MeasurementVariance;
        model.EstimateMeasurementNoise = 0;
    else
        model.NoiseA = MeasurementVariance(2);
        model.NoiseB = MeasurementVariance(3);
        model.MeasurementVariance = MeasurementVariance(1);
        model.EstimateMeasurementNoise = 1;
    end
    
    model.MeasurementDensity = MeasurementDensity;
    model.MeasurementTime = MeasurementTime;
    model.T = max(MeasurementTime);
    
    model.RandomEffectIndex = randomEffectParams;
    model.RandomEffectA = randomEffectA;
    model.RandomEffectB = randomEffectB;
    model.RandomEffect = randomEffect;
    model.InputFun = inputFun;
    model.InputParams = inputParams;
    
    model.HyperPriorMu = HyperPriorMu;
    model.HyperPriorSigma = HyperPriorSigma;
    
    if (~isempty(MorphHyper))
        model.MorphologicalFeatures = 1;
        model.MorphA = MorphHyper(1);
        model.MorphB = MorphHyper(2);
        model.MorphPriorMu = MorphPriorMu;
        model.MorphPriorSigma = MorphPriorSigma;
    else
        model.MorphologicalFeatures = 0;
        model.MorphA = [];
        model.MorphB = [];
        model.MorphPriorMu = [];
        model.MorphPriorSigma = [];
    end
   
end