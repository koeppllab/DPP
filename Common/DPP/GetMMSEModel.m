function modelMMSE = GetMMSEModel(pDist, model)

    sigma = mean(pDist.MeasurementVarianceChain);
    a = mean(pDist.RandomEffectHyperChain(1, :));
    b = mean(pDist.RandomEffectHyperChain(2, :));
    u = mean(pDist.MorphHyperChain(1, :));
    v = mean(pDist.MorphHyperChain(2, :));
    
    kineticRateIndex = model.RandomEffectIndex(model.RandomEffect == 0);
    kineticRateA = model.RandomEffectA(model.RandomEffect == 0);
    kineticRateB = model.RandomEffectB(model.RandomEffect == 0);
    populationRateIndex = model.RandomEffectIndex(model.RandomEffect == 1);
    
    [Rates] = DrawParameterPosteriorSamples(pDist, kineticRateA, kineticRateB, ...
    kineticRateIndex, 5*pDist.NumParticles);
    
    modelMMSE = model;
    modelMMSE.c(kineticRateIndex) = mean(Rates, 2);
    modelMMSE.InputParams.inputLevels = ...
        modelMMSE.c(modelMMSE.InputParams.inputRateIndex) * ...
        modelMMSE.InputParams.inputLevels;
    modelMMSE.RandomEffectA(populationRateIndex) = a;
    modelMMSE.RandomEffectB(populationRateIndex) = b;
    modelMMSE.MorphA = u;
    modelMMSE.MorphB = v;
    modelMMSE.MeasurementVariance = sigma;
end