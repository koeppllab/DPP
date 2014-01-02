function [Stat, modelMMSE] = GetParameterStatistics(pDist, model)


%% Plot Kinetic Parameter Results

kineticRateIndex = model.RandomEffectIndex(model.RandomEffect == 0);
kineticRateA = model.RandomEffectA(model.RandomEffect == 0);
kineticRateB = model.RandomEffectB(model.RandomEffect == 0);

heterogeneousRateIndex = model.RandomEffectIndex(model.RandomEffect == 1);

n = length(kineticRateIndex);


modelMMSE = model;

[Rates, HyperParams, MorphParams] = DrawParameterPosteriorSamples(pDist, kineticRateA, kineticRateB, ...
    kineticRateIndex, 3);

    
quantiles = [0.05, 0.95];

modelMMSE.LogEvidence = sum(log(pDist.ModelEvidence));

for l=1:n
    data = Rates(l, :);
       

    
    a = kineticRateA(l);
    b = kineticRateB(l);
    
    meanPrior = a/b;
    mmseEst = mean(data);
    
    dataPrior = gamrnd(a, 1/b, 10000, 1);
    
    qPrior = quantile(dataPrior, quantiles);
    qPosterior = quantile(data, quantiles);
    
    
    Stat.KineticParameters{l}.Name = ...
        model.KineticModel.ReactionNames{kineticRateIndex(l)};
    Stat.KineticParameters{l}.QuantilesPrior = qPrior;
    Stat.KineticParameters{l}.MeanPrior = meanPrior;
    Stat.KineticParameters{l}.QuantilesPosterior = qPosterior;
    Stat.KineticParameters{l}.MeanPosterior = mmseEst;
    
    modelMMSE.c(kineticRateIndex(l)) = mmseEst;
end



if (sum(model.RandomEffect) > 0)
    
    %% Plot HyperParameters
    
    data = HyperParams([1, 2], :);
    data = data';
    
    
    mmseValue(1) = mean(data(:, 1));
    mmseValue(2) = mean(data(:, 2));
    
    
    dataPrior = DrawLogNormal(model.HyperPriorMu, model.HyperPriorSigma, 10000);
    
    
    qPrior = quantile(dataPrior', quantiles);
    qPosterior = quantile(data, quantiles);
    meanPrior = mean(dataPrior');
    
    
    Stat.ExtrinsicStatistics.A.QuantilesPrior = qPrior(:, 1)';
    Stat.ExtrinsicStatistics.A.MeanPrior = meanPrior(1);
    Stat.ExtrinsicStatistics.A.QuantilesPosterior = qPosterior(:, 1)';
    Stat.ExtrinsicStatistics.A.MeanPosterior = mmseValue(1);
    
    
    Stat.ExtrinsicStatistics.B.QuantilesPrior = qPrior(:, 2)';
    Stat.ExtrinsicStatistics.B.MeanPrior = meanPrior(2);
    Stat.ExtrinsicStatistics.B.QuantilesPosterior = qPosterior(:, 2)';
    Stat.ExtrinsicStatistics.B.MeanPosterior = mmseValue(2);
    
    modelMMSE.RandomEffectA(heterogeneousRateIndex) = mmseValue(1);
    modelMMSE.RandomEffectB(heterogeneousRateIndex) = mmseValue(2);
    
end


if (model.MorphologicalFeatures == 1)
    
    %% Plot HyperParameters
    
    data = MorphParams([1, 2], :);
    data = data';
    
    mmseValue(1) = mean(data(:, 1));
    mmseValue(2) = mean(data(:, 2));
    

    dataPrior = DrawLogNormal(model.MorphPriorMu, model.MorphPriorSigma, 10000);
    
    
    qPrior = quantile(dataPrior', quantiles);
    qPosterior = quantile(data, quantiles);
    meanPrior = mean(dataPrior');
    
    
    Stat.MorphologicalParameters.A.QuantilesPrior = qPrior(:, 1)';
    Stat.MorphologicalParameters.A.MeanPrior = meanPrior(1);
    Stat.MorphologicalParameters.A.QuantilesPosterior = qPosterior(:, 1)';
    Stat.MorphologicalParameters.A.MeanPosterior = mmseValue(1);
    
    
    Stat.MorphologicalParameters.B.QuantilesPrior = qPrior(:, 2)';
    Stat.MorphologicalParameters.B.MeanPrior = meanPrior(2);
    Stat.MorphologicalParameters.B.QuantilesPosterior = qPosterior(:, 2)';
    Stat.MorphologicalParameters.B.MeanPosterior = mmseValue(2);
    
    modelMMSE.MorphA = Stat.MorphologicalParameters.A.MeanPosterior;
    modelMMSE.MorphB = Stat.MorphologicalParameters.B.MeanPosterior;

end



%% Plot Measurement Variance

data = pDist.MeasurementVarianceChain;


mmseEst = mean(data);
modelMMSE.MeasurementVariance = mmseEst;

a = model.NoiseA;
b = model.NoiseB;
      
dataPrior = 1 ./ sqrt(gamrnd(a, 1/b, 10000, 1));

qPrior = quantile(dataPrior, quantiles);
qPosterior = quantile(data, quantiles);
meanPrior = mean(dataPrior);

Stat.MeasurementError.QuantilesPrior = qPrior;
Stat.MeasurementError.MeanPrior = meanPrior;
Stat.MeasurementError.QuantilesPosterior = qPosterior;
Stat.MeasurementError.MeanPosterior = mmseEst;



