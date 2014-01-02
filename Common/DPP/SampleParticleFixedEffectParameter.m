function Rates = SampleParticleFixedEffectParameter(particle, model)

    A = zeros(length(model.UnknownRateIndex), length(particle.Paths));
    B = A;

    
    ATotal = particle.A + model.APrior;
    BTotal = particle.B + model.BPrior;

    if (sum(ATotal == 0) == 0 && sum(BTotal == 0) == 0)
        Rates = SampleParameterConditional(ATotal, BTotal);
    else
        Rates = particle.Rates;
    end

end