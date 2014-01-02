function Rates = SampleParameterConditional(a, b)

    Rates = zeros(length(a), 1);

    for l=1:length(a)
        Rates(l) = ...
            gamrnd(a(l), 1 / b(l));
    end

end