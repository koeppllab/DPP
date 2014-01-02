function PlotPropagatedPriorHistograms(pDist, model, options)



if (options.plotRefValues == 0)
    c = Inf(size(pDist.Paths{1}.a));
    measurementVariance = Inf;
    hyperA = Inf;
    hyperB = Inf;
else
    c = model.c;
    measurementVariance = model.MeasurementVariance;
    rdIdx = model.RandomEffect == 1;
    hyperA = model.RandomEffectA(rdIdx);
    hyperB = model.RandomEffectB(rdIdx);
end

numBins = options.binVec(1);
num2DBins = options.binVec(2);

logTC = options.logTransform(1);
logTHyper = options.logTransform(2);
logTMorph = options.logTransform(3);
logTMV = options.logTransform(4);
%% Plot Kinetic Parameter Results

kineticRateIndex = model.RandomEffectIndex(model.RandomEffect == 0);
kineticRateA = model.RandomEffectA(model.RandomEffect == 0);
kineticRateB = model.RandomEffectB(model.RandomEffect == 0);

if (sum(model.RandomEffect) == 0)
   extr = 0; 
else
   extr = 1;
end

if (model.MorphologicalFeatures == 1)
   morph = 1; 
else
   morph = 0; 
end

cols = options.Columns;
n = length(kineticRateIndex);
if (options.PlotAllVSAll == 1)
    plotCount = (n^2 - n) / 2 + n + 3;
else
    plotCount = n + extr + morph + 1;
end
rows = ceil(plotCount / cols);



[Rates, HyperParams, MorphParams] = DrawParameterPosteriorSamples(pDist, kineticRateA, kineticRateB, ...
    kineticRateIndex, options.SampleMultiplier);

    
quantiles = [0.05, 0.95];


if (isfield(options, 'InputRateMultiplier'))
    for i=1:size(options.InputRateMultiplier, 1)
        idx = options.InputRateMultiplier(i, 1);
        factor = options.InputRateMultiplier(i, 2);
        Rates(idx, :) = Rates(idx, :) * factor;
    end
end

figure(options.FigureHandle);

for l=1:n
    subplot(rows, cols,l);
    data = Rates(l, :);
        
    mmseEst = mean(data);
    
    a = kineticRateA(l);
    b = kineticRateB(l);
    
    
    if (logTC == 1)
        data = log10(data);
    end
    
    [h, b] = EvaluateMarkovChain(data, 1, numBins);
    a = area(b, h);
    
    set(a, 'FaceColor', [0.3, 0.6, 0.9]);
    
    ylims = ylim;
    hold on;
    trueVal = c(kineticRateIndex(l));
    if (logTC == 1)
        trueVal = log10(trueVal);
        mmseEst = log10(mmseEst);
    end
    
    plot([trueVal, trueVal], ylims, 'r', 'LineWidth', 2);
    
    plot([mmseEst, mmseEst], ylims, '--', 'LineWidth', 2, 'Color', [0.5, 0.5, 0.5]);
    hold off;
    box off;
    xlabel(['c_' num2str(kineticRateIndex(l))]);
    %set(gca, 'PlotBoxAspectRatio', [1, 1, 1]);


end


if (n>1 && options.PlotAllVSAll == 1)
    
    perms = nchoosek([1:n], 2);
    
    numPerms = size(perms, 1);
    
    for k=1:numPerms
        l = perms(k, 1);
        i = perms(k, 2);
        
        subplot(rows, cols, n+k);
        
        data = Rates([i, l], :)';
        
        mmseValue(1) = mean(data(:, 1));
        mmseValue(2) = mean(data(:, 2));
        
        trueValue(1) = c(kineticRateIndex(i));
        trueValue(2) = c(kineticRateIndex(l));
        
        if (logTC == 1)
            data = log10(data);
            
            mmseValue = log10(mmseValue);
            trueValue = log10(trueValue);
        end
        
        [counts, centers] = hist3(data, [num2DBins, num2DBins]);
        contourf(centers{2}, centers{1}, counts, 30, 'LineColor', 'none'); hold on;
        CM = colormap('gray');
        CM = 1-CM;
        colormap(CM);
        xlabel(['c_' num2str(kineticRateIndex(l))]);
        ylabel(['c ' num2str(kineticRateIndex(i))]);
        
        hold on;
        
        
        xlims = xlim;
        ylims = ylim;
        
        plot([trueValue(2), trueValue(2)], ylims, 'r', 'LineWidth', 2);
        plot([mmseValue(2), mmseValue(2)], ylims, '--', ...
            'Color', [0.5, 0.5, 0.5], 'LineWidth', 2);
        
        plot(xlims, [trueValue(1), trueValue(1)], 'r', 'LineWidth', 2);
        plot(xlims, [mmseValue(1), mmseValue(1)], '--', ...
            'Color', [0.5, 0.5, 0.5], 'LineWidth', 2);
        box off;
        %set(gca, 'PlotBoxAspectRatio', [1, 1, 1]);
        hold off;
    end
    
    plotOffset = (n^2 - n) / 2 + n;
else
    plotOffset = n;
end



if (sum(model.RandomEffect) > 0)
    
    %% Plot HyperParameters
    
    subplot(rows, cols, plotOffset + 1);
    
    data = HyperParams([2, 1], :);
    
    if (isfield(options, 'TransformHyperParameters'))
        if (options.TransformHyperParameters == 1)
            a = data(2, :);
            b = data(1, :);
            
            meanV = a./b;
            CV = sqrt(a./(b.^2)) ./ meanV;
            
            data(2, :) = meanV;
            data(1, :) = CV;
        end
    end
    
    data = data';
    
    
    mmseValue(1) = mean(data(:, 1));
    mmseValue(2) = mean(data(:, 2));
    
    trueValue(1) = hyperB;
    trueValue(2) = hyperA;
    
    
    
    if (logTHyper == 1)
        data = log10(data);
        
        mmseValue = log10(mmseValue);
        trueValue = log10(trueValue);
    end
    
    [counts, centers] = hist3(data, [num2DBins, num2DBins]);
    contourf(centers{2}, centers{1}, counts, 30, 'LineColor', 'none'); hold on;
    CM = colormap('gray');
    CM = 1-CM;
    colormap(CM);
    xlabel('\alpha');
    ylabel('\beta');
    
    hold on;
    
    
    xlims = xlim;
    ylims = ylim;
    
    plot([trueValue(2), trueValue(2)], ylims, 'r', 'LineWidth', 2);
    plot([mmseValue(2), mmseValue(2)], ylims, '--', ...
        'Color', [0.5, 0.5, 0.5], 'LineWidth', 2);
    
    plot(xlims, [trueValue(1), trueValue(1)], 'r', 'LineWidth', 2);
    plot(xlims, [mmseValue(1), mmseValue(1)], '--', ...
        'Color', [0.5, 0.5, 0.5], 'LineWidth', 2);
    box off;
    %set(gca, 'PlotBoxAspectRatio', [1, 1, 1]);
    hold off;
   
    
end


if (model.MorphologicalFeatures == 1)
    
    %% Plot HyperParameters
    
    subplot(rows, cols, plotOffset + 2);
    data = MorphParams([2, 1], :);%pDist.RandomEffectHyperChain([2, 1], :)';
    
    if (isfield(options, 'TransformHyperParameters'))
        if (options.TransformHyperParameters == 1)
            a = data(2, :);
            b = data(1, :);
            
            meanV = a./b;
            CV = sqrt(a./(b.^2)) ./ meanV;
            
            data(2, :) = meanV;
            data(1, :) = CV;
        end
    end
    
    data = data';
    
    mmseValue(1) = mean(data(:, 1));
    mmseValue(2) = mean(data(:, 2));
    
    trueValue(1) = model.MorphB;
    trueValue(2) = model.MorphA;
    
    

    if (logTMorph == 1)
        data = log10(data);
        
        mmseValue = log10(mmseValue);
        trueValue = log10(trueValue);
    end
    
    [counts, centers] = hist3(data, [num2DBins, num2DBins]);
    contourf(centers{2}, centers{1}, counts, 30, 'LineColor', 'none'); hold on;
    CM = colormap('gray');
    CM = 1-CM;
    colormap(CM);
    xlabel('u');
    ylabel('v');
    
    hold on;
    
    
    xlims = xlim;
    ylims = ylim;
    
    plot([trueValue(2), trueValue(2)], ylims, 'r', 'LineWidth', 2);
    plot([mmseValue(2), mmseValue(2)], ylims, '--', ...
        'Color', [0.5, 0.5, 0.5], 'LineWidth', 2);
    
    plot(xlims, [trueValue(1), trueValue(1)], 'r', 'LineWidth', 2);
    plot(xlims, [mmseValue(1), mmseValue(1)], '--', ...
        'Color', [0.5, 0.5, 0.5], 'LineWidth', 2);
    box off;
    %set(gca, 'PlotBoxAspectRatio', [1, 1, 1]);
    hold off;
   
    
end



%% Plot Measurement Variance
if (model.EstimateMeasurementNoise == 1)
    subplot(rows, cols, plotOffset + extr + morph + 1);
    data = pDist.MeasurementVarianceChain;

    mmseEst = mean(data);


    a = model.NoiseA;
    b = model.NoiseB;


    if (logTC == 1)
        data = log10(data);
    end

    [h, b] = EvaluateMarkovChain(data, 1, numBins);
    a = area(b, h);

    set(a, 'FaceColor', [0.3, 0.6, 0.9]);

    ylims = ylim;
    hold on;
    trueVal = measurementVariance;
    if (logTMV == 1)
        trueVal = log10(trueVal);
        mmseEst = log10(mmseEst);
    end

    plot([trueVal, trueVal], ylims, 'r', 'LineWidth', 2);

    plot([mmseEst, mmseEst], ylims, '--', 'LineWidth', 2, 'Color', [0.5, 0.5, 0.5]);
    hold off;
    box off;
    xlabel('\omega');

end

drawnow;


end