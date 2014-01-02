function PrintLatexTable(ParameterStats, model, plotRef)

if (plotRef == 1)
   startIdx = 1; 
   refStr = '& Reference';
   colStr = 'c';
else
   startIdx = 0; 
   refStr = '';
   colStr = '';
end


kineticRateIndex = model.RandomEffectIndex(model.RandomEffect == 0);
kineticRateA = model.RandomEffectA(model.RandomEffect == 0);
kineticRateB = model.RandomEffectB(model.RandomEffect == 0);




fprintf('\\begin{tabular}{@{\\vrule height 10.5pt depth4pt  width0pt}lc%s|ccc|ccc|}\n', colStr);
if (plotRef == 1)
    fprintf('&  & \\multicolumn{1}{c}{} &  \\multicolumn3c{Prior} &\\multicolumn3c{Posterior}\\\\ \n');
else
    fprintf('&  \\multicolumn{1}{c}{} &  \\multicolumn3c{Prior} &\\multicolumn3c{Posterior}\\\\ \n'); 
end
fprintf('\\cline{%d-%d} \\cline{%d-%d}', 3+startIdx, 3+startIdx+2, 6+startIdx, 6+startIdx+2);
fprintf('\\vrule depth 6pt width 0pt Name %s & Unit &$q_{5}$&$q_{95}$ &Mean &$q_{5}$&$q_{95}$& Mean\\\\ \n', refStr);
fprintf('\\hline \n');

for l=1:length(ParameterStats.KineticParameters);
    
    param = ParameterStats.KineticParameters{l};
    
    if (plotRef == 1)
        refVal = model.c(kineticRateIndex(l));
        fprintf('$\\parameters_{%d}$ &$%3.2e$ &$s^{-1}$ & $%13.2e$ & $%13.2e$ & $%13.2e$ &$%13.2e$ & $%13.2e$ & $%13.2e$ \\\\ \n', ...
            kineticRateIndex(l), refVal, param.QuantilesPrior(1), param.QuantilesPrior(2), param.MeanPrior,...
            param.QuantilesPosterior(1), param.QuantilesPosterior(2), param.MeanPosterior);
        
    else
        fprintf('$\\parameters_{%d}$ &$s^{-1}$ & $%13.2e$ & $%13.2e$ & $%13.2e$ &$%13.2e$ & $%13.2e$ & $%13.2e$ \\\\ \n', ...
            kineticRateIndex(l), param.QuantilesPrior(1), param.QuantilesPrior(2), param.MeanPrior,...
            param.QuantilesPosterior(1), param.QuantilesPosterior(2), param.MeanPosterior);
    end

end


if (sum(model.RandomEffectIndex) > 0)
    
    param = ParameterStats.ExtrinsicStatistics.A;
    refVal = model.RandomEffectA(find(model.RandomEffect));
    
    if (plotRef == 1)
        fprintf('$\\alpha$ &$%3.2e$ &$-$ & $%13.2e$ & $%13.2e$ & $%13.2e$ &$%13.2e$ & $%13.2e$ & $%13.2e$ \\\\ \n', ...
        refVal, param.QuantilesPrior(1), param.QuantilesPrior(2), param.MeanPrior,...
        param.QuantilesPosterior(1), param.QuantilesPosterior(2), param.MeanPosterior);
    else
        fprintf('$\\alpha$ &$-$ & $%13.2e$ & $%13.2e$ & $%13.2e$ &$%13.2e$ & $%13.2e$ & $%13.2e$ \\\\ \n', ...
        param.QuantilesPrior(1), param.QuantilesPrior(2), param.MeanPrior,...
        param.QuantilesPosterior(1), param.QuantilesPosterior(2), param.MeanPosterior);
    end
    
    param = ParameterStats.ExtrinsicStatistics.B;
    refVal = model.RandomEffectB(find(model.RandomEffect));
    
    if (plotRef == 1)
        fprintf('$\\beta$ &$%3.2e$ &$-$ & $%13.2e$ & $%13.2e$ & $%13.2e$ &$%13.2e$ & $%13.2e$ & $%13.2e$ \\\\ \n', ...
        refVal, param.QuantilesPrior(1), param.QuantilesPrior(2), param.MeanPrior,...
        param.QuantilesPosterior(1), param.QuantilesPosterior(2), param.MeanPosterior);
    else
        fprintf('$\\beta$ &$-$ & $%13.2e$ & $%13.2e$ & $%13.2e$ &$%13.2e$ & $%13.2e$ & $%13.2e$ \\\\ \n', ...
        param.QuantilesPrior(1), param.QuantilesPrior(2), param.MeanPrior,...
        param.QuantilesPosterior(1), param.QuantilesPosterior(2), param.MeanPosterior);
    end
    
end

if (model.MorphologicalFeatures > 0)
    
   param = ParameterStats.MorphologicalParameters.A;
    refVal = model.MorphA;
    
    if (plotRef == 1)
        fprintf('$\\rho$ &$%3.2e$ &$-$ & $%13.2e$ & $%13.2e$ & $%13.2e$ &$%13.2e$ & $%13.2e$ & $%13.2e$ \\\\ \n', ...
        refVal, param.QuantilesPrior(1), param.QuantilesPrior(2), param.MeanPrior,...
        param.QuantilesPosterior(1), param.QuantilesPosterior(2), param.MeanPosterior);
    else
        fprintf('$\\rho$ &$-$ & $%13.2e$ & $%13.2e$ & $%13.2e$ &$%13.2e$ & $%13.2e$ & $%13.2e$ \\\\ \n', ...
        param.QuantilesPrior(1), param.QuantilesPrior(2), param.MeanPrior,...
        param.QuantilesPosterior(1), param.QuantilesPosterior(2), param.MeanPosterior);
    end
    
    param = ParameterStats.MorphologicalParameters.B;
    refVal = model.MorphB;
    
    if (plotRef == 1)
        fprintf('$\\phi$ &$%3.2e$ &$-$ & $%13.2e$ & $%13.2e$ & $%13.2e$ &$%13.2e$ & $%13.2e$ & $%13.2e$ \\\\ \n', ...
        refVal, param.QuantilesPrior(1), param.QuantilesPrior(2), param.MeanPrior,...
        param.QuantilesPosterior(1), param.QuantilesPosterior(2), param.MeanPosterior);
    else
        fprintf('$\\phi$ &$-$ & $%13.2e$ & $%13.2e$ & $%13.2e$ &$%13.2e$ & $%13.2e$ & $%13.2e$ \\\\ \n', ...
        param.QuantilesPrior(1), param.QuantilesPrior(2), param.MeanPrior,...
        param.QuantilesPosterior(1), param.QuantilesPosterior(2), param.MeanPosterior);
    end
    
end

if (isfield(ParameterStats, 'MeasurementError'))
    
    param = ParameterStats.MeasurementError;
    refVal = model.MeasurementVariance;
    
    if (plotRef == 1)
        fprintf('$\\omega$ &$%3.2e$ &$-$ & $%13.2e$ & $%13.2e$ & $%13.2e$ &$%13.2e$ & $%13.2e$ & $%13.2e$ \\\\ \n', ...
        refVal, param.QuantilesPrior(1), param.QuantilesPrior(2), param.MeanPrior,...
        param.QuantilesPosterior(1), param.QuantilesPosterior(2), param.MeanPosterior);
    else
        fprintf('$\\omega$ &$-$ & $%13.2e$ & $%13.2e$ & $%13.2e$ &$%13.2e$ & $%13.2e$ & $%13.2e$ \\\\ \n', ...
        param.QuantilesPrior(1), param.QuantilesPrior(2), param.MeanPrior,...
        param.QuantilesPosterior(1), param.QuantilesPosterior(2), param.MeanPosterior);
    end
end


fprintf('\\hline \n');
fprintf('\\end{tabular} \n');

% 
% if (n>1)
%     
%     perms = nchoosek([1:n], 2);
%     
%     numPerms = size(perms, 1);
%     
%     for k=1:numPerms
%         l = perms(k, 1);
%         i = perms(k, 2);
%         
%         subplot(rows, cols, n+k);
%         
%         data = Rates([i, l], :)';
%         
%         mmseValue(1) = mean(data(:, 1));
%         mmseValue(2) = mean(data(:, 2));
%         
%         trueValue(1) = c(kineticRateIndex(i));
%         trueValue(2) = c(kineticRateIndex(l));
%         
%         if (logTC == 1)
%             data = log10(data);
%             
%             mmseValue = log10(mmseValue);
%             trueValue = log10(trueValue);
%         end
%         
%         [counts, centers] = hist3(data, [num2DBins, num2DBins]);
%         contourf(centers{2}, centers{1}, counts, 200, 'LineColor', 'none'); hold on;
%         CM = colormap('gray');
%         CM = 1-CM;
%         colormap(CM);
%         xlabel(['c_' num2str(kineticRateIndex(l))]);
%         ylabel(['c ' num2str(kineticRateIndex(i))]);
%         
%         hold on;
%         
%         
%         xlims = xlim;
%         ylims = ylim;
%         
%         plot([trueValue(2), trueValue(2)], ylims, 'r', 'LineWidth', 2);
%         plot([mmseValue(2), mmseValue(2)], ylims, '--', ...
%             'Color', [0.5, 0.5, 0.5], 'LineWidth', 2);
%         
%         plot(xlims, [trueValue(1), trueValue(1)], 'r', 'LineWidth', 2);
%         plot(xlims, [mmseValue(1), mmseValue(1)], '--', ...
%             'Color', [0.5, 0.5, 0.5], 'LineWidth', 2);
%         box off;
%         %set(gca, 'PlotBoxAspectRatio', [1, 1, 1]);
%         hold off;
%     end
% end
% 
% plotOffset = (n^2 - n) / 2 + n;
% 
% if (sum(model.RandomEffect) > 0)
%     
%     %% Plot HyperParameters
%     
%     subplot(rows, cols, plotOffset + 1);
%     
%     data = HyperParams([2, 1], :);
%     
%     if (isfield(options, 'TransformHyperParameters'))
%         if (options.TransformHyperParameters == 1)
%             a = data(2, :);
%             b = data(1, :);
%             
%             meanV = a./b;
%             CV = sqrt(a./(b.^2)) ./ meanV;
%             
%             data(2, :) = meanV;
%             data(1, :) = CV;
%         end
%     end
%     
%     data = data';
%     
%     
%     mmseValue(1) = mean(data(:, 1));
%     mmseValue(2) = mean(data(:, 2));
%     
%     trueValue(1) = hyperB;
%     trueValue(2) = hyperA;
%     
%     
%     dataPrior = DrawLogNormal(model.HyperPriorMu, model.HyperPriorSigma, 10000);
%     
%     
%     qPrior = quantile(dataPrior', quantiles);
%     qPosterior = quantile(data, quantiles);
%     meanPrior = mean(dataPrior');
%     
%     fprintf('$\\alpha$ &$-$ & $%13.2e$ & $%13.2e$ & $%13.2e$ &$%13.2e$ & $%13.2e$ & $%13.2e$ \\\\ \n', ...
%         qPrior(1, 1), qPrior(2, 1), meanPrior(1),...
%         qPosterior(1, 2), qPosterior(2, 2), mmseValue(2));
%     fprintf('$\\beta$ &$-$ & $%13.2e$ & $%13.2e$ & $%13.2e$ &$%13.2e$ & $%13.2e$ & $%13.2e$ \\\\ \n', ...
%         qPrior(1, 2), qPrior(2, 2), meanPrior(2),...
%         qPosterior(1, 1), qPosterior(2, 1), mmseValue(1));
%     
%     if (logTHyper == 1)
%         data = log10(data);
%         
%         mmseValue = log10(mmseValue);
%         trueValue = log10(trueValue);
%     end
%     
%     [counts, centers] = hist3(data, [num2DBins, num2DBins]);
%     contourf(centers{2}, centers{1}, counts, 200, 'LineColor', 'none'); hold on;
%     CM = colormap('gray');
%     CM = 1-CM;
%     colormap(CM);
%     xlabel('\alpha');
%     ylabel('\beta');
%     
%     hold on;
%     
%     
%     xlims = xlim;
%     ylims = ylim;
%     
%     plot([trueValue(2), trueValue(2)], ylims, 'r', 'LineWidth', 2);
%     plot([mmseValue(2), mmseValue(2)], ylims, '--', ...
%         'Color', [0.5, 0.5, 0.5], 'LineWidth', 2);
%     
%     plot(xlims, [trueValue(1), trueValue(1)], 'r', 'LineWidth', 2);
%     plot(xlims, [mmseValue(1), mmseValue(1)], '--', ...
%         'Color', [0.5, 0.5, 0.5], 'LineWidth', 2);
%     box off;
%     %set(gca, 'PlotBoxAspectRatio', [1, 1, 1]);
%     hold off;
%    
%     
% end
% 
% 
% if (model.MorphologicalFeatures == 1 && sum(sum(model.MorphPriorSigma))>0)
%     
%     %% Plot HyperParameters
%     
%     subplot(rows, cols, plotOffset + 2);
%     data = MorphParams([2, 1], :);%pDist.RandomEffectHyperChain([2, 1], :)';
%     
%     if (isfield(options, 'TransformHyperParameters'))
%         if (options.TransformHyperParameters == 1)
%             a = data(2, :);
%             b = data(1, :);
%             
%             meanV = a./b;
%             CV = sqrt(a./(b.^2)) ./ meanV;
%             
%             data(2, :) = meanV;
%             data(1, :) = CV;
%         end
%     end
%     
%     data = data';
%     
%     mmseValue(1) = mean(data(:, 1));
%     mmseValue(2) = mean(data(:, 2));
%     
%     trueValue(1) = model.MorphB;
%     trueValue(2) = model.MorphA;
%     
%     
%     dataPrior = DrawLogNormal(model.MorphPriorMu, model.MorphPriorSigma, 10000);
%     
%     
%     qPrior = quantile(dataPrior', quantiles);
%     qPosterior = quantile(data, quantiles);
%     meanPrior = mean(dataPrior');
%     
%     fprintf('$\\rho$ &$-$ & $%13.2e$ & $%13.2e$ & $%13.2e$ &$%13.2e$ & $%13.2e$ & $%13.2e$ \\\\ \n', ...
%         qPrior(1, 2), qPrior(2, 2), meanPrior(2),...
%         qPosterior(1, 2), qPosterior(2, 2), mmseValue(2));
%     fprintf('$\\phi$ &$-$ & $%13.2e$ & $%13.2e$ & $%13.2e$ &$%13.2e$ & $%13.2e$ & $%13.2e$ \\\\ \n', ...
%         qPrior(1, 1), qPrior(2, 1), meanPrior(1),...
%         qPosterior(1, 1), qPosterior(2, 1), mmseValue(1));
%     
%     if (logTMorph == 1)
%         data = log10(data);
%         
%         mmseValue = log10(mmseValue);
%         trueValue = log10(trueValue);
%     end
%     
%     [counts, centers] = hist3(data, [num2DBins, num2DBins]);
%     contourf(centers{2}, centers{1}, counts, 200, 'LineColor', 'none'); hold on;
%     CM = colormap('gray');
%     CM = 1-CM;
%     colormap(CM);
%     xlabel('u');
%     ylabel('v');
%     
%     hold on;
%     
%     
%     xlims = xlim;
%     ylims = ylim;
%     
%     plot([trueValue(2), trueValue(2)], ylims, 'r', 'LineWidth', 2);
%     plot([mmseValue(2), mmseValue(2)], ylims, '--', ...
%         'Color', [0.5, 0.5, 0.5], 'LineWidth', 2);
%     
%     plot(xlims, [trueValue(1), trueValue(1)], 'r', 'LineWidth', 2);
%     plot(xlims, [mmseValue(1), mmseValue(1)], '--', ...
%         'Color', [0.5, 0.5, 0.5], 'LineWidth', 2);
%     box off;
%     %set(gca, 'PlotBoxAspectRatio', [1, 1, 1]);
%     hold off;
%    
%     
% end
% 
% 
% 
% %% Plot Measurement Variance
% 
% subplot(rows, cols, plotOffset + 3);
% data = pDist.MeasurementVarianceChain;
% 
% mmseEst = mean(data);
% 
% 
% a = model.NoiseA;
% b = model.NoiseB;
%     
%     
% dataPrior = 1 ./ sqrt(gamrnd(a, 1/b, 10000, 1));
% 
% qPrior = quantile(dataPrior, quantiles);
% qPosterior = quantile(data, quantiles);
% meanPrior = mean(dataPrior);
% 
% fprintf('$\\omega$ &$-$ & $%13.2e$ & $%13.2e$ & $%13.2e$ &$%13.2e$ & $%13.2e$ & $%13.2e$ \\\\ \n', ...
%     qPrior(1), qPrior(2), meanPrior,...
%     qPosterior(1), qPosterior(2), mmseEst);
% 
% 
% 
% if (logTC == 1)
%     data = log10(data);
% end
% 
% [h, b] = EvaluateMarkovChain(data, 1, numBins);
% a = area(b, h);
% 
% set(a, 'FaceColor', [0.3, 0.6, 0.9]);
% 
% ylims = ylim;
% hold on;
% trueVal = measurementVariance;
% if (logTMV == 1)
%     trueVal = log10(trueVal);
%     mmseEst = log10(mmseEst);
% end
% 
% plot([trueVal, trueVal], ylims, 'r', 'LineWidth', 2);
% 
% plot([mmseEst, mmseEst], ylims, '--', 'LineWidth', 2, 'Color', [0.5, 0.5, 0.5]);
% hold off;
% box off;
% xlabel('\omega');
% 
% 
% drawnow;
% 
% end