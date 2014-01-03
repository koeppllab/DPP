function [pDistOut] = RunDPP(pDist, model, cells, options, histPlotOptions, ...
    qPlotOptions)


    callbackProvided = ~isempty(options.CallBackTime);

    if (isfield(options, 'NumMeasurements'))
       numMeasurements = options.NumMeasurements; 
    else
        numMeasurements = length(cells{1}.MeasurementTime);
    end

    options.N = 30000;
    
    
    for k=1:floor(numMeasurements/v)
        options.TimeIndex = k;
        qPlotOptions.numPoints = k*10;
        
        
        pDist = UpdatePosterior(pDist, model, cells, options);
        
        if (callbackProvided)
            options.CallBackTime(pDist, k);
        end
        
        
        fprintf('\nFinished time point %d.\n',k);
        
        if (options.Plot > 0)
            fprintf('Plotting intermediate results...\n');
            PlotPropagatedPriorHistograms(pDist, model, histPlotOptions);
        end
        
        if (options.StorePaths == 1)
            if (options.Plot > 0)
                PlotStateStatistics(pDist, qPlotOptions, cells);
            end
        elseif (options.StorePaths == 2)
            save(['distTmp_' num2str(k) '.mat'], 'pDist');
            if (options.Plot > 0)
                PlotStateStatistics(pDist, qPlotOptions, cells);
            end
        end
        
    end
    
    
    % If paths were stored on the hard drive (instead of the memory), load data
    % and reconstruct the full path distribution.
    if (options.StorePaths == 2)
        pDist = ReconstructFullPaths(1:options.NumMeasurements, model);
        qPlotOptions.numPoints = options.NumMeasurements * ...
            options.StorePathNumPoints;
        PlotStateStatistics(pDist, qPlotOptions, cells);
    end
    
    pDistOut = pDist;

end