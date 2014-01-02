function [pDistOut] = RunDPP_SingleCell(pDist, model, cells, options, histPlotOptions, ...
    qPlotOptions)


    callbackProvided = ~isempty(options.CallBackTime);

    if (isfield(options, 'NumMeasurements'))
       numMeasurements = options.NumMeasurements; 
    else
        numMeasurements = length(cells{1}.MeasurementTime);
    end
    
    numCells = length(cells);
    options.N = 20000;
    
    for l=1:numCells
    
        pDist = InitializePriorParticleDistribution(pDist, model, l);
        pDist.CurrentTime = 0;
 
        v = 1;
        for k=1:floor(numMeasurements/v)
            options.TimeIndex = [(k-1)*v+1:k*v]; 
            qPlotOptions.numPoints = k*8*v;

            pDist = UpdatePosteriorSingleCell(pDist, model, cells, l, options);
           
            if (callbackProvided)
               options.CallBackTime(pDist, k); 
            end


            fprintf('\nFinished cell %d time point %d.\n',l, k);

            if (options.Plot > 0)
                fprintf('Plotting intermediate results...\n');
                PlotPropagatedPriorHistograms(pDist, model, histPlotOptions);
            end

            if (options.StorePaths == 1)
                if (options.Plot > 0)
                    PlotStateStatistics(pDist, qPlotOptions, cells); 
                end
            elseif (options.StorePaths == 2)
                save(['distTmp_' num2str(l) '_' num2str(k) '.mat'], 'pDist');
                if (options.Plot > 0)
                    PlotStateStatistics(pDist, qPlotOptions, cells);
                end
            end

        end
    
    end
    
    
    
    pDistOut = pDist;

end