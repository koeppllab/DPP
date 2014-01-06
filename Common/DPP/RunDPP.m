%% %%%%%%%%%%%%%%%%%%%%%%%%%   RunDPP.m   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Objective: Main function of the DPP algorithm.
%--------------------------------------------------------------------------
% Description: Takes all the options, model specifications and data and
% performs the DPP algorithm.
% IMPORTANT NOTE: Up to now, RunDPP only implements the "simultaneous"-mode
% (see tutorial for a description). The "cells-first" and "timepoint-first"
% modes will be implemented in a future version.
%--------------------------------------------------------------------------
% DPP - A MATLAB toolbox for the inference of heterogeneous reaction 
% kinetics.
%
% Copyright (C) 2014  Christoph Zechner
% 
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
%--------------------------------------------------------------------------
% For feedback and questions please contact me at:
%   zechner@control.ee.ethz.ch
%--------------------------------------------------------------------------


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