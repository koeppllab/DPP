%%%%%%%%%%%%%%%%%%%%% ForwardSimulateSingleCellFD.m %%%%%%%%%%%%%%%%%%%%%%%
% Objective: Extend a sample path until the next time point.
%--------------------------------------------------------------------------
% Description: Takes a particle and extends the sample path of the
% specified cell ("targetCell") until "endTime".
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

function simulatedCell = ForwardSimulateSingleCellFD(particle, model, ... 
    numCells, c, aPrior, bPrior, startTime, endTime, options, ...
    targetCellIdx, targetCell)

N = options.N;


if (nargin == 10)
    targetCell = {};
    morph = 0;
elseif (nargin == 11)
    if (model.MorphologicalFeatures == 1)
        morph = 1;
    else
        morph = 0;
    end
end

CurrentState = particle.CurrentState;

sampleGrid = linspace(startTime, endTime, options.StorePathNumPoints);

Pre = model.Pre;
Post = model.Post;


X0 = CurrentState(:, 1, targetCellIdx);

currentA = particle.Paths{targetCellIdx}.a;
currentB = particle.Paths{targetCellIdx}.b;


aPriorTmp = aPrior;
bPriorTmp = bPrior;

kineticRateIdx = ~model.RandomEffect;
kineticRateTotIdx = model.RandomEffectIndex(kineticRateIdx);
randomEffectIdx = logical(model.RandomEffect);

for k=1:numCells
    if (k==targetCellIdx)
        continue;
    end
    aPriorTmp(kineticRateIdx) = aPriorTmp(kineticRateIdx) + ...
        particle.Paths{k}.a(kineticRateTotIdx);
    bPriorTmp(kineticRateIdx) = bPriorTmp(kineticRateIdx) + ...
        particle.Paths{k}.b(kineticRateTotIdx);
end


if (morph==1)
    aPriorTmp(randomEffectIdx) = aPriorTmp(randomEffectIdx) + particle.MorphHyper(1);
    bPriorTmp(randomEffectIdx) = bPriorTmp(randomEffectIdx) + ...
        particle.MorphHyper(2) / targetCell.MorphologicalFeatureValue;
end


[xPath, tPath, rPath, G] = ...
    SimulateMarginalProcess(X0, c, Pre,...
    Post, N, endTime, model.InputParams.inputRateIndex, ...
    model.InputParams.inputTimes, model.InputParams.inputLevels, ...
    model.InputParams.inputTimes, ...
    startTime, model.RandomEffectIndex, ...
    aPriorTmp, bPriorTmp,...
    currentA(model.RandomEffectIndex), ...
    currentB(model.RandomEffectIndex));

% It might be the case that more reaction fire than the specified N-value.
% In this case, the sample is rejected.
if (tPath(end) ~= endTime)
    simulatedCell.Reject = 1;
else
    simulatedCell.Reject = 0;
end


if (options.StorePaths > 0)
    
    if (options.StorePathNumPoints > 0)
        [xSampled] = ...
            SampleCTMPPathGrid(xPath, tPath, sampleGrid);
        
        tSampled = sampleGrid;
    else
        xSampled = xPath;
        tSampled = tPath;
    end
    
    if (sum(options.StorePathCellIdx == targetCellIdx) == 1)
        if (options.StorePaths == 1)
            tSampled = ...
                [particle.Paths{targetCellIdx}.tPath(1:end-1), tSampled];
            xSampled = ...
                [particle.Paths{targetCellIdx}.xPath(:, 1:end-1), xSampled];
        end
    end
else
    xSampled = xPath(:, end);
    tSampled = tPath(:, end);
end


xSampledMeasurement = SampleCTMPPathGrid(xPath, tPath, model.MeasurementTime(options.TimeIndex));


simulatedCell.StateMeasurements = ...
    [particle.Paths{targetCellIdx}.StateMeasurements, xSampledMeasurement];

simulatedCell.xPath = xSampled;
simulatedCell.tPath = tSampled;

simulatedCell.EndPoint = xSampled;
simulatedCell.xSampled = xSampledMeasurement;


[a, b] = CalculatePathStatistics(tPath, rPath, G);
simulatedCell.LastA = a;
simulatedCell.LastB = b;
simulatedCell.a = a + currentA;
simulatedCell.b = b + currentB;

end

