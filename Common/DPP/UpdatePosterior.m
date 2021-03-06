%% %%%%%%%%%%%%%%%%%%%%%%%%% UpdatePosterior.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Objective: Takes the current posterior distribution and updates it using
% the next measurements.
%--------------------------------------------------------------------------
% Description: The current belief is passed to the function as a
% pDist-structure and updated using the measurements at the next time
% point (the time index is specified in the "options"-structure). This
% function implements the "simultaneous"-mode of DPP, which means that all
% cells are processed at once. In particular, each cell is updated using an
% individual Metropolis-Hastings sampling step. 
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

function [pPosterior] = UpdatePosterior(pDist, model, targetCells, options)

callbackProvided = ~isempty(options.CallBackIteration);


M = options.M(1);


if (options.TimeIndex > 1)
   M = ceil(M * options.ParticleReduce); 
end


maxTime = max(targetCells{1}.MeasurementTime(max(options.TimeIndex)));

numCells = length(targetCells);

pPosterior = InitializeParticleDistributionME(M, model, numCells);
currentParticle = GetParticleSample(pDist, 1);


Paths = cell(1);
c = model.c;


aPrior = model.RandomEffectA;
bPrior = model.RandomEffectB;

measurementVariance = model.MeasurementVariance;

randomEffectIdx = model.RandomEffect == 1;

blockSize = options.BlockSize;
numBlocks = ceil(M / blockSize);
rest = mod(M, blockSize);
currentBlockSize = blockSize;

oldPSize = pDist.NumParticles;
particleIdx = randperm(oldPSize);


timeOffset = pDist.CurrentTime;

for u=1:numBlocks
    
    if (u==numBlocks && rest > 0)
        currentBlockSize = rest;
    end
    
    if (callbackProvided == 1)
        options.CallBackIteration();
    end
    
    [particles, particleIdx] = DrawFromParticleDistribution(pDist, currentBlockSize);
    
    %% Run MCMC scheme
    for l=1:currentBlockSize

        currParticleIdx = (u - 1) * blockSize + l;
        
        optionsTmp = options;
        modelTmp = model;
        
        particle = particles{l};
        prevIdx = particleIdx(l); %particleIdx(l);
        
        %% Resample static parameters.
        particle = ApplyInvariantKernelToParticleFD(particle, modelTmp, targetCells, optionsTmp.SubMCMCRuns);
        

        %% Resample Paths
        
        aPriorTmp = aPrior;
        bPriorTmp = bPrior;
        
        
        aPriorTmp(randomEffectIdx) = particle.RandomEffectHyper(1);
        bPriorTmp(randomEffectIdx) = particle.RandomEffectHyper(2);
        
        rndIdx = randperm(numCells);
        for i=1:numCells
            targetCellIdx = rndIdx(i);
            
            simulatedCell = ForwardSimulateSingleCellFD(particle, model, 1, c, ...
                aPriorTmp, bPriorTmp, timeOffset, maxTime, options, targetCellIdx, targetCells{targetCellIdx});


            if (model.EstimateMeasurementNoise == 1)
                measurementVariance = particle.MeasurementVariance;
            end

            %xSampled = xPath(:, end);
            xSampled = simulatedCell.xSampled;

            prediction = model.OutputWeightMatrix' * xSampled;
            measurement = targetCells{targetCellIdx}.Measurement(:, options.TimeIndex);

            [L, erf] = EvaluateMeasurementLikelihood(model, ...
                measurement, prediction, measurementVariance);

            LNew = sum(sum(log(L)));



            MeasurementLikelihoodPrior(targetCellIdx) = LNew;

            if (currParticleIdx==1)
                alpha = 1;
            elseif (simulatedCell.Reject == 1)
                alpha = 0;
            else
                LOld = MeasurementLikelihood(targetCellIdx);
                alpha = min(1, exp(LNew - LOld));
            end

            if (rand<=alpha)
                MeasurementLikelihood(targetCellIdx) = LNew;

            
                CurrentState(:, 1, targetCellIdx) = xSampled(:, end);

                Paths{targetCellIdx}.a = simulatedCell.a;
                Paths{targetCellIdx}.b = simulatedCell.b;
                Paths{targetCellIdx}.LastA = simulatedCell.LastA;
                Paths{targetCellIdx}.LastB = simulatedCell.LastB;
                Paths{targetCellIdx}.PriorCellsA = particle.Paths{targetCellIdx}.PriorCellsA;
                Paths{targetCellIdx}.PriorCellsB = particle.Paths{targetCellIdx}.PriorCellsB;
                Paths{targetCellIdx}.xPath = simulatedCell.xPath;
                Paths{targetCellIdx}.tPath = simulatedCell.tPath;

                Paths{targetCellIdx}.Erf = [particle.Paths{targetCellIdx}.Erf, erf];
                Paths{targetCellIdx}.StateMeasurements = simulatedCell.StateMeasurements;

                particleIndex = particle.Index;
                Paths{targetCellIdx}.PreviousPathIndex = prevIdx;

                RandomEffectHyper = particle.RandomEffectHyper;
                MorphHyper = particle.MorphHyper;
                MeasurementVariance = particle.MeasurementVariance;
                
                %fprintf('+');
            else
                %fprintf('-');
            end
        
        particle.Paths{targetCellIdx} = Paths{targetCellIdx};
        
        end
        %fprintf('\n');
        
        %% Set newly sampled states and parameters
        currentParticle.MeasurementLikelihood = MeasurementLikelihood;
        currentParticle.MeasurementLikelihoodPrior = MeasurementLikelihoodPrior;
        currentParticle.CurrentState = CurrentState;
        currentParticle.Paths = Paths;
        currentParticle.Index = particleIndex;
        currentParticle.RandomEffectHyper = RandomEffectHyper;
        currentParticle.MorphHyper = MorphHyper;
        currentParticle.MeasurementVariance = MeasurementVariance;
        
        %% Store particle in the new posterior distribution
        idx = currParticleIdx;
        pPosterior.IndexChain(idx) = currentParticle.Index;

        pPosterior.MeasurementLikelihoodChain(:, idx) = currentParticle.MeasurementLikelihood;
        pPosterior.MeasurementLikelihoodChainPrior(:, idx) = currentParticle.MeasurementLikelihoodPrior;
        pPosterior.CurrentStates(:, idx, :) = currentParticle.CurrentState;
        pPosterior.Weights(idx) = currentParticle.Weight;
        pPosterior.Paths(:, idx) = currentParticle.Paths;

        pPosterior.RandomEffectHyperChain(:, idx) = currentParticle.RandomEffectHyper;
        pPosterior.MorphHyperChain(:, idx) = currentParticle.MorphHyper;

        pPosterior.MeasurementVarianceChain(:, idx) = currentParticle.MeasurementVariance;
        pPosterior.ErrorFunctionChain(:, idx) = currentParticle.ErrorFunction;
        
        
        %% Data Plotting
        
        if (mod(currParticleIdx, options.UpdateLength) == 0 && options.Plot == 2)
            axes(options.FigureHandle(1));
            p = plot(1:currParticleIdx, model.OutputWeightMatrix' * ...
                pPosterior.CurrentStates(:, 1:currParticleIdx, targetCellIdx));
            
            axes(options.FigureHandle(2));
            p = plot(1:currParticleIdx, pPosterior.MeasurementVarianceChain(1:currParticleIdx));
            
            axes(options.FigureHandle(3));
            p = plot(1:currParticleIdx, pPosterior.RandomEffectHyperChain(1, 1:currParticleIdx));
            
            axes(options.FigureHandle(4));
            p = plot(1:currParticleIdx, pPosterior.RandomEffectHyperChain(2, 1:currParticleIdx));
            
            axes(options.FigureHandle(5));
            p = plot(1:currParticleIdx, pPosterior.MorphHyperChain(1, 1:currParticleIdx));
            
            axes(options.FigureHandle(6));
            p = plot(1:currParticleIdx, pPosterior.MorphHyperChain(2, 1:currParticleIdx));
            
            
            drawnow;
            %fprintf('\n');
        end
    end
    
    fprintf('Finished block %d (%d particles) out of %d.\n', u, ...
        currentBlockSize, numBlocks);
end

if (options.BurnIn < 1)
   burnIn = ceil(M * options.BurnIn); 
else
   burnIn = options.BurnIn; 
end

pPosterior = FinalizeParticleDistribution(pPosterior, maxTime, burnIn);
pPosterior.ModelEvidence = [pDist.ModelEvidence, mean(exp(sum(pPosterior.MeasurementLikelihoodChainPrior, 1)))];

fprintf('\n');


end
