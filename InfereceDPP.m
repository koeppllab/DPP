%% %%%%%%%%%%%%%%%%%%%%%%%%%   InerenceDPP.m   %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Objective: Example script for state- and parameter inference using DPP.
%--------------------------------------------------------------------------
% Description: The two-state gene expression model is used to infer the
% heterogeneous gene expression kinetics using synthetic data of multiple
% cells.
%--------------------------------------------------------------------------
% Author: Christoph Zechner, BISON Group, ETH Zurich
%--------------------------------------------------------------------------
% For feedback and questions please contact me at:
%   zechner@control.ee.ethz.ch
%--------------------------------------------------------------------------

clear all;
close all;

% Add required library paths
addpath('Common/');
addpath('Common/Statistics');
addpath('Common/StochasticSimulation');
addpath('Common/DPP');
addpath('Models');

%% Model and input specification

% Load two-state gene expression model
load TwoStateModel.mat;


% Specify input signal (e.g., TF abundance). The code below simulates a
% double-pulse.
T = 100*60;
numInputSteps = 4;
J = 0.005;

inputParams.inputLevels = [1, 0, 1, 0, 0];
inputParams.inputLevels = J * inputParams.inputLevels;
inputParams.inputTimes = T / numInputSteps * [1:numInputSteps];
inputParams.inputRateIndex = kineticModel.InputRateIdx;

%% Initialize state-space model
% Output-Weight matrix. The measurements are assumed to be a noisy readout
% of linear combinations of the state, i.e., Y ~ p(. | W'*X).
W = [0 0 0 1]'; %Set fourth species (i.e., protein) to be read out.
 
% For instance, if mRNA is measured as well, the Output-Weight matrix would
% be W = [0, 0, 1, 0;
%         0, 0, 0, 1];

% Specify the measurement time points. If experimental data is used, those
% must be extracted from the dataset.
MeasurementTime = linspace(0, T, 25);
MeasurementTime = MeasurementTime(2:end);

% Specify reactions rates which are to be estimated from the data. All other
% reaction rates will be set to the value set in the kineticModel
% structure.
targetReactionIndex = [4; 5];


% Specify prior parameters. Priors are assumed to be Gamma distributed, 
% i.e., Gamma(a, b). Currently, only Gamma-type priors are supported. 
% Further priors will be included in a future release.
% Hint: For the gamma distribution it holds that mean = a/b, var = a/b^2 
% and CV = 1/sqrt(a). For highly informative priors 'a' needs to be large 
% (and the other way round.
aPrior = [4; 3];
bPrior = [ 200; 300];

% Specify extrinsic reaction. Currently, only one extrinsic reactions is
% supported.
randomEffect = [ 0; 1]; 

% Specify hyperparameters for the extrinsic statistics. The hyperparameters
% correspond to a bivariate log-normal distribution LN(MuHyper, SigmaHyper).
MuHyper = [log(1);log(1000)];
SigmaHyper = [2, 0; 
              0, 2];
          
% Specify hyperparameters for the morphological parameters. The
% hyperparameters correspond to a bivariate log-normal distribution 
% LN(MuMorph, SigmaMorph). MorphParams will be ignored if the morphological
% parameters are estimated from the data. If MuMorph and SigmaMorph are
% empty, MorphParams will be used and assumed to be known in the inference.
% In this example, MorphParams will also be used to simulate the
% morpholical feature.
MorphParams = [5, 10000]; 
MuMorph = [log(10); log(5000)];
SigmaMorph = [3, 0; 0, 3];

% Specify the scaling parameter and the corresponding hyperparameters of
% the measurement noise model. Note that the scaling parameter will be
% ignored during inference but is used for generating the synthetic
% trajectories.
measurementParams = [0.15, 2, 0.03];
mDens = 'logn';

model = InitializeMEGSM(kineticModel, MeasurementTime, ...
    W, measurementParams, mDens, targetReactionIndex, randomEffect, ...
    aPrior, bPrior, MuHyper, SigmaHyper, MorphParams, ...
    MuMorph, SigmaMorph, @GetPieceWiseConstantInput, inputParams);


%% Simulate reference data. 
% This will simulate the time-lapse single cell and the corresponding
% morphological features (using MorphParams).
numCells = 3;
cells = SimulateCells(model, numCells);

% Plot cells
figure(1);
PlotCells(cells);


%% Set options for inference and plotting

% specify which cell to monitor
cellIdx = [1,2];

% specify which states to monitor
stateIdx = [3, 4]; %mRNA+Protein


% Create axes handle for plotting.
figure(2);
stateHandles = [subplot(1,2,1); subplot(1,2,2)];
qPlotOptions = CreatePlotOptions(stateHandles);
qPlotOptions.cellIdx = cellIdx;
qPlotOptions.stateIdx = stateIdx;


options = CreateDPPOptions(figure);
options.M = 1000;
options.BurnIn = 0.03;
options.ParticleReduce = 0.5;
options.StorePathCellIdx = cellIdx;
options.BlockSize = 500;
options.Plot = 1;


options.UpdateLength = 100;

% specifies the number of time points used for inference. If this field is 
% left empty, all time points will be used.
options.NumMeasurements = 5; 

% This will swap intermediate distributions to the hard drive. As a 
% consequence, the state distribution is only plotted for the last time 
% iteration during simulation. Subsequently, the intermediate files will be 
% loaded to reconstruct the full state distribution. 
%   0: don't store paths (parameter inference only)
%   1: completely store paths (may need a huge amount of memory)
%   2: swap paths (writes intermediate distributions to disk)
options.StorePaths = 2; 

histPlotOptions = CreateHistogramOptions;
histPlotOptions.plotRefValues = 1; %plot the underlying ture value;

% Initialize particle distribution at time zero.
pDist = InitializeParticleDistributionME(options.M, model, numCells);


%% Perform state- and parameter inference.
% Run the dynamic prior propagation algorithm. The final posterior
% distribution at time T is returned in pDist.
pDist = RunDPP(pDist, model, cells, options, histPlotOptions, qPlotOptions);


% Get inferred parameter statistics and save the results
ParameterStats = GetParameterStatistics(pDist, model);
save results/Experiment_TwoState.mat pDist ParameterStats model cells;