%% %%%%%%%%%%%%%%%%%%%%%%%%%   Install.m   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Objective: Compilation of the required mex files.
%--------------------------------------------------------------------------
% Author: Christoph Zechner
%--------------------------------------------------------------------------

clear all;

% Delete old mex files
evalc('!rm -f Common/StochasticSimulation/*.mex*');

% Compile *.c files
mex -output Common/StochasticSimulation/SimulateMarginalProcess ...
    MexFiles/SimulateMarginalProcess.c;
mex -output Common/StochasticSimulation/SampleCTMPPathGrid ...
    MexFiles/SampleCTMPPathGrid.c
mex -output Common/StochasticSimulation/CountMRNA ...
    MexFiles/CountMRNA.c;

