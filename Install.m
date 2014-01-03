%% %%%%%%%%%%%%%%%%%%%%%%%%%   Install.m   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Objective: Compilation of the required mex files.
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

