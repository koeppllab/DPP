%% Create transcriptional input function (double pulse)


T = 100*60;

numInputSteps = 4;

J = 0.1;

inputParams.inputLevels = [1, 0, 1, 0];
inputParams.inputLevels = J * inputParams.inputLevels;
inputParams.inputTimes = T / numInputSteps * [1:numInputSteps];


save InputSimpleModel.mat inputParams;
