clear all;
close all;

numSpecies = 1;
numEvents = 30;

x = unidrnd(3, 1, numEvents);
tau = exprnd(10, 1, numEvents-1);
tPath = [0, cumsum(tau)];

stairs(tPath, x');

grid = linspace(0, max(tPath), 100);

xSampled = SampleCTMPPathGrid(x, tPath, grid);


hold on;
plot(grid, xSampled, 'or');