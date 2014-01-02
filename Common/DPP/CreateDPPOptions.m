function options = CreateDPPOptions(handle)

    if (nargin==0)
       options.Plot = 0; 
    else
       options.Plot = 1;
       figure(handle);
    end

    options.M = 100;
    options.BlockSize = 20;
    options.UpdateLength = 100;
    options.BurnIn = round(options.M * 0.3);
    if (options.Plot == 1)
        options.FigureHandle(1) = subplot(2, 3, 1);
        options.FigureHandle(2) = subplot(2, 3, 2);
        options.FigureHandle(3) = subplot(2, 3, 3);
        options.FigureHandle(4) = subplot(2, 3, 4);
        options.FigureHandle(5) = subplot(2, 3, 5);
        options.FigureHandle(6) = subplot(2, 3, 6);
    end
    options.TimeIndex = 1;
    options.StorePaths = 1;
    options.StorePathNumPoints = 20;
    options.StorePathCellIdx = 1;
    options.SubMCMCRuns = 5;
    options.CallBackTime = [];
    options.CallBackIteration = [];
    options.ParticleReduce = 1;

end