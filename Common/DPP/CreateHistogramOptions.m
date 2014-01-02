function histPlotOptions = CreateHistogramOptions(handle)

    if (nargin==0)
       handle = figure; 
    end

    %indicates which parameters shall be shown logarithmically:
    %   1: kinetic rates
    %   2: extrinsic statistics
    %   3: morphological parameters
    %   4: measurement noise parameters
    histPlotOptions.logTransform = [0, 1, 1, 0];
    
    histPlotOptions.binVec = [25, 50];
    histPlotOptions.plotRefValues = 0;
    histPlotOptions.Columns = 3;
    histPlotOptions.SampleMultiplier = 3; %Subsampling
    histPlotOptions.PlotAllVSAll = 0;
    histPlotOptions.FigureHandle = handle;
end