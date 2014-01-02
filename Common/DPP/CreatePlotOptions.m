function qPlotOptions = CreatePlotOptions(handles)

    if (nargin==0)
        handles(1) = figure;
    end

    qPlotOptions.quantile = [0.05, 0.95];
    qPlotOptions.plotReferenceCells = 1;
    qPlotOptions.stateIdx = [1, 2];
    qPlotOptions.cellIdx = [1];
    qPlotOptions.FigureHandlesCells = handles;
    qPlotOptions.plotPopulation = 0;
    qPlotOptions.plotConfidence = 1;
    qPlotOptions.numPoints = 10;

end