function [p] = PlotStateStatistics(pDist, options, cells, targetIdx)

numPoints = options.numPoints;
q = options.quantile;
stateIdx = options.stateIdx;

if (options.plotReferenceCells > 0)
    targetCells = cells;
end

cellIdx = options.cellIdx;

if (nargin == 4)
    if (sum(cellIdx == targetIdx)>0)
        cellIdx = targetIdx;
    else
        return;
    end
end



cols = min(3, length(cellIdx));
rows = ceil(length(cellIdx) / cols);


for k=1:length(cellIdx)
    j = cellIdx(k);
    
    if (sum(pDist.Paths{j, 1}.a)==0)
        return;
    end
    
    for i=1:pDist.NumParticles
        
        path = pDist.Paths{j, i};
        
        xPath = path.xPath;
        tPath = path.tPath;
        
        grid = linspace(min(tPath), max(tPath), numPoints);
        [dPath] = SampleCTMPPathGrid(xPath, tPath, grid);
        
        dTime = grid;
        
        for l=1:length(stateIdx)
            cellsQuantile{k}.States{l}.Data(i, :) = dPath(stateIdx(l), :);
        end
    end
    
    
    h = options.FigureHandlesCells(k);
    title(h, ['Cell ' num2str(j)]);
    
    deltaCol = 1/length(stateIdx);
    
    for l=1:length(stateIdx)
        if (options.plotConfidence == 1)
            cellsQuantile{k}.States{l}.Quantiles = ...
                quantile(cellsQuantile{k}.States{l}.Data, q);
            
            Q = cellsQuantile{k}.States{l}.Quantiles;
            
            for i=2:length(q)
                
                Q(i, :) = cellsQuantile{k}.States{l}.Quantiles(i, :) ...
                    - cellsQuantile{k}.States{l}.Quantiles(i-1, :);
                
            end
            
            p = area(h, dTime / 60, Q');
            hold(h, 'on');
            set(p(1), 'FaceColor', 'none');
            set(p(1), 'LineWidth', 2);
            set(p(1), 'LineStyle', 'none');
            
            for i=2:length(q)
                factor = ...
                    2/length(q) * (min(i-1, length(q)-i+1) -1);
                factorVec(i-1) = factor;
                
                
                val = (1-factor);
                colVec = [val*0.7, val*0.4, val*0.4];
                set(p(i), 'FaceColor', val*[0.7, l*deltaCol*0.8, l*deltaCol*0.8]);
                set(p(i), 'LineWidth', 2);
                set(p(i), 'LineStyle', 'none');
            end
        end
        
        
        cellsQuantile{k}.States{l}.Mean = ...
            mean(cellsQuantile{k}.States{l}.Data);
        
        p = plot(h, dTime / 60, cellsQuantile{k}.States{l}.Mean);
        set(p, 'LineWidth', 2);
        set(p, 'Color', [0.3, l*deltaCol*0.3, l*deltaCol*0.3]);
        
        xlabel(h, 'Time in min');
        ylabel(h, 'Molecule Count');
    end
    
    
    %% Compute True Path
    
    if (options.plotReferenceCells > 0)
        hold(h, 'on');
        
        if (options.plotReferenceCells == 2)
            xPathTarget = targetCells{j}.xPath;
            tPathTarget = targetCells{j}.tPath;
            
            timeIdx = and(tPathTarget >= min(dTime), tPathTarget <= max(dTime));
            
            for l=1:length(stateIdx)
                p = stairs(h, tPathTarget(timeIdx) / 60, xPathTarget(stateIdx(l), timeIdx), '--');
                set(p, 'LineWidth', 2);
                set(p, 'Color', [0.3, l*deltaCol*0.5, l*deltaCol*0.5]);
            end
            
            hold(h, 'on');
            
        end
        
        measIdx = targetCells{j}.MeasurementTime <= max(dTime);
        
        p = plot(h, targetCells{j}.MeasurementTime(measIdx) / 60, ...
            targetCells{j}.Measurement(:, measIdx), 'o');
        
        set(p, 'MarkerSize', 8);
        set(p, 'LineWidth', 1);
        set(p, 'MarkerFaceColor', [0.5, 0.5, 0.5]);
        set(p, 'MarkerEdgeColor', 'k');
        
        hold(h, 'off');
    end
    
    lims = xlim(h);
    xlim(h, [0, lims(2)]);

    box(h, 'off');
    drawnow;
    
end



end