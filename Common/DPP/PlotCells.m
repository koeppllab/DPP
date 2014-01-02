function [p] = PlotCells(cells)

    for k=1:length(cells)
        if (cells{k}.Simulated == 1)
            len = size(cells{k}.xPath, 2);
            z = repmat(k, 1, len);

            p = plot3(cells{k}.tPath / 60, z, cells{k}.xPath);
            set(p, 'LineWidth', 2); 
        
        end
        hold on;
        len = size(cells{k}.MeasurementTime, 2);
        z = repmat(k, 1, len);
        
        p = plot3(cells{k}.MeasurementTime / 60, z, cells{k}.Measurement, 'o');
        set(p, 'MarkerSize', 8);
        set(p, 'MarkerFaceColor', 'r');
        set(p, 'MarkerEdgeColor', 'r');
        
        if (cells{k}.Simulated == 0)
           set(p, 'LineStyle', '--'); 
           set(p, 'Color', 'k');
           set(p, 'LineWidth', 1);
        end
        
    end
    
    xlabel('Time in min');
    ylabel('Cell Index');
    zlabel('Molecule Count');

    set(gca, 'CameraPosition', [850, -10, 860]);
    set(gca, 'CameraViewAngle', 11.141988726609318);
    drawnow;
    hold off;
end