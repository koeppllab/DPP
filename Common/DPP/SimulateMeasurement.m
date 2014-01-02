function Measurement = SimulateMeasurement(model, xPath, tPath)

    xSampled = SampleCTMPPathGrid(xPath, tPath, model.MeasurementTime);


    if (strcmp(model.MeasurementDensity, 'logn'))
        mean = log(model.OutputWeightMatrix'*xSampled + 1);
        Measurement = random(model.MeasurementDensity, ...
            mean, model.MeasurementVariance, ...
            model.OutputDim, length(model.MeasurementTime));
    elseif (strcmp(model.MeasurementDensity, 'logn0'))
        mean = log(model.OutputWeightMatrix'*xSampled);
        Measurement = random('logn', ...
            mean, model.MeasurementVariance, ...
            model.OutputDim, length(model.MeasurementTime));
    else
        Measurement = model.OutputWeightMatrix'*xSampled + ...
            random(model.MeasurementDensity, 0, model.MeasurementVariance, ...
            model.OutputDim, length(model.MeasurementTime));
    end
    


end