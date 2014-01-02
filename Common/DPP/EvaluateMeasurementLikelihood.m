function [L, erf] = EvaluateMeasurementLikelihood(model, measurement, prediction,...
    scalingParameter)

    if (nargin < 4)
       scalingParameter = model.MeasurementVariance; 
    end
    
    if (strcmp(model.MeasurementDensity, 'normal'))
        deviations = measurement - prediction;
        
        L = prod(normpdf(deviations, 0, scalingParameter), 1);
        erf = power(deviations, 2);
    elseif (strcmp(model.MeasurementDensity, 'logn'))

     	%L = prod(lognpdf(measurement, ...
        %        log(prediction+1), scalingParameter), 1);
       
        L = 1;
        erf(length(measurement)) = 0;
        
        for l=1:length(measurement)
            [LTmp, e] = lognpdf_fast(measurement(l), log(prediction(l)+1), scalingParameter);
            L = L*LTmp;
            erf(l) = e;
        end
    end
       

end