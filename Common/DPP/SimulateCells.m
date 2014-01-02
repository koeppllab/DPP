function [cells] = SimulateCells(model, numCells, numRuns)

    

    cells = cell(numCells, 1);

    N = 200000;
    
    extrIdx = model.RandomEffectIndex(logical(model.RandomEffect));
    APrior = model.RandomEffectA(find(model.RandomEffect));
    BPrior = model.RandomEffectB(find(model.RandomEffect));
    
    if (nargin==2)

        parfor k=1:numCells

            c = model.c;

            %% Simulate parameters that are subject to cell-2-cell variability
            for l=1:length(model.RandomEffect)
                if (model.RandomEffect(l) == 1)
                    c(model.RandomEffectIndex(l)) = gamrnd(model.RandomEffectA(l), 1 / model.RandomEffectB(l));
                    if (model.MorphologicalFeatures == 1)
                        cells{k}.MorphologicalFeatureValue = 1 / gamrnd(model.MorphA, ...
                            1 /(model.MorphB * c(model.RandomEffectIndex(l))));
                    end
                end
            end

%             [xPath, tPath, rPath, h, g, cVarying] = ...
%                 SimulateGillespieTimeVaryingSync(model.X0, c, model.Pre, ...
%                 model.Post, 300000, model.T, 1, ...
%                 model.InputFun, model.InputParams, model.InputParams.inputTimes);

            
            [xPath, tPath, rPath] = ... 
            SimulateMarginalProcess(model.X0, c, model.Pre,...
                   model.Post, N, model.T, model.InputParams.inputRateIndex, ...
                   model.InputParams.inputTimes, model.InputParams.inputLevels, ...
                   model.InputParams.inputTimes, ...
                   0, [], ...
                   [], [],...
                   0, ...
                   0);

            cells{k}.xPath = xPath;
            cells{k}.tPath = tPath;
            cells{k}.rPath = rPath;

            cells{k}.Simulated = 1;
            cells{k}.c = c;


            cells{k}.Measurement = SimulateMeasurement(model, xPath, tPath);
            cells{k}.MeasurementTime = model.MeasurementTime;

            fprintf('Simulated cell %d.\n', k);
        end


    else
        geneIdx = 3;
        mRNAIdx = 4;
        
        
        parfor k=1:numCells
            
            % required to remove PARFOR warnings
            startTime = -1;
            
            c = model.c;

            %% Simulate parameters that are subject to cell-2-cell variability
            for l=1:length(model.RandomEffect)
                if (model.RandomEffect(l) == 1)
                    c(model.RandomEffectIndex(l)) = gamrnd(model.RandomEffectA(l), 1 / model.RandomEffectB(l));
                
                    if (model.MorphologicalFeatures == 1)
                        cells{k}.MorphologicalFeatureValue = 1 / gamrnd(model.MorphA, ...
                            1 /(model.MorphB * c(model.RandomEffectIndex(l))));
                    end
                end
            end



            for l=1:numRuns
                %[xPath, tPath] = ...
                %    SimulateGillespieTimeVaryingSync(model.X0, c, model.Pre, ...
                %    model.Post, 500000, model.T, 1, ...
                %    model.InputFun, model.InputParams, model.InputParams.inputTimes);



                [xPath, tPath] = ...
                    SimulateMarginalProcess(model.X0, c, model.Pre,...
                    model.Post, N, model.T, model.InputParams.inputRateIndex, ...
                    model.InputParams.inputTimes, model.InputParams.inputLevels, ...
                    model.InputParams.inputTimes, ...
                    0, [], ...
                    [], [],...
                    0, ...
                    0);

                %cells{k}.Runs{l}.xPath = xPath;
                %cells{k}.Runs{l}.tPath = tPath;
                %cells{k}.Runs{l}.rPath = rPath;
                [xSampled, tSampled] = SampleCTMPPathGrid(xPath, tPath, model.MeasurementTime);

                cells{k}.Runs{l}.xPath = xSampled;
                cells{k}.Runs{l}.tPath = tSampled;
                cells{k}.Runs{l}.Measurement = SimulateMeasurement(model, xPath, tPath);
                cells{k}.Runs{l}.MeasurementTime = model.MeasurementTime;

                fprintf('Simulated cell %d - Run: %d.\n', k, l);
            end

            cells{k}.Simulated = 1;
            cells{k}.c = c;

        end
    end

end