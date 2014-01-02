function fullDist = ReconstructFullPaths(distIndex, model)


   
    numDists = length(distIndex);
    
    %% load temporary files
    
    dists = cell(numDists, 1);
    
    for i=1:numDists
       load(['distTmp_' num2str(distIndex(i)) '.mat']);
       dists{i} = pDist;
       clear pDist;
    end

    numCells = dists{end}.NumCells;
    numParticles = dists{end}.NumParticles;
    
    Paths = cell(numCells, numParticles);
    fullDist = dists{end};
    
    
    for l=1:numParticles
        

        for k=1:numCells
            
            currParticle = GetParticleSample(dists{end}, l);
            Paths{k, l} = currParticle.Paths{k};
            Paths{k, l}.Probability = 1;
            
            for i=2:numDists
                
                idx = numDists - i + 1;
                dist = dists{idx};

                
                currParticle = GetParticleSample(dist, currParticle.Paths{k}.PreviousPathIndex);
            
                
                if (isfield(currParticle.Paths{k}, 'xPath'))
                    if (~isempty(currParticle.Paths{k}.xPath))
                        Paths{k, l}.xPath = [currParticle.Paths{k}.xPath, Paths{k, l}.xPath(:, 2:end)];
                        Paths{k, l}.tPath = [currParticle.Paths{k}.tPath, Paths{k, l}.tPath(:, 2:end)];
                    end
                end
            end
            
        end
        
        fprintf('Reconstructed path %d of %d...\n', l, numParticles);
    end
    
    
    for l=1:numParticles
        for k=1:numCells
            fullDist.Paths{k, l} = Paths{k, l};
        end
    end
end