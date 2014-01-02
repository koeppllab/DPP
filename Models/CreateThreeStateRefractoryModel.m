%% Create kinetic model

	X0 = [1       0         0   0         0     ]';

    GeneOn = 0.006;
    GeneRef = 0.005;
    GeneOff = 0.001;
    Transcription = 1;
    mRNADegradation = 0.02;
    Translation = 0.004;
    ProteinDegradation = 0.0004;
    
    
    c = [
        
         GeneOn
         GeneRef
         GeneOff
         
         Transcription
         mRNADegradation
         Translation
         ProteinDegradation

         ];

           %geneOff geneOn  geneRef     mRNA        P  
    Pre = [ 1       0       0            0           0   
            0       1       0            0           0   
            
            0       0       1            0           0
            
            0       1       0            0           0   
            0       0       0            1           0   
            
            0       0       0            1           0   
            0       0       0            0           1     
           ];
        
    Post = [0       1       0            0           0   
            0       0       1            0           0   
       
            1       0       0            0           0
            
            0       1       0            1           0   
            0       0       0            0           0   
            
            0       0       0            1           1   
            0       0       0            0           0   
            ];

S = Post - Pre;
[NumReactions, NumSpecies] = size(S);

speciesNames{1} = 'GeneOff';
speciesNames{2} = 'GeneOn';
speciesNames{3} = 'GeneRef';
speciesNames{4} = 'mRNA';
speciesNames{5} = 'Protein';

reactionNames{1} = 'Gene-On';
reactionNames{2} = 'Gene-On-Ref';
reactionNames{3} = 'Gene-Off';
reactionNames{4} = 'Transc.';
reactionNames{5} = 'Deg. (mRNA)';
reactionNames{6} = 'Transl.';
reactionNames{7} = 'Deg. (Protein)';

kineticModel.Name = 'Three-State-Refractory Model';
kineticModel.Pre = Pre;
kineticModel.Post = Post;
kineticModel.S = S;
kineticModel.c = c;
kineticModel.X0 = X0;
kineticModel.NumSpecies = NumSpecies;
kineticModel.NumReactions = NumReactions;
kineticModel.SpeciesNames = speciesNames;
kineticModel.ReactionNames = reactionNames;
kineticModel.InputRateIdx = 1; %Gene-On Rate

save ThreeStateRefractoryModel.mat kineticModel;
