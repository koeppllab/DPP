%% Create kinetic model

	X0 = [1       0       0         0     ]';

    GeneOn = 0.006;
    GeneOff = 0.005;
    Transcription = 1;
    mRNADegradation = 0.02;
    Translation = 0.004;
    ProteinDegradation = 0.0004;
    
    
    c = [
        
         GeneOn
         GeneOff
    
         Transcription
         mRNADegradation
         Translation
         ProteinDegradation

         ];

           %geneOff geneOn  mRNA        P  
    Pre = [ 1       0       0           0   
            0       1       0           0   
            
            0       1       0           0   
            0       0       1           0   
            
            0       0       1           0   
            0       0       0           1     
           ];
        
    Post = [0       1       0           0   
            1       0       0           0   
       
            0       1       1           0   
            0       0       0           0   
            
            0       0       1           1   
            0       0       0           0   
            ];

S = Post - Pre;
[NumReactions, NumSpecies] = size(S);

speciesNames{1} = 'GeneOff';
speciesNames{2} = 'GeneOn';
speciesNames{3} = 'mRNA';
speciesNames{4} = 'Protein';

reactionNames{1} = 'Gene-On';
reactionNames{2} = 'Gene-Off';
reactionNames{3} = 'Transc.';
reactionNames{4} = 'Deg. (mRNA)';
reactionNames{5} = 'Transl.';
reactionNames{6} = 'Deg. (Protein)';

kineticModel.Name = 'Two-State Model';
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

save TwoStateModel.mat kineticModel;
