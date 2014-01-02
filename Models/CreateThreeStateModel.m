%% Create kinetic model

	X0 = [1       0         0   0         0     ]';

    GeneOn = 0.006;
    GeneOnR = 0.001;
    GeneOff = 0.005;
    GeneOffR = 0.001;
    Transcription = 1;
    mRNADegradation = 0.02;
    Translation = 0.004;
    ProteinDegradation = 0.0004;
    
    
    c = [
        
         GeneOn
         GeneOnR
         GeneOff
         GeneOffR
         
         Transcription
         mRNADegradation
         Translation
         ProteinDegradation

         ];

           %geneOff geneOn  geneOnR     mRNA        P  
    Pre = [ 1       0       0            0           0   
            0       1       0            0           0   
            
            0       1       0            0           0
            0       0       1            0           0
            
            0       0       1            0           0   
            0       0       0            1           0   
            
            0       0       0            1           0   
            0       0       0            0           1     
           ];
        
    Post = [0       1       0            0           0   
            1       0       0            0           0   
       
            0       0       1            0           0
            0       1       0            0           0
            
            0       0       1            1           0   
            0       0       0            0           0   
            
            0       0       0            1           1   
            0       0       0            0           0   
            ];

S = Post - Pre;
[NumReactions, NumSpecies] = size(S);

speciesNames{1} = 'GeneOff';
speciesNames{2} = 'GeneOn';
speciesNames{3} = 'GeneOnR';
speciesNames{4} = 'mRNA';
speciesNames{5} = 'Protein';

reactionNames{1} = 'Gene-On';
reactionNames{2} = 'Gene-Off';
reactionNames{3} = 'Gene-On-R';
reactionNames{4} = 'Gene-Off-R';
reactionNames{5} = 'Transc.';
reactionNames{6} = 'Deg. (mRNA)';
reactionNames{7} = 'Transl.';
reactionNames{8} = 'Deg. (Protein)';

kineticModel.Name = 'Three-State Model';
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

save ThreeStateModel.mat kineticModel;
