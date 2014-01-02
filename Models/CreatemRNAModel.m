%% Create mRNA model

      %mRNA      P
X0 = [1       0       0         ]';

    GeneOn = 0.1;
    GeneOff = 0.1;
    Transcription = 0.05;
    mRNADegradation = 0.001;
    
    
    c = [
        
         GeneOn
         GeneOff
    
         Transcription
         mRNADegradation

         ];

           %geneOff geneOn  mRNA       
    Pre = [ 1       0       0         
            0       1       0          
            
            0       1       0         
            0       0       1         
           ];
        
    Post = [0       1       0          
            1       0       0          
        
            0       1       1        
            0       0       0          
            ];


S = Post - Pre;
[NumReactions, NumSpecies] = size(S);

speciesNames{1} = 'GeneOff';
speciesNames{2} = 'GeneOn';
speciesNames{3} = 'mRNA';

reactionNames{1} = 'Gene-On';
reactionNames{2} = 'Gene-Off';
reactionNames{3} = 'Transc.';
reactionNames{4} = 'Deg. (mRNA)';

kineticModel.Name = 'mRNA Model';
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

save mRNAModel.mat kineticModel;
