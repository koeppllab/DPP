%% Create kinetic model

      %mRNA      P
X0 = [0         0     ]';

Transcription = 1; % will be modulated by external function
mRNADegradation = 0.005;
Translation = 0.01;
ProteinDegradation = 0.001;


c = [
    Transcription
    mRNADegradation
    Translation
    ProteinDegradation
    
    ];

      %mRNA        P
Pre = [
    0           0
    1           0
    
    1           0
    0           1
    ];

Post = [
    1           0
    0           0
    
    1           1
    0           0
    
    ];

S = Post - Pre;
[NumReactions, NumSpecies] = size(S);

speciesNames{1} = 'mRNA';
speciesNames{2} = 'Protein';

reactionNames{1} = 'Transc.';
reactionNames{2} = 'Deg. (mRNA)';
reactionNames{3} = 'Transl.';
reactionNames{4} = 'Deg. (Protein)';

kineticModel.Name = 'Simple Expression Model';
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

save SimpleModel.mat kineticModel;
