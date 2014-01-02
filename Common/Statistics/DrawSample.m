function X = DrawSample(M, probabilities)
    
    edges = [0; cumsum(probabilities)];
    U = rand(1, M);
    
    [~, X] = histc(U, edges);

end