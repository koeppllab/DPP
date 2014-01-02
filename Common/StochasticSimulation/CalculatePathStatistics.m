function [R, GInt] = CalculatePathStatistics(tPath, rPath, g)

    [NumReactions, NumEvents] = size(g);

    R = zeros(NumReactions, 1);    
    GInt = zeros(NumReactions, 1);
    
    if (NumEvents == 0)
       return; 
    end
    
    GInt = sum(g.*repmat(diff(tPath), NumReactions, 1), 2);


    for k=1:NumReactions
        R(k) = sum(rPath == k);
    end

    
end