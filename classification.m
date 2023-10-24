function pred = classification(x,fs,LPC_train)
    LPC = mean(feature_extraction(x,fs),2);
    
    dist = zeros(1,3);
    for i = 1:length(dist)
        dist(i) = sum((LPC - LPC_train(:,i)).^2);
    end
    [~,pred] = min(dist);
    if pred == 3
        pred = 8;
    end
end