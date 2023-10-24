function preds = testing(set)
    p = 15;
    LPC_train = training(p);
    
    folder = set+"\";
    f = ["jedan","dva","osam"];
    preds = zeros(1,15);
    for k = 1:3
        path = folder + f(k);
        d = dir(path);
        sequences = d([d.isdir]==0);
        for i = 1:length(sequences)
            seq = sequences(i).name;
            seq = folder + f(k) + "\" + seq;
            [x,fs] = audioread(seq);
            rec = preprocessing(x,fs);
            preds(i+(k-1)*length(sequences)) = classification(rec,fs,LPC_train);
        end

    end
end