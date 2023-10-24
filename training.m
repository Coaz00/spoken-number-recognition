function LPC_train = training(p)
    LPC_train = zeros(p+1,3);
    folder = "train\";
    f = ["jedan","dva","osam"];
    for k = 1:3

        path = folder + f(k);
        d = dir(path);
        sequences = d([d.isdir]==0);

        LPC_new = zeros(p+1,1);
        for i = 1:length(sequences)
            seq = sequences(i).name;
            seq = folder + f(k) + "\" + seq;
            [x,fs] = audioread(seq);

            rec = preprocessing(x,fs);
            LPC = mean(feature_extraction(rec,fs),2);

            LPC_new = LPC_new + LPC;
        end
        LPC_train(:,k) = LPC_new/length(sequences); % centri klasa u p+1-dimenzionalnom prostoru
    end

end