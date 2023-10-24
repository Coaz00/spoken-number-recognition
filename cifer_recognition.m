function pred = cifer_recognition()
    fs = 8000;
    p = 15;
    speech_record(fs,16,1,1.5,'test.wav');

    [x,fs] = audioread('test.wav');
    rec = preprocessing(x,fs);
    LPC_train = training(p);
    pred = classification(rec,fs,LPC_train);

end