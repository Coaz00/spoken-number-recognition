function speech_record(fs,nbits,nchans,duration,filename)

% Snimanje govornog signala
x = audiorecorder(fs,nbits,nchans); 
disp('start speaking');
recordblocking(x,duration);
disp('end of recording');

% Cuvanje snimljenog govornog signala
y = getaudiodata(x);
audiowrite(filename,y,fs);

end