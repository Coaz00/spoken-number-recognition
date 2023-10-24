clear 
close all
clc

%% Ucitavanje signala
[x, fs] = audioread('jovan.wav');

figure(1) 
plot(1/fs:1/fs:length(x)/fs,x)
title('Vremenski oblik signala')
xlabel('t[s]')


% Prvo cemo segmentirati reci zatim cemo na svakoj proceniti pitch periodu

%% Segmentacija

Wn = [60 3500]/(fs/2);
[B,A] = butter(6,Wn,'bandpass'); 

xf = filter(B,A,x);

wl = 20e-3*fs;

E = zeros(size(xf)); 
Z = zeros(size(xf)); 

for i = wl:length(xf)
    rng = (i - wl + 1):i-1;
    E(i) = sum(xf(rng).^2);
    Z(i) = sum(abs(sign(xf(rng + 1)) - sign(xf(rng))));
end

Z = Z/wl/2;

ITU = 0.1*max(E);
ITL = 0.00001*max(E);

pocetak = [];
kraj = [];

% poredjenje sa vecim pragom ITU
for i = 2:length(E)
    if E(i) > ITU && E(i-1) < ITU
        pocetak = [pocetak i];
    end
    if E(i) < ITU && E(i-1) > ITU
        kraj = [kraj i];
    end
end

% Poredjenje sa nizim pragom
for i = 1:length(pocetak)
    pomeraj = pocetak(i);
    while E(pomeraj)>ITL
        pomeraj = pomeraj - 1;
    end
    pocetak(i) = pomeraj; 
end

for i = 1:length(kraj)
    pomeraj = kraj(i);
    while E(pomeraj)>ITL
        pomeraj = pomeraj + 1;
    end
    kraj(i) = pomeraj; 
end

% uklanjanje duplikata
pocetak = unique(pocetak);
kraj = unique(kraj);

IZCT = 0.2;

for i = 1:length(pocetak)
    if sum((Z(pocetak(i)-25:pocetak(i)-1)>IZCT)>3)
        pocetak(i) = pocetak(i) - find(Z(pocetak(i)-25:pocetak(i)-1)>IZCT,'first');
    end
end

for i = 1:length(kraj)
    if sum((Z(kraj(i)+1:kraj(i)+25)>IZCT)>3)
        kraj(i) = kraj(i) + find(Z(kraj(i)+1:kraj(i)+25)>IZCT,'last');
    end
end

rec=zeros(length(xf),1);
for i=1:length(pocetak)
    rec(pocetak(i):kraj(i),1)=max(E)*ones(kraj(i)-pocetak(i)+1,1);
end

t = 1/fs:1/fs:length(E)/fs;
figure(2)
plot(t,E,t,rec);
xlabel('t[s]')
title('Finalna segmentacija')

for i = 1:length(pocetak)
    reci{i} = xf(pocetak(i):kraj(i));
end 

%% Preslusavanje segmentiranih reci
for i = 1:length(pocetak)
    sound(reci{i},fs)
    pause()
end

%% Metod paralelnog procesiranja
pts = zeros(1,length(pocetak));

for i = 1:length(pocetak)
    x = reci{i};

    % filtriranje pred procenu pitch periode
    Wn = [70 250]/(fs/2); % muski glas -> oko 100Hz 
    [B, A] = butter(6,Wn,'bandpass');
    x = filter(B,A,x);
    
    [m1, m2, m3, m4, m5, m6] = sekvence(x);

    figure(i+2)

    subplot(6,1,1)
    stem(m1(1:1000));

    subplot(6,1,2)
    stem(m2(1:1000));

    subplot(6,1,3)
    stem(m3(1:1000));

    subplot(6,1,4)
    stem(m4(1:1000));

    subplot(6,1,5)
    stem(m5(1:1000));

    subplot(6,1,6)
    stem(m6(1:1000));

    N = length(x);
    pt = procena_periode(fs, N, m1, m2, m3, m4, m5, m6);

    pts(i) = median(pt);

    figure(i+5)
    plot(1./pt)
    xlabel('#procena')
    ylabel('f[Hz]')
    
end

pitch = zeros(1,length(pocetak));

for i = 1:length(pitch)
    pitch(i) = median(1./pts(i));
end

pitch_parallel = median(pitch);

%% Autokorelaciona metoda
pitch = zeros(1,length(pocetak));
for i = 1:length(pocetak)
    x = reci{i};


    N = length(x);

    Wn = [70 200]/(fs/2); % muski glas -> oko 100Hz 
    [B, A] = butter(6,Wn,'bandpass');
    x = filter(B,A,x);

    cl = 0.2*max(x);

    x_clip = zeros(1,N);
    
    % 3-level clipping
    for j = 1:N
        if(x(j)>cl) 
            x_clip(j)=1;
        end
        if(x(j)<-cl)
            x_clip(j)=-1;
        end
    end
    
    w_l = 50e-3*fs;
    Rx = zeros(1,N);
    win = 1:w_l;
    
    % pomerena procena autokorelacione funkcije
    for k = 1: N-w_l
        Rx(k) = 1/(N-w_l)*sum(x_clip(win).*x_clip(win+k));
    end
    
    figure(i+8)
    plot(Rx);
    title('Procena autokorelacione funkcije')
    xlabel('k[odb]','Interpreter','Latex')
    ylabel('$\hat{R}(k)$','Interpreter','Latex')
    
    % pronalazenje razmaka izmedju pikova
    [~,idx] = findpeaks(Rx(1:1000));
    pitch(i) = fs./median(diff(idx));
end

