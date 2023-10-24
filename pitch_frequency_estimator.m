function [pt1, pt2, pt3, pt4, pt5, pt6, pt] = procena_periode(fs, N, m1, m2, m3, m4, m5, m6)

win = round(15e-3*fs);
lambda = 120/fs;
tau = 4e-3*fs; % blanking time
Nw = floor(N/(win/2)); % broj prozora, ujedno i duzina pt1...pt6

pt1 = zeros(1,Nw); pt2 = zeros(1,Nw); pt3 = zeros(1,Nw);
pt4 = zeros(1,Nw); pt5 = zeros(1,Nw); pt6 = zeros(1,Nw);
pt = zeros(1,Nw);

i = 2;
for k = 1:win/2:N-win+1
    pt1(i) = estimator(m1(k:k+win-1), lambda, tau, win, fs);
    pt2(i) = estimator(m2(k:k+win-1), lambda, tau, win, fs);
    pt3(i) = estimator(m3(k:k+win-1), lambda, tau, win, fs);
    pt4(i) = estimator(m4(k:k+win-1), lambda, tau, win, fs);
    pt5(i) = estimator(m5(k:k+win-1), lambda, tau, win, fs);
    pt6(i) = estimator(m6(k:k+win-1), lambda, tau, win, fs);
    pt(i) = nanmedian([pt1(i), pt2(i), pt3(i), pt4(i), pt5(i), pt6(i), pt(i-1)]);
    i = i+1;
end