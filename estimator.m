function Tp = estimator(x, lambda, tau, win, fs)

idx = find(x,1);
A = x(idx);
pocetak = idx + tau; % blanking time
kraj = 0;

for i = pocetak:win
    if (x(i)>=A*exp(-lambda*(i-pocetak))) % prvi stubic koji preseca eksp. f-ju
        kraj = i;
        break;
    end
end

% ako ne presece nikad Tp postvljamo na duzinu prozora
if (kraj == 0)
    Tp = win/fs;
else
    Tp = (kraj - idx) / fs;
end

end